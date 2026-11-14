"""
view3d.py - the 3D tab: a structure sent from the drawing, shown as ball and stick, sticks, wireframe or space-filling.
Drawn with QPainter (depth-sorted shaded spheres and cylinders; space-filling per pixel with a depth buffer, as its
spheres cut into each other), so it needs no graphics driver and no extra library.

    left drag   rotate            click an atom    pick it: 2 atoms distance, 3 angle, 4 dihedral
    right drag  pan               click elsewhere  clear the picks
    wheel       zoom              double-click     the starting view again
"""
import math

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtGui import QColor, QPainter, QPen
from rdkit import Chem

import md_export
import mol3d
import theme

MODES = (('ball', 'Ball and stick'), ('stick', 'Sticks'), ('wire', 'Wireframe'), ('fill', 'Space-filling'))
SURFACE_STYLES = (('solid', 'Solid'), ('translucent', 'Translucent'), ('mesh', 'Mesh'), ('dots', 'Dots'))
SURFACE_COLORS = (('uniform', 'One colour'), ('element', 'Element'), ('charge', 'Partial charge'),
                  ('hydrophobicity', 'Hydrophobicity'))
SURFACE_WHITE = (242, 242, 242)
NEGATIVE, POSITIVE = (230, 60, 50), (50, 90, 230)          # partial charge: red negative, blue positive
HYDROPHOBIC, HYDROPHILIC = (176, 110, 40), (60, 150, 220)
HBOND = '#1FA7B5'
MESH_STEP = 0.75                                            # Å between mesh lines
COLORS = {'H': '#E8E8E8', 'C': '#707070', 'N': '#3050F8', 'O': '#FF0D0D', 'F': '#90E050', 'Cl': '#1FD01F',
          'Br': '#A62929', 'I': '#940094', 'S': '#E6C229', 'P': '#FF8000', 'B': '#FFB5B5', 'Si': '#C8A070',
          'Na': '#AB5CF2', 'K': '#8F40D4', 'Mg': '#8AFF00', 'Ca': '#3DFF00', 'Fe': '#E06633', 'Zn': '#7D80B0',
          'Cu': '#C88033', 'Se': '#FFA100', 'Li': '#CC80FF'}
OTHER_COLOR = '#B080B0'
STICK = 0.13                                                # Å: bond radius
BALL = 0.24                                                 # ball radius as a part of the van der Waals radius
DRAG_PX = 4
_LIGHT = np.array((-0.35, -0.45, 0.82)) / np.linalg.norm((-0.35, -0.45, 0.82))     # from the upper left, screen y down
_HALF = (_LIGHT + (0, 0, 1)) / np.linalg.norm(_LIGHT + (0, 0, 1))


def _shade(c, k):
    """The colour c towards white (k > 0) or black (k < 0)."""
    t, to = (k, 255) if k > 0 else (-k, 0)
    return QColor(int(c.red() + (to - c.red()) * t), int(c.green() + (to - c.green()) * t),
                  int(c.blue() + (to - c.blue()) * t))


def _fog(c, bg, t):
    return QColor(int(c.red() + (bg.red() - c.red()) * t), int(c.green() + (bg.green() - c.green()) * t),
                  int(c.blue() + (bg.blue() - c.blue()) * t))


def start_rotation(xyz):
    """A proper rotation (no mirror image) with the molecule's long axis across the screen and its flat side facing
    the viewer."""
    if len(xyz) < 3:
        return np.eye(3)
    _u, _s, vt = np.linalg.svd(xyz - xyz.mean(axis=0))
    r = vt.copy()
    if np.linalg.det(r) < 0:
        r[2] = -r[2]
    return r


class View3D(QtWidgets.QWidget):
    picks_changed = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.model = None
        self.conf = None
        self.mode = 'ball'
        self.show_h = True
        self.show_labels = False
        self.surface = None                                 # None, 'sas' or 'vdw' (mol3d.SURFACES)
        self.surface_style = 'translucent'
        self.surface_color = 'uniform'
        self.show_hbonds = False
        self.depth_fade = True
        self.perspective = True
        self.show_axes = False
        self.background = None                              # a QColor, or None: the theme's
        self._cache = {}
        self.picks = []                                     # atom indices, in the order clicked
        self.rot = np.eye(3)
        self.zoom = 1.0
        self.pan = QPointF(0, 0)
        self._xyz = None                                    # the conformer, centred on the molecule
        self._radius = 1.0
        self._press = None
        self.setMinimumSize(320, 240)
        self.setMouseTracking(False)

    # ------------------------------------------------------------------------------------------- content
    def set_model(self, model):
        self.model = model
        self.picks = []
        self.conf = model.confs[0][0] if model else None
        self._load()
        self.reset_view()
        self.picks_changed.emit()

    def set_conformer(self, conf_id):
        self.conf = conf_id
        self._load()
        self.update()
        self.picks_changed.emit()                           # the measurement is of this conformer

    def _load(self):
        if not self.model:
            self._xyz = None
            return
        xyz = self.model.positions(self.conf)
        self._xyz = xyz - xyz.mean(axis=0)
        self._cache = {}
        self._radius = max(float(np.linalg.norm(self._xyz, axis=1).max()), 1.0)

    def reset_view(self):
        self.rot = start_rotation(self._xyz) if self._xyz is not None else np.eye(3)
        self.zoom, self.pan = 1.0, QPointF(0, 0)
        self.update()

    def set_mode(self, mode):
        self.mode = mode
        self.update()

    def set_show_h(self, on):
        self.show_h = on
        if not on and self.model:
            self.picks = [a for a in self.picks if self.model.mol.GetAtomWithIdx(a).GetAtomicNum() != 1]
            self.picks_changed.emit()
        self.update()

    def set_show_labels(self, on):
        self.show_labels = on
        self.update()

    def set_option(self, name, value):
        """surface, surface_style, surface_color, show_hbonds, depth_fade, perspective, show_axes or background."""
        setattr(self, name, value)
        self.update()

    def _cached(self, key, make):
        if key not in self._cache:
            self._cache[key] = make()
        return self._cache[key]

    def _surface_colors(self):
        """An (n, 3) array of RGB per atom for the surface."""
        mol = self.model.mol
        n = mol.GetNumAtoms()
        scheme = self.surface_color
        if scheme == 'element':
            return np.array([QColor(COLORS.get(a.GetSymbol(), OTHER_COLOR)).getRgb()[:3] for a in mol.GetAtoms()], float)
        vals = self._cached(('values', scheme), lambda: mol3d.atom_values(self.model, scheme))
        white = np.array(SURFACE_WHITE, float)
        if vals is None:
            return np.tile(white, (n, 1))
        lo, hi = (NEGATIVE, POSITIVE) if scheme == 'charge' else (HYDROPHILIC, HYDROPHOBIC)
        t = np.clip(vals / (0.35 if scheme == 'charge' else 0.5), -1, 1)[:, None]
        return np.where(t < 0, white + (np.array(lo, float) - white) * -t, white + (np.array(hi, float) - white) * t)

    # ------------------------------------------------------------------------------------------- geometry
    def _visible(self):
        return [a.GetIdx() for a in self.model.mol.GetAtoms() if self.show_h or a.GetAtomicNum() != 1]

    def _scale(self, w, h):
        return 0.44 * min(w, h) / (self._radius + 1.7) * self.zoom

    def _atom_radius(self, sym):
        vdw = Chem.GetPeriodicTable().GetRvdw(sym)
        return {'ball': BALL * vdw, 'stick': STICK, 'wire': 0.0, 'fill': vdw}[self.mode]

    def _projector(self, w, h):
        """A function from rotated coordinates (x, y, z towards the viewer; Å) to (screen x, screen y, px per Å)."""
        s = self._scale(w, h)
        camera = 3.5 * (self._radius + 1.7)                 # perspective: the near side a little larger
        cx, cy = w / 2 + self.pan.x(), h / 2 + self.pan.y()

        def project(p):
            f = s * camera / (camera - p[2]) if self.perspective else s
            return cx + p[0] * f, cy - p[1] * f, f
        return project

    def atom_at(self, pos):
        """The front-most visible atom under the point pos, or None."""
        if self._xyz is None:
            return None
        project = self._projector(self.width(), self.height())
        rotated = self._xyz @ self.rot.T
        best, best_z = None, -math.inf
        for i in self._visible():
            x, y, f = project(rotated[i])
            r = max(self._atom_radius(self.model.mol.GetAtomWithIdx(i).GetSymbol()) * f, 8.0)
            if math.hypot(pos.x() - x, pos.y() - y) <= r and rotated[i][2] > best_z:
                best, best_z = i, rotated[i][2]
        return best

    # ------------------------------------------------------------------------------------------- painting
    def paintEvent(self, _e):
        p = QPainter(self)
        self.paint(p, self.width(), self.height(), QColor(theme.T['fig_bg']), QColor(theme.T['fig_fg']))
        p.end()

    def paint(self, p, w, h, bg, fg):
        if self.background is not None:
            bg = QColor(self.background)
            fg = QColor('black' if bg.lightness() > 128 else 'white')
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        self._paint_model(p, w, h, bg, fg)
        if self._xyz is not None:
            if self.surface:
                self._paint_surface(p, w, h, bg)
            self._paint_marks(p, self._xyz @ self.rot.T, self._projector(w, h), fg, set(self._visible()))
            if self.show_axes:
                self._paint_axes(p, h)

    def _paint_model(self, p, w, h, bg, fg):
        p.fillRect(QtCore.QRectF(0, 0, w, h), bg)
        if self._xyz is None:
            p.setPen(fg)
            p.drawText(QtCore.QRectF(0, 0, w, h), Qt.AlignmentFlag.AlignCenter | Qt.TextFlag.TextWordWrap,
                       'Select a molecule in the Drawing tab, then Structure > Send to 3D (Ctrl+3).')
            return
        mol = self.model.mol
        project = self._projector(w, h)
        rotated = self._xyz @ self.rot.T
        visible = set(self._visible())
        if self.surface and self.surface_style == 'solid':
            return                                          # the surface hides it all
        if self.mode == 'fill':
            self._paint_fill(p, w, h, rotated, project, bg, visible)
            return
        zmax = self._radius or 1.0
        fog = (lambda z: 0.45 * (zmax - z) / (2 * zmax)) if self.depth_fade else (lambda z: 0.0)
        items = []                                          # (depth, kind, data), far ones are painted first
        for i in visible:
            items.append((rotated[i][2], 'atom', i))
        for b in mol.GetBonds():
            a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if a in visible and c in visible:
                mid = (rotated[a] + rotated[c]) / 2
                items.append(((rotated[a][2] + mid[2]) / 2, 'half', (a, c, mid, b)))
                items.append(((rotated[c][2] + mid[2]) / 2, 'half', (c, a, mid, b)))
        items.sort(key=lambda it: it[0])
        for _z, kind, data in items:
            if kind == 'atom':
                self._paint_atom(p, data, rotated, project, bg, fog)
            else:
                self._paint_half(p, data, rotated, project, bg, fog)

    def _paint_fill(self, p, w, h, rotated, project, bg, visible):
        """Space-filling: every sphere per pixel, the nearest surface wins."""
        mol = self.model.mol
        atoms = sorted(visible)
        radii = np.array([self._atom_radius(mol.GetAtomWithIdx(i).GetSymbol()) for i in atoms])
        colors = np.array([QColor(COLORS.get(mol.GetAtomWithIdx(i).GetSymbol(), OTHER_COLOR)).getRgb()[:3]
                           for i in atoms], float)
        rgb, _mask, _where, ss = self._spheres(w, h, rotated[atoms], radii, colors, project, bg)
        self._draw_buffer(p, w, h, rgb, None, ss)

    def _spheres(self, w, h, centres, radii, colors, project, bg, coords=False):
        """Shaded spheres per pixel, the nearest surface winning, as (rgb, covered mask, rotated xyz per pixel or None,
        supersampling). Twice the resolution for smooth edges, except while the mouse is turning the model."""
        ss = 1 if self._press else 2
        W, H = int(w * ss), int(h * ss)
        depth = np.full((H, W), -np.inf)
        rgb = np.empty((H, W, 3), dtype=np.uint8)
        rgb[:] = (bg.red(), bg.green(), bg.blue())
        where = np.zeros((H, W, 3)) if coords else None
        for centre, ra, base in zip(centres, radii, colors):
            x, y, f = project(centre)
            r = ra * f * ss
            x0, x1 = max(int(x * ss - r), 0), min(int(x * ss + r) + 2, W)
            y0, y1 = max(int(y * ss - r), 0), min(int(y * ss + r) + 2, H)
            if x0 >= x1 or y0 >= y1:
                continue
            dx = ((np.arange(x0, x1) + 0.5 - x * ss) / r)[None, :]
            dy = ((np.arange(y0, y1) + 0.5 - y * ss) / r)[:, None]
            rr = dx * dx + dy * dy
            nz = np.sqrt(np.clip(1.0 - rr, 0.0, None))
            z = centre[2] + nz * ra
            near = (rr <= 1.0) & (z > depth[y0:y1, x0:x1])
            if not near.any():
                continue
            depth[y0:y1, x0:x1][near] = z[near]
            diffuse = np.clip(dx * _LIGHT[0] + dy * _LIGHT[1] + nz * _LIGHT[2], 0.0, 1.0)
            shine = np.clip(dx * _HALF[0] + dy * _HALF[1] + nz * _HALF[2], 0.0, 1.0) ** 40 * 0.55
            shade = (0.3 + 0.7 * diffuse)[..., None] * base + shine[..., None] * 255.0
            rgb[y0:y1, x0:x1][near] = np.clip(shade, 0, 255)[near].astype(np.uint8)
            if coords:
                here = np.stack(np.broadcast_arrays(centre[0] + dx * ra, centre[1] - dy * ra, z), axis=-1)
                where[y0:y1, x0:x1][near] = here[near]
        return rgb, depth > -np.inf, where, ss

    def _draw_buffer(self, p, w, h, rgb, alpha, ss):
        """Draws an RGB buffer, with alpha (0-255 per pixel) when given, scaled back from its supersampling."""
        H, W = rgb.shape[:2]
        if alpha is None:
            data = np.ascontiguousarray(rgb)
            img = QtGui.QImage(data.data, W, H, 3 * W, QtGui.QImage.Format.Format_RGB888).copy()
        else:
            data = np.ascontiguousarray(np.dstack((rgb, alpha.astype(np.uint8))))
            img = QtGui.QImage(data.data, W, H, 4 * W, QtGui.QImage.Format.Format_RGBA8888).copy()
        p.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        p.drawImage(QtCore.QRectF(0, 0, w, h), img)

    def _paint_surface(self, p, w, h, bg):
        """The molecular surface over the model: solid, translucent, a mesh of lines fixed to the molecule, or dots.
        Every atom counts, shown or not."""
        project = self._projector(w, h)
        rotated = self._xyz @ self.rot.T
        radii = mol3d.surface_radii(self.model, self.surface)
        colors = self._surface_colors()
        if self.surface_style == 'dots':
            pts, owners, _area = self._cached(('dots', self.surface), lambda: mol3d.surface_dots(self._xyz, radii, 160))
            pr = pts @ self.rot.T
            zmax = self._radius + radii.max()
            p.setPen(Qt.PenStyle.NoPen)
            for k in np.argsort(pr[:, 2]):
                x, y, _f = project(pr[k])
                c = QColor(*[int(v) for v in colors[owners[k]]])
                if self.depth_fade:
                    c = _fog(c, bg, 0.5 * (zmax - pr[k][2]) / (2 * zmax))
                p.setBrush(c)
                p.drawEllipse(QPointF(x, y), 1.3, 1.3)
            return
        mesh = self.surface_style == 'mesh'
        rgb, mask, where, ss = self._spheres(w, h, rotated, radii, colors, project, bg, coords=mesh)
        if self.surface_style == 'solid':
            self._draw_buffer(p, w, h, rgb, None, ss)
            return
        if self.surface_style == 'translucent':
            self._draw_buffer(p, w, h, rgb, mask * 120, ss)
            return
        obj = (where.reshape(-1, 3) @ self.rot).reshape(where.shape) / MESH_STEP      # back to the molecule's own axes
        frac = np.abs(obj - np.round(obj))
        line = mask & (frac.min(axis=-1) < 0.1)
        self._draw_buffer(p, w, h, rgb, line * 230, ss)

    def _paint_axes(self, p, h):
        """The molecule's x, y and z axes, small, in the lower left corner."""
        ox, oy, length = 34.0, h - 34.0, 22.0
        font = p.font(); font.setPointSizeF(8); p.setFont(font)
        for k, (name, color) in enumerate((('x', '#D03030'), ('y', '#20A040'), ('z', '#3050E0'))):
            v = self.rot[:, k]
            ex, ey = ox + v[0] * length, oy - v[1] * length
            p.setPen(QPen(QColor(color), 2))
            p.drawLine(QPointF(ox, oy), QPointF(ex, ey))
            p.drawText(QPointF(ex + (3 if v[0] >= 0 else -9), ey + (-3 if v[1] >= 0 else 11)), name)

    def _paint_atom(self, p, i, rotated, project, bg, fog):
        sym = self.model.mol.GetAtomWithIdx(i).GetSymbol()
        x, y, f = project(rotated[i])
        c = _fog(QColor(COLORS.get(sym, OTHER_COLOR)), bg, fog(rotated[i][2]))
        r = self._atom_radius(sym) * f
        p.setPen(Qt.PenStyle.NoPen)
        if self.mode == 'wire':                             # a bond-less atom must still show; others get a dot
            r = max(1.5, 0.05 * f)
            p.setBrush(c)
        else:
            g = QtGui.QRadialGradient(QPointF(x, y), r, QPointF(x - r * 0.35, y - r * 0.4))
            g.setColorAt(0.0, _shade(c, 0.65)); g.setColorAt(0.45, c); g.setColorAt(1.0, _shade(c, -0.55))
            p.setBrush(QtGui.QBrush(g))
        p.drawEllipse(QPointF(x, y), r, r)

    def _paint_half(self, p, data, rotated, project, bg, fog):
        """Half a bond, from atom a to the middle, in a's colour; a multiple bond as parallel lines."""
        a, c, mid, bond = data
        mol = self.model.mol
        sym = mol.GetAtomWithIdx(a).GetSymbol()
        col = QColor(COLORS.get(sym, OTHER_COLOR))
        if self.mode == 'wire' and sym == 'H' and bg.lightness() > 128:
            col = QColor('#A0A0A0')                         # white hydrogen lines vanish on a white page
        col = _fog(col, bg, fog((rotated[a][2] + mid[2]) / 2))
        start = rotated[a]
        d = rotated[c] - rotated[a]
        length = float(np.linalg.norm(d))
        if self.mode == 'ball' and length:                  # from the ball's surface: never drawn over its own ball
            start = rotated[a] + d / length * min(self._atom_radius(sym) * 0.8, length / 2)
        x1, y1, f1 = project(start)
        x2, y2, _f2 = project(mid)
        dx, dy = x2 - x1, y2 - y1
        n = math.hypot(dx, dy)
        if n < 0.5:
            return
        ox, oy = -dy / n, dx / n                            # across the line, on the screen
        wire = self.mode == 'wire'
        order = bond.GetBondType()
        lines = {Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3, Chem.BondType.AROMATIC: 2}.get(order, 1)
        if self.mode == 'stick':
            lines = 1                                       # Sticks: only the skeleton, as is usual
        radius = (0.5 / f1 if wire else STICK if lines == 1 else STICK * 0.62) * f1
        gap = (0.09 * f1 if wire else radius * 1.35)
        for k in range(lines):
            off = (k - (lines - 1) / 2) * 2 * gap
            ax, ay, bx, by = x1 + ox * off, y1 + oy * off, x2 + ox * off, y2 + oy * off
            dashed = order == Chem.BondType.AROMATIC and k == 1
            if wire:
                pen = QPen(col, 1.5)
            else:
                g = QtGui.QLinearGradient(QPointF(ax + ox * radius, ay + oy * radius),
                                          QPointF(ax - ox * radius, ay - oy * radius))
                bright, dark = (0.0, 1.0) if oy < 0 else (1.0, 0.0)     # lit from the top
                g.setColorAt(bright, _shade(col, 0.45)); g.setColorAt(dark, _shade(col, -0.5)); g.setColorAt(0.5, col)
                pen = QPen(QtGui.QBrush(g), 2 * radius)
            pen.setCapStyle(Qt.PenCapStyle.FlatCap)
            if dashed:
                pen.setStyle(Qt.PenStyle.DashLine)
            p.setPen(pen)
            p.drawLine(QPointF(ax, ay), QPointF(bx, by))

    def _paint_marks(self, p, rotated, project, fg, visible):
        mol = self.model.mol
        font = p.font(); font.setPointSizeF(9); p.setFont(font)
        if self.show_hbonds:
            p.setPen(QPen(QColor(HBOND), 1.6, Qt.PenStyle.DashLine))
            for hyd, acc in self._cached(('hbonds',), lambda: mol3d.hydrogen_bonds(self.model, self.conf)):
                start = hyd if hyd in visible else mol.GetAtomWithIdx(hyd).GetNeighbors()[0].GetIdx()
                x1, y1, _f = project(rotated[start])
                x2, y2, _g = project(rotated[acc])
                p.drawLine(QPointF(x1, y1), QPointF(x2, y2))
        if self.show_labels:
            p.setPen(fg)
            for i in sorted(visible, key=lambda k: rotated[k][2]):
                x, y, f = project(rotated[i])
                r = self._atom_radius(mol.GetAtomWithIdx(i).GetSymbol()) * f
                p.drawText(QPointF(x + r * 0.7 + 2, y - r * 0.7 - 2), mol3d.atom_label(mol, i))
        accent = QColor(theme.current_accent)
        points = [project(rotated[i]) for i in self.picks]
        p.setPen(QPen(accent, 1.4, Qt.PenStyle.DashLine))
        for (x1, y1, _f), (x2, y2, _g) in zip(points, points[1:]):
            p.drawLine(QPointF(x1, y1), QPointF(x2, y2))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.setPen(QPen(accent, 2))
        for n, (i, (x, y, f)) in enumerate(zip(self.picks, points)):
            r = max(self._atom_radius(mol.GetAtomWithIdx(i).GetSymbol()) * f, 6.0) + 3
            p.drawEllipse(QPointF(x, y), r, r)
            p.drawText(QPointF(x + r + 2, y - r), str(n + 1))

    def image(self, scale=2):
        """The view as a picture on white, at scale times the size on the screen."""
        w, h = self.width(), self.height()
        img = QtGui.QImage(w * scale, h * scale, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        p = QPainter(img)
        p.scale(scale, scale)
        self.paint(p, w, h, QColor('white'), QColor('black'))     # a chosen background wins over the white
        p.end()
        return img

    # ------------------------------------------------------------------------------------------- mouse
    def mousePressEvent(self, e):
        self._press = (e.button(), e.position(), e.position(), False)

    def mouseMoveEvent(self, e):
        if not self._press:
            return
        btn, start, last, moved = self._press
        pos = e.position()
        moved = moved or (pos - start).manhattanLength() > DRAG_PX
        if moved:
            dx, dy = pos.x() - last.x(), pos.y() - last.y()
            if btn == Qt.MouseButton.LeftButton:
                self.rotate(dx * 0.01, dy * 0.01)
            else:
                self.pan += QPointF(dx, dy)
                self.update()
        self._press = (btn, start, pos, moved)

    def rotate(self, about_y, about_x):
        """Turns the molecule around the screen's axes (radians)."""
        cy, sy, cx, sx = math.cos(about_y), math.sin(about_y), math.cos(about_x), math.sin(about_x)
        ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
        rx = np.array([[1, 0, 0], [0, cx, sx], [0, -sx, cx]])
        self.rot = rx @ ry @ self.rot
        self.update()

    def mouseReleaseEvent(self, e):
        if not self._press:
            return
        btn, _start, _last, moved = self._press
        self._press = None
        if btn == Qt.MouseButton.LeftButton and not moved and self.model:
            self.pick(self.atom_at(e.position()))

    def pick(self, atom):
        """Adds an atom to the measurement (a fifth starts again; clicking the last one or nothing removes picks)."""
        if atom is None:
            self.picks = []
        elif atom in self.picks:
            self.picks = self.picks[:self.picks.index(atom)]
        else:
            self.picks = (self.picks if len(self.picks) < 4 else []) + [atom]
        self.update()
        self.picks_changed.emit()

    def mouseDoubleClickEvent(self, _e):
        self.reset_view()

    def wheelEvent(self, e):
        self.zoom = min(max(self.zoom * 1.15 ** (e.angleDelta().y() / 120), 0.1), 30.0)
        self.update()


class _BuildWorker(QtCore.QThread):
    done = QtCore.pyqtSignal(object, str)

    def __init__(self, mol, name, n_confs, force_field, parent):
        super().__init__(parent)
        self.mol, self.name, self.n_confs, self.force_field = mol, name, n_confs, force_field

    def run(self):
        try:
            self.done.emit(mol3d.build(self.mol, self.name, self.n_confs, force_field=self.force_field), '')
        except ValueError as e:
            self.done.emit(None, str(e))
        except Exception as e:                              # shown in the tab, also logged
            import logging
            logging.getLogger('kekulex').warning('3D model failed', exc_info=True)
            self.done.emit(None, str(e))


class Panel3D(QtWidgets.QWidget):
    """The 3D tab: the view and, beside it, what is shown and the conformers."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.worker = None
        self.view = View3D()
        self.view.picks_changed.connect(self._show_measure)
        side = QtWidgets.QWidget()
        col = QtWidgets.QVBoxLayout(side)
        self.title = QtWidgets.QLabel('Nothing sent yet.')
        self.title.setWordWrap(True)
        self.title.setStyleSheet('font-weight: bold')
        col.addWidget(self.title)
        self.notes = QtWidgets.QLabel()
        self.notes.setWordWrap(True)
        col.addWidget(self.notes)
        col.addWidget(QtWidgets.QLabel('Display:'))
        self.mode = QtWidgets.QComboBox()
        for key, text in MODES:
            self.mode.addItem(text, key)
        self.mode.currentIndexChanged.connect(lambda _i: self.view.set_mode(self.mode.currentData()))
        col.addWidget(self.mode)
        self.hydrogens = QtWidgets.QCheckBox('Hydrogens')
        self.hydrogens.setChecked(True)
        self.hydrogens.toggled.connect(self.view.set_show_h)
        col.addWidget(self.hydrogens)
        self.labels = QtWidgets.QCheckBox('Atom labels')
        self.labels.toggled.connect(self.view.set_show_labels)
        col.addWidget(self.labels)
        grid = QtWidgets.QGridLayout()
        for n, (attr, text, on, tip) in enumerate((
                ('show_hbonds', 'H-bonds', False, 'Hydrogen bonds inside the molecule, dashed'),
                ('depth_fade', 'Depth fading', True, 'Atoms further away fade into the background'),
                ('perspective', 'Perspective', True, 'Nearer parts drawn larger'),
                ('show_axes', 'Axes', False, "The molecule's x, y and z axes in the corner"))):
            box = QtWidgets.QCheckBox(text)
            box.setChecked(on)
            box.setToolTip(tip)
            box.toggled.connect(lambda v, a=attr: self.view.set_option(a, v))
            grid.addWidget(box, n // 2, n % 2)
        col.addLayout(grid)
        self.background_button = QtWidgets.QPushButton('Background colour…')
        self.background_button.setToolTip('Right-click: back to the theme\'s background')
        self.background_button.clicked.connect(self.choose_background)
        self.background_button.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.background_button.customContextMenuRequested.connect(lambda _p: self.view.set_option('background', None))
        col.addWidget(self.background_button)
        surf = QtWidgets.QFormLayout()
        self.surface = QtWidgets.QComboBox()
        self.surface.addItem('None', None)
        for key, text in mol3d.SURFACES:
            self.surface.addItem(text, key)
        self.surface.setToolTip(f'Solvent accessible: where the centre of a water molecule ({mol3d.PROBE} Å) can go. '
                                'Van der Waals: the outside of the atoms\' spheres.')
        self.surface_style = QtWidgets.QComboBox()
        self.surface_color = QtWidgets.QComboBox()
        for box, items, attr in ((self.surface, None, 'surface'), (self.surface_style, SURFACE_STYLES, 'surface_style'),
                                 (self.surface_color, SURFACE_COLORS, 'surface_color')):
            for key, text in items or ():
                box.addItem(text, key)
            box.currentIndexChanged.connect(lambda _i, b=box, a=attr: self.view.set_option(a, b.currentData()))
        self.surface_style.setCurrentIndex(1)
        self.surface_color.setToolTip('Partial charge (Gasteiger): red negative, blue positive. Hydrophobicity '
                                      '(Crippen logP per atom): brown fat-loving, blue water-loving.')
        surf.addRow('Surface:', self.surface)
        surf.addRow('Style:', self.surface_style)
        surf.addRow('Colour:', self.surface_color)
        col.addLayout(surf)
        form = QtWidgets.QFormLayout()
        self.force_field = QtWidgets.QComboBox()
        self.force_field.addItems(mol3d.FORCE_FIELDS)
        self.force_field.setToolTip('The force field for new conformers. MMFF94s suits crystal-like geometries; '
                                    'UFF covers more elements but is rougher.')
        form.addRow('Force field:', self.force_field)
        self.n_confs = QtWidgets.QSpinBox()
        self.n_confs.setRange(1, 200)
        self.n_confs.setValue(20)
        form.addRow('Search:', self.n_confs)
        col.addLayout(form)
        self.search_button = QtWidgets.QPushButton('Search conformers again')
        self.search_button.setToolTip('New conformers of this molecule with the force field and number above')
        self.search_button.clicked.connect(self.search_again)
        col.addWidget(self.search_button)
        col.addWidget(QtWidgets.QLabel('Conformers (energy, kcal/mol):'))
        self.confs = QtWidgets.QListWidget()
        self.confs.currentRowChanged.connect(self._conformer_chosen)
        col.addWidget(self.confs, 1)
        self.minimize_button = QtWidgets.QPushButton('Minimize this one')
        self.minimize_button.clicked.connect(self.minimize)
        self.minimize_all_button = QtWidgets.QPushButton('Minimize all')
        self.minimize_all_button.clicked.connect(self.minimize_all)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.minimize_button, 2)
        row.addWidget(self.minimize_all_button, 1)
        col.addLayout(row)
        self.energy = QtWidgets.QLabel()
        self.energy.setWordWrap(True)
        col.addWidget(self.energy)
        col.addWidget(QtWidgets.QLabel('Measure:'))
        self.measure = QtWidgets.QLabel()
        self.measure.setWordWrap(True)
        self.measure.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        col.addWidget(self.measure)
        col.addWidget(QtWidgets.QLabel('Properties of this conformer:'))
        self.properties = QtWidgets.QLabel()
        self.properties.setWordWrap(True)
        self.properties.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        col.addWidget(self.properties)
        row = QtWidgets.QHBoxLayout()
        for text, slot in (('Reset view', self.view.reset_view), ('Copy', self.copy_picture), ('PNG…', self.export_png)):
            b = QtWidgets.QPushButton(text)
            b.clicked.connect(lambda _c=False, f=slot: f())
            row.addWidget(b)
        col.addLayout(row)
        self.export_button = QtWidgets.QPushButton('Export for MD programs…')
        self.export_button.setToolTip(md_export.NOTE)
        self.export_button.clicked.connect(self.export_md)
        col.addWidget(self.export_button)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(side)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setFixedWidth(272)
        self.confs.setMinimumHeight(110)
        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.view, 1)
        lay.addWidget(scroll)
        self._show_measure()
        self._enable(False)

    def _enable(self, on):
        for w in (self.minimize_button, self.minimize_all_button):
            w.setEnabled(on and bool(self.view.model and self.view.model.method))
        for w in (self.search_button, self.export_button):
            w.setEnabled(on)

    # ------------------------------------------------------------------------------------------- building
    def busy(self):
        return self.worker is not None and self.worker.isRunning()

    def search_again(self):
        model = self.view.model
        if model:
            self.send(model.source, model.name.split('  (')[0])

    def send(self, mol, name):
        """Starts building the model of mol; False if one is still being built."""
        if self.busy():
            return False
        self.title.setText(f'Building {name} …')
        self.notes.setText('')
        self._enable(False)
        self.worker = _BuildWorker(mol, name, self.n_confs.value(), self.force_field.currentText(), self)
        self.worker.done.connect(self._built)
        self.worker.start()
        return True

    def _built(self, model, error):
        if model is None:
            self.title.setText('No 3D model')
            self.notes.setText(error)
            return
        self.title.setText(model.name + (f'  ({model.method})' if model.method else ''))
        self.notes.setText('\n'.join(model.notes))
        self.view.set_model(model)
        self._fill_confs(0)
        self._show_energy()
        self._enable(True)

    def _fill_confs(self, row):
        confs = self.view.model.confs
        lowest = confs[0][1]
        self.confs.blockSignals(True)
        self.confs.clear()
        for n, (_c, e) in enumerate(confs):
            self.confs.addItem(f'{n + 1}:  not minimised' if e is None else f'{n + 1}:  {e:.2f}   (+{e - lowest:.2f})')
        self.confs.setCurrentRow(row)
        self.confs.blockSignals(False)

    def _conformer_chosen(self, row):
        model = self.view.model
        if model and 0 <= row < len(model.confs):
            self.view.set_conformer(model.confs[row][0])
            self._show_energy()

    def minimize(self):
        model = self.view.model
        if not model:
            return
        conf = self.view.conf
        before = dict(model.confs)[conf]
        e, converged = mol3d.minimize(model, conf)
        self.view.set_conformer(conf)
        self._fill_confs(next(n for n, (c, _x) in enumerate(model.confs) if c == conf))
        self._show_energy(f'{before:.2f} → {e:.2f} kcal/mol' + ('' if converged else ', not converged: press it again'))

    def minimize_all(self):
        model = self.view.model
        if not model:
            return
        conf = self.view.conf
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            left = mol3d.minimize_all(model)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self.view.set_conformer(conf)
        self._fill_confs(next(n for n, (c, _x) in enumerate(model.confs) if c == conf))
        self._show_energy('All minimized' + (f', {left} not converged' if left else ''))

    def _show_energy(self, first=''):
        """What the last minimization did and how far this conformer is from a minimum; and its properties."""
        v = self.view
        g = mol3d.rms_gradient(v.model, v.conf) if v.model else None
        self.energy.setText('\n'.join(t for t in (first, '' if g is None else f'RMS force: {g:.4f} kcal/mol/Å') if t))
        self.properties.setText('\n'.join(f'{name}: {text}' for name, text in mol3d.shape_properties(v.model, v.conf))
                                if v.model else '')

    def choose_background(self):
        c = QtWidgets.QColorDialog.getColor(self.view.background or QColor(theme.T['fig_bg']), self, 'Background colour')
        if c.isValid():
            self.view.set_option('background', c)

    # ------------------------------------------------------------------------------------------- output
    def _show_measure(self):
        v = self.view
        if not v.model or not v.picks:
            self.measure.setText('Click 2 atoms for a distance, 3 for an angle, 4 for a dihedral.')
            return
        names = ' – '.join(mol3d.atom_label(v.model.mol, a) for a in v.picks)
        m = mol3d.measure(v.model, v.conf, v.picks)
        self.measure.setText(f'{names}\n{m[0]}: {m[1]:.3f} {m[2]}' if m and m[2] == 'Å'
                             else f'{names}\n{m[0]}: {m[1]:.1f} {m[2]}' if m else names)

    def copy_picture(self):
        QtWidgets.QApplication.clipboard().setImage(self.view.image())

    def export_png(self):
        path, _f = QtWidgets.QFileDialog.getSaveFileName(self, 'Save picture', 'structure_3d.png', 'PNG picture (*.png)')
        if path:
            self.view.image().save(path)

    def export_md(self):
        v = self.view
        if not v.model:
            return
        filters = ';;'.join(f for f, _k in md_export.FORMATS)
        path, chosen = QtWidgets.QFileDialog.getSaveFileName(self, 'Export for MD programs', 'structure.mol2', filters,
                                                             md_export.FORMATS[1][0])
        if not path:
            return
        fmt = dict(md_export.FORMATS)[chosen]
        try:
            md_export.write(v.model, v.conf, fmt, path)
        except OSError as e:
            QtWidgets.QMessageBox.warning(self, 'Export', f'Could not write the file: {e}')
            return
        QtWidgets.QMessageBox.information(self, 'Export', md_export.NOTE)

    def wait(self):
        """Lets a running build end (the program is closing: it cannot be stopped half-way)."""
        if self.worker is not None:
            self.worker.wait()
