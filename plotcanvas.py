"""
plotcanvas.py - a matplotlib plot inside Qt with the viewer's mouse handling, the same on every plot:

    drag (left button)      zoom to the box               click (no drag)   on_click(x)  (e.g. move the lambda cursor)
    wheel                   zoom x around the mouse       Ctrl+wheel        zoom y
    right drag              pan                           double-click / F  whole range again, y automatic

The figure holds one or more panels stacked top to bottom on a shared x axis (layout(n)), like the PDA and MS
chromatograms of LC viewer: x zooms and pans move every panel, y only the panel under the mouse. A panel may carry a
right-hand axis (twin(ax)), which is styled with it and lets the mouse through to the panel.
"""
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from PyQt6 import QtCore

import theme

DRAG_PX = 5                                                 # less movement than this is a click


class PlotCanvas(FigureCanvasQTAgg):
    x_only = False                                          # True: drag zoom and pan leave y alone (Ctrl+wheel still zooms it)

    def __init__(self, parent=None):
        self.fig = Figure(figsize=(7, 4), layout='constrained')
        super().__init__(self.fig)
        self.setParent(parent)
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.ClickFocus)
        self.panels = []                                    # main axes, top to bottom
        self.twins = {}                                     # panel -> its right-hand axis
        self.layout(1)
        self.view = None                                    # {'x': xlim, 'y': {panel index: ylim}} or None = automatic
        self.on_click = None                                # callables set by the owner: on_click(x)
        self.on_click_xy = None                             # on_click_xy(panel index, x, y), called instead when set
        self.on_hover = None                                # on_hover(x, y)  (None, None outside)
        self.on_hover_xy = None                             # on_hover_xy(panel index or None, x, y)
        self.on_key = None                                  # on_key(key name): keys other than F
        self.on_reset = None
        self.on_view_changed = None                         # after a mouse zoom / pan (the axis sliders follow)
        self._press = None                                  # (button, panel index, x px, y px, xdata, ydata, xlim, ylim)
        self._box = None
        for ev, fn in (('button_press_event', self._down), ('motion_notify_event', self._move),
                       ('button_release_event', self._up), ('scroll_event', self._wheel),
                       ('key_press_event', self._key)):
            self.mpl_connect(ev, fn)

    @property
    def ax(self):
        return self.panels[0]

    @property
    def ax2(self):
        """The first panel's right-hand axis, or None."""
        return self.twins.get(self.panels[0])

    @ax2.setter
    def ax2(self, value):
        if value is None:
            self.twins.pop(self.panels[0], None)
        else:
            self.twins[self.panels[0]] = value

    def layout(self, n, heights=None):
        """Clear the figure and make n panels on a shared x axis (top to bottom); returns them. A different number of
        panels than before starts with automatic limits (the zoom is per panel)."""
        if n != len(self.panels):
            self.view = None
        self.fig.clear()
        self.twins = {}
        gs = self.fig.add_gridspec(n, 1, height_ratios=heights or [1] * n)
        self.panels = []
        for i in range(n):
            ax = self.fig.add_subplot(gs[i, 0], sharex=self.panels[0] if self.panels else None)
            if i < n - 1:
                ax.tick_params(labelbottom=False)
            self.panels.append(ax)
        return self.panels

    def twin(self, ax):
        """The right-hand axis of panel ax (made on first use)."""
        if ax not in self.twins:
            self.twins[ax] = ax.twinx()
        return self.twins[ax]

    # -------------------------------------------------------------------------------------------- look
    def style(self):
        """Theme colours on the figure (call after drawing into the axes)."""
        T = theme.T
        self.fig.set_facecolor(T['fig_bg'])
        for ax in self.panels:
            ax.set_facecolor(T['fig_bg'])
            ax.grid(True, color=T['fig_grid'], lw=0.6)
            ax.set_axisbelow(True)
        for ax in self.panels + list(self.twins.values()):
            for s in ax.spines.values():
                s.set_color(T['fig_fg'])
            ax.tick_params(colors=T['fig_fg'], which='both')
            ax.xaxis.label.set_color(T['fig_fg']); ax.yaxis.label.set_color(T['fig_fg']); ax.title.set_color(T['fig_fg'])
            leg = ax.get_legend()
            if leg is not None:
                leg.get_frame().set_facecolor(T['fig_bg']); leg.get_frame().set_edgecolor(T['fig_grid'])
                for t in leg.get_texts():
                    t.set_color(T['fig_fg'])

    def apply_view(self):
        """After a redraw: keep the user's zoom, or let matplotlib fit the data."""
        if not self.view:
            return
        if self.view.get('x'):
            self.ax.set_xlim(self.view['x'])                # shared: every panel follows
        for i, yl in self.view.get('y', {}).items():
            if i < len(self.panels):
                self.panels[i].set_ylim(yl)

    def reset_view(self):
        self.view = None
        if self.on_reset:
            self.on_reset()

    def _remember(self, i=None):
        v = self.view or {'x': None, 'y': {}}
        v['x'] = self.ax.get_xlim()
        if i is not None:
            v['y'][i] = self.panels[i].get_ylim()
        self.view = v
        if self.on_view_changed:
            self.on_view_changed()

    def set_limits(self, x=None, y=None, panel=0):
        """Limits from outside (the axis sliders): x for every panel, y for one; kept like a mouse zoom, drawn at once
        without redrawing the data."""
        v = self.view or {'x': None, 'y': {}}
        if x is not None:
            self.ax.set_xlim(x); v['x'] = tuple(x)
        if y is not None and panel < len(self.panels):
            self.panels[panel].set_ylim(y); v['y'][panel] = tuple(y)
        self.view = v
        self.draw_idle()

    def auto_y(self):
        """Forget the y zoom of every panel (x stays)."""
        if self.view:
            self.view['y'] = {}

    # -------------------------------------------------------------------------------------------- mouse
    def _in(self, e):
        """The panel index under the event (None outside), with e.xdata / e.ydata in that panel's units (also over its
        right-hand axis, which lies on top)."""
        ax = e.inaxes
        for p, t in self.twins.items():
            if ax is t:
                ax = p
        if ax is None or ax not in self.panels:
            e.xdata = e.ydata = None
            return None
        if ax is not e.inaxes:
            e.xdata, e.ydata = ax.transData.inverted().transform((e.x, e.y))
        return self.panels.index(ax)

    def _down(self, e):
        i = self._in(e)
        if i is None:
            return
        if e.dblclick:
            self.reset_view(); return
        ax = self.panels[i]
        self._press = (e.button, i, e.x, e.y, e.xdata, e.ydata, ax.get_xlim(), ax.get_ylim())

    def _move(self, e):
        i = self._in(e)
        if self.on_hover:
            self.on_hover(e.xdata, e.ydata)
        if self.on_hover_xy:
            self.on_hover_xy(i, e.xdata, e.ydata)
        if not self._press:
            return
        btn, pi, px, py, x0, y0, xl, yl = self._press
        ax = self.panels[pi]
        if btn == 1 and i == pi and (abs(e.x - px) > DRAG_PX or abs(e.y - py) > DRAG_PX):
            if self._box is None:
                self._box = Rectangle((x0, y0), 0, 0, fill=True, alpha=0.12, lw=1, ls='--',
                                      ec=theme.T['fig_fg'], fc=theme.current_accent, zorder=10)
                ax.add_patch(self._box)
            if self.x_only:
                self._box.set_bounds(min(x0, e.xdata), min(yl), abs(e.xdata - x0), abs(yl[1] - yl[0]))
            else:
                self._box.set_bounds(min(x0, e.xdata), min(y0, e.ydata), abs(e.xdata - x0), abs(e.ydata - y0))
            self.draw_idle()
        elif btn == 3:                                      # pan: move the limits by the mouse travel (in data units)
            inv = ax.transData.inverted()
            (ax0, ay0), (ax1, ay1) = inv.transform((px, py)), inv.transform((e.x, e.y))
            dx, dy = ax1 - ax0, ay1 - ay0
            ax.set_xlim(xl[0] - dx, xl[1] - dx)
            if not self.x_only:
                ax.set_ylim(yl[0] - dy, yl[1] - dy)
            self._remember(pi)
            self.draw_idle()

    def _up(self, e):
        if not self._press:
            return
        btn, pi, px, py, x0, y0, *_ = self._press
        self._press = None
        if self._box is not None:
            x, y, w, h = self._box.get_x(), self._box.get_y(), self._box.get_width(), self._box.get_height()
            self._box.remove(); self._box = None
            if w > 0 and h > 0:
                ax = self.panels[pi]
                ax.set_xlim((x, x + w) if ax.get_xlim()[0] < ax.get_xlim()[1] else (x + w, x))
                if not self.x_only:
                    ax.set_ylim(y, y + h)
                self._remember(pi)
            self.draw_idle()
        elif btn == 1 and abs(e.x - px) <= DRAG_PX and abs(e.y - py) <= DRAG_PX and x0 is not None:
            if self.on_click_xy:
                self.on_click_xy(pi, x0, y0)
            elif self.on_click:
                self.on_click(x0)

    def _wheel(self, e):
        i = self._in(e)
        if i is None:
            return
        ax = self.panels[i]
        f = 0.8 if e.button == 'up' else 1.25
        ctrl = 'control' in (e.key or '') or 'ctrl' in (e.key or '')
        if ctrl:
            lo, hi = ax.get_ylim(); c = e.ydata
            ax.set_ylim(c - (c - lo) * f, c + (hi - c) * f)
            self._remember(i)
        else:
            lo, hi = ax.get_xlim(); c = e.xdata
            ax.set_xlim(c - (c - lo) * f, c + (hi - c) * f)
            self._remember()
        self.draw_idle()

    def _key(self, e):
        if e.key in ('f', 'F'):
            self.reset_view()
        elif self.on_key and e.key:
            self.on_key(e.key)
