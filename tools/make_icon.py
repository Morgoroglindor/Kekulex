"""Draws structureeditor.ico (the program / installer icon). Several designs are kept so they can be compared; the
one in use is CHOSEN. LC viewer is white peaks on magenta and Spectra viewer a coloured band on near-black, so these
avoid both.

    python tools\\make_icon.py                 (writes structureeditor.ico next to the sources)
    python tools\\make_icon.py --sheet out.png (all designs side by side at 256 / 48 / 32 / 16 px)
"""
import math
import os
import sys

from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N = 1024                                                    # drawn large and scaled down, for smooth edges
CHOSEN = 'nitro'                                            # the user's pick (0.2.0); 'teal' was 0.1.0's
RED, BLUE, WHITE, INK = (240, 70, 60), (70, 130, 245), (255, 255, 255), (40, 30, 20)


def tile(top, bottom, radius=200):
    """A rounded square with a vertical gradient."""
    bg = Image.new('RGBA', (N, N))
    d = ImageDraw.Draw(bg)
    for y in range(N):
        t = y / N
        d.line((0, y, N, y), fill=tuple(int(top[k] + (bottom[k] - top[k]) * t) for k in range(3)))
    mask = Image.new('L', (N, N), 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, N - 1, N - 1), radius=radius, fill=255)
    img = Image.new('RGBA', (N, N), (0, 0, 0, 0))
    img.paste(bg, (0, 0), mask)
    return img


def hexagon(cx, cy, r, start=-90):
    return [(cx + r * math.cos(math.radians(start + 60 * k)), cy + r * math.sin(math.radians(start + 60 * k))) for k in range(6)]


def benzene(d, cx, cy, r, w, color=WHITE, inner=True):
    pts = hexagon(cx, cy, r)
    d.line(pts + [pts[0], pts[1]], fill=color, width=w, joint='curve')
    if inner:
        ins = hexagon(cx, cy, 0.72 * r)
        for k in (0, 2, 4):                                 # the inner lines of the double bonds
            (x1, y1), (x2, y2) = ins[k], ins[(k + 1) % 6]
            d.line((x1 + (x2 - x1) * 0.14, y1 + (y2 - y1) * 0.14, x2 - (x2 - x1) * 0.14, y2 - (y2 - y1) * 0.14),
                   fill=color, width=w)
    return pts


def dot(d, x, y, r, color):
    d.ellipse((x - r, y - r, x + r, y + r), fill=color)


def font(size, bold=True):
    for name in (('arialbd.ttf' if bold else 'arial.ttf'), 'DejaVuSans-Bold.ttf'):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def teal():
    """Benzene with a bond to a red O (phenol), on teal: the first icon."""
    img = tile((18, 110, 120), (8, 62, 74))
    d = ImageDraw.Draw(img)
    pts = benzene(d, 430, 590, 250, 44)
    x0, y0 = pts[1]
    x1, y1 = x0 + 190 * math.cos(math.radians(-30)), y0 + 190 * math.sin(math.radians(-30))
    d.line((x0, y0, x1, y1), fill=WHITE, width=44)
    dot(d, x1 + 40, y1 - 22, 95, RED)
    return img


def nitro():
    """Nitrobenzene (the product of the scheme in the user's ChemDraw file): blue N, two red O, on navy."""
    img = tile((32, 44, 92), (16, 22, 52))
    d = ImageDraw.Draw(img)
    pts = benzene(d, 512, 650, 220, 42)
    x0, y0 = pts[0]                                         # the top vertex, a bond straight up to N
    ny = y0 - 175
    d.line((x0, y0, x0, ny), fill=WHITE, width=42)
    for ang in (-150, -30):                                 # N=O and N-O
        ox, oy = x0 + 165 * math.cos(math.radians(ang)), ny + 165 * math.sin(math.radians(ang))
        d.line((x0, ny, ox, oy), fill=WHITE, width=42)
        dot(d, ox, oy, 66, RED)
    dot(d, x0, ny, 74, BLUE)
    return img


def scheme():
    """A ring, a reaction arrow and a second (gold) ring: a reaction scheme, on deep violet."""
    img = tile((84, 48, 140), (44, 22, 84))
    d = ImageDraw.Draw(img)
    benzene(d, 240, 512, 150, 34)
    d.line((430, 512, 640, 512), fill=WHITE, width=34)
    d.polygon([(710, 512), (610, 455), (635, 512), (610, 569)], fill=WHITE)
    p = hexagon(860, 512, 120)
    d.line(p + [p[0], p[1]], fill=(255, 205, 80), width=34, joint='curve')
    return img


def pen():
    """A hexagon being drawn: the last side a wedge with a pen nib at its open end, on amber."""
    img = tile((246, 170, 40), (214, 112, 18))
    d = ImageDraw.Draw(img)
    p = hexagon(450, 580, 260)
    d.line(p[:6], fill=INK, width=46, joint='curve')        # five sides
    a, b = p[5], p[0]                                       # the sixth as a wedge, wide end at the top vertex
    dx, dy = b[0] - a[0], b[1] - a[1]
    L = math.hypot(dx, dy)
    nx, ny = -dy / L * 40, dx / L * 40
    d.polygon([a, (b[0] + nx, b[1] + ny), (b[0] - nx, b[1] - ny)], fill=INK)
    tx, ty = p[0][0] + 70, p[0][1] - 40                     # a pen nib touching down near the top vertex
    d.polygon([(tx, ty), (tx + 150, ty - 250), (tx + 280, ty - 170)], fill=WHITE)
    d.line((tx + 30, ty - 30, tx + 120, ty - 150), fill=(214, 112, 18), width=16)
    return img


def mz():
    """A ring over 'm/z': drawing plus masses, white on green."""
    img = tile((30, 140, 90), (14, 88, 56))
    d = ImageDraw.Draw(img)
    benzene(d, 512, 390, 210, 42)
    d.text((512, 820), 'm/z', fill=WHITE, font=font(250), anchor='mm')
    return img


DESIGNS = {'teal': teal, 'nitro': nitro, 'scheme': scheme, 'pen': pen, 'mz': mz}


def draw(name=CHOSEN):
    return DESIGNS[name]()


def sheet(path):
    """Every design at 256 px (on light) and at 48, 32 and 16 px (on light and on a dark taskbar strip)."""
    col_w = 300
    out = Image.new('RGB', (col_w * len(DESIGNS), 470), (245, 246, 248))
    d = ImageDraw.Draw(out)
    d.rectangle((0, 400, out.width, 470), fill=(32, 34, 38))
    for k, (name, fn) in enumerate(DESIGNS.items()):
        img = fn()
        big = img.resize((256, 256), Image.LANCZOS)
        out.paste(big, (k * col_w + 22, 20), big)
        d.text((k * col_w + col_w // 2, 300), f'{k + 1}. {name}', fill=(30, 30, 30), font=font(26), anchor='mm')
        sx = k * col_w + 40
        for s in (48, 32, 16):
            small = img.resize((s, s), Image.LANCZOS)
            out.paste(small, (sx, 330 + (48 - s) // 2), small)
            out.paste(small, (sx, 411 + (48 - s) // 2), small)
            sx += s + 40
    out.save(path)


def main():
    if len(sys.argv) > 2 and sys.argv[1] == '--sheet':
        sheet(sys.argv[2])
        return
    img = draw()
    sizes = [(s, s) for s in (16, 20, 24, 32, 40, 48, 64, 128, 256)]
    img.resize((256, 256), Image.LANCZOS).save(os.path.join(HERE, 'structureeditor.ico'), sizes=sizes)


if __name__ == '__main__':
    main()
