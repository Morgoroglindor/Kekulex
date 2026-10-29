"""Draws structureeditor.ico (the program / installer icon): a white benzene ring with a bond to a red oxygen on a deep
teal rounded square (LC viewer is magenta, Spectra viewer near-black, so the three are told apart at a glance).

    python tools\\make_icon.py          (writes structureeditor.ico next to the sources)
"""
import math
import os

from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N = 1024                                                    # drawn large and scaled down, for smooth edges
TOP, BOTTOM = (18, 110, 120), (8, 62, 74)                   # background gradient


def draw():
    bg = Image.new('RGBA', (N, N))
    d = ImageDraw.Draw(bg)
    for y in range(N):
        t = y / N
        d.line((0, y, N, y), fill=tuple(int(TOP[k] + (BOTTOM[k] - TOP[k]) * t) for k in range(3)))
    mask = Image.new('L', (N, N), 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, N - 1, N - 1), radius=200, fill=255)
    img = Image.new('RGBA', (N, N), (0, 0, 0, 0))
    img.paste(bg, (0, 0), mask)
    d = ImageDraw.Draw(img)
    cx, cy, r, w = 430, 590, 250, 44
    pts = [(cx + r * math.cos(math.radians(-90 + 60 * k)), cy + r * math.sin(math.radians(-90 + 60 * k))) for k in range(6)]
    d.line(pts + [pts[0], pts[1]], fill='white', width=w, joint='curve')
    inner = [(cx + 0.72 * r * math.cos(math.radians(-90 + 60 * k)), cy + 0.72 * r * math.sin(math.radians(-90 + 60 * k)))
             for k in range(6)]
    for k in (0, 2, 4):                                     # the inner lines of the double bonds
        (x1, y1), (x2, y2) = inner[k], inner[(k + 1) % 6]
        d.line((x1 + (x2 - x1) * 0.14, y1 + (y2 - y1) * 0.14, x2 - (x2 - x1) * 0.14, y2 - (y2 - y1) * 0.14),
               fill='white', width=w)
    x0, y0 = pts[1]                                         # a bond up-right to an OH
    x1, y1 = x0 + 190 * math.cos(math.radians(-30)), y0 + 190 * math.sin(math.radians(-30))
    d.line((x0, y0, x1, y1), fill='white', width=w)
    d.ellipse((x1 - 95 + 40, y1 - 95 - 22, x1 + 95 + 40, y1 + 95 - 22), fill=(240, 70, 60))
    return img


def main():
    img = draw()
    sizes = [(s, s) for s in (16, 20, 24, 32, 40, 48, 64, 128, 256)]
    img.resize((256, 256), Image.LANCZOS).save(os.path.join(HERE, 'structureeditor.ico'), sizes=sizes)


if __name__ == '__main__':
    main()
