"""Android launcher icons from the desktop program's kekulex.ico.   ..\\.venv\\Scripts\\python tools\\make_android_icons.py"""
import os

from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
MOBILE = os.path.dirname(HERE)
RES = os.path.join(MOBILE, 'android', 'app', 'src', 'main', 'res')
ico = Image.open(os.path.join(os.path.dirname(MOBILE), 'kekulex.ico'))
ico.size = max(ico.ico.sizes())
ico.load()
img = ico.convert('RGBA')
for name, size in (('mdpi', 48), ('hdpi', 72), ('xhdpi', 96), ('xxhdpi', 144), ('xxxhdpi', 192)):
    folder = os.path.join(RES, f'mipmap-{name}')
    square = img.resize((size, size), Image.LANCZOS)
    square.save(os.path.join(folder, 'ic_launcher.png'))
    square.save(os.path.join(folder, 'ic_launcher_round.png'))
    # the adaptive icon: 108 dp with the icon in the safe middle 72 dp, on the icon's own background colour
    full = int(size * 108 / 48 + 0.5)
    fg = Image.new('RGBA', (full, full), img.getpixel((img.width // 2, 2)))
    inner = int(size * 72 / 48 + 0.5)
    icon = img.resize((inner, inner), Image.LANCZOS)
    fg.alpha_composite(icon, ((full - inner) // 2, (full - inner) // 2))
    fg.save(os.path.join(folder, 'ic_launcher_foreground.png'))
print('icons written')
