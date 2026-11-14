"""Pictures (PNG, JPEG, BMP, GIF ...) placed on the drawing: a Shape of kind 'image' whose params['png'] holds the
picture as base64 PNG, drawn into its box (x1, y1)-(x2, y2) - flipped when the box is. Pasted from the clipboard
(Ctrl+V: a copied picture, or a copied image file), dropped on the window or inserted from File > Insert picture."""
import base64
import hashlib

from PyQt6.QtCore import QBuffer, QByteArray, QIODevice
from PyQt6.QtGui import QImage

IMAGE_EXTS = ('.png', '.jpg', '.jpeg', '.bmp', '.gif', '.tif', '.tiff', '.webp')
LONGEST_PX = 4000                                          # bigger pictures are scaled down when they come in
LARGEST_PT = 320.0                                         # and placed at most this wide or high (pt)
_cache = {}


def encode(img):
    """base64 PNG of a QImage (scaled down past LONGEST_PX)."""
    if max(img.width(), img.height()) > LONGEST_PX:
        from PyQt6.QtCore import Qt
        img = img.scaled(LONGEST_PX, LONGEST_PX, Qt.AspectRatioMode.KeepAspectRatio,
                         Qt.TransformationMode.SmoothTransformation)
    data = QByteArray()
    buf = QBuffer(data)
    buf.open(QIODevice.OpenModeFlag.WriteOnly)
    img.save(buf, 'PNG')
    buf.close()
    return base64.b64encode(bytes(data)).decode('ascii')


def decode(b64):
    """The QImage of a stored picture (cached)."""
    key = hashlib.sha1(b64.encode('ascii')).hexdigest()
    if key not in _cache:
        img = QImage()
        img.loadFromData(base64.b64decode(b64))
        if len(_cache) > 64:
            _cache.clear()
        _cache[key] = img
    return _cache[key]


def natural_size(img, bond_length):
    """(width, height) in bond lengths at the picture's own resolution (96 dpi when it has none), at most LARGEST_PT."""
    dpm = img.dotsPerMeterX()
    dpi = dpm * 0.0254 if dpm and 30 < dpm * 0.0254 < 1200 else 96.0
    w, h = img.width() * 72.0 / dpi, img.height() * 72.0 / dpi
    k = min(1.0, LARGEST_PT / max(w, h, 1e-6))
    return w * k / bond_length, h * k / bond_length


def add(doc, img, x, y, bond_length):
    """Adds picture img with its top left at (x, y) (bond lengths); returns the shape id."""
    w, h = natural_size(img, bond_length)
    i = doc.add_shape(x, y, x + w, y + h, 'image')
    doc.shapes[i].params = {'png': encode(img), 'px': [img.width(), img.height()]}
    return i


def from_mime(md):
    """A QImage from clipboard / drop data: a picture itself, or the first image file among copied files; or None."""
    if md.hasImage():
        img = QImage(md.imageData())
        if not img.isNull():
            return img
    if md.hasUrls():
        for url in md.urls():
            if url.isLocalFile() and url.toLocalFile().lower().endswith(IMAGE_EXTS):
                img = QImage(url.toLocalFile())
                if not img.isNull():
                    return img
    return None


def paint(rend, p, i):
    """Draws picture shape i into its box (pt)."""
    from PyQt6.QtCore import QRectF
    sh = rend.doc.shapes[i]
    b64 = sh.params.get('png')
    if not b64:
        return
    img = decode(b64)
    if img.isNull():
        return
    if sh.x2 < sh.x1 or sh.y2 < sh.y1:                      # flipped with the drawing
        img = img.mirrored(sh.x2 < sh.x1, sh.y2 < sh.y1)
    BL = rend.BL
    x1, x2 = sorted((sh.x1 * BL, sh.x2 * BL))
    y1, y2 = sorted((sh.y1 * BL, sh.y2 * BL))
    from PyQt6.QtGui import QPainter
    p.save()
    p.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    p.drawImage(QRectF(x1, y1, x2 - x1, y2 - y1), img)
    p.restore()
