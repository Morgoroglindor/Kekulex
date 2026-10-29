"""
theme.py - light / dark look of the Structure editor (copied from the Spectra viewer's theme.py).

    theme.T            current colours (a dict that is updated in place when the theme changes)
    theme.set('dark')  switch; theme.listeners are called by theme.apply() afterwards
    theme.make_app()   QApplication with the saved theme (Light, Dark or Follow system)

The look is Qt's Fusion style with a complete palette and one stylesheet (stylesheet()) built from T, set on each
top-level window (register_window), so every control - buttons, fields, drop-downs, tabs, scroll bars, menus, tables -
follows the theme. The small arrow and
check-mark images the stylesheet needs are drawn at start-up (no image files ship with the program).
The choice is kept per user (QSettings: the registry on Windows).
"""
import os
import sys
import tempfile

# Surfaces: window = the page behind everything, card = the panels on it (and the plot background, so plots sit on their
# card seamlessly), base = fields and tables, alt = secondary panels (control bars), button = push buttons.
# light / midlight / mid / dark / shadow: the bevel shades Fusion still uses for a few frames.
# accent_fg: the accent as text / marks on the page; accent2: the second colour (green) for structure - the card
# heading bars, the selected tab's underline, sliders; grad: a green curve; btn_off: a disabled coloured button;
# sel: selected rows / menu items; header: table headers; zebra: every other table row.
# The accent (blue, current_accent) marks what can be clicked: check boxes, focus, primary buttons, selections.
LIGHT = dict(name='light', window='#F1F4F6', card='#FFFFFF', base='#FFFFFF', alt='#F6F8F9', text='#1C2226',
             button='#FFFFFF', hover='#E9EEF2', pressed='#DCE4EA', border='#C7D0D8', muted='#66707A', rule='#DDE3E8',
             err='#B3261E', ok='#2E7D32', groove='#D5DDE3', tip='#FFFFFF', disabled='#A0A8AF',
             light='#FFFFFF', midlight='#EFF3F5', mid='#B2BCC4', dark='#97A3AD', shadow='#68737C',
             header='#F3F6F8', zebra='#F8FAFB', sel='#E2EEF8', scroll='#C3CDD5', scroll_hover='#9FADB8',
             fig_bg='white', fig_fg='black', fig_grid='0.88', ink='#222222', link='#1B6CA8',
             accent_fg='#1B6CA8', accent2='#2E8B57', grad='#2E7D32', btn_off='#D3DBE1', btn_off_text='#F6F8F9')
DARK = dict(name='dark', window='#121518', card='#1A1E22', base='#101316', alt='#1F2429', text='#E4E8EB',
            button='#262C32', hover='#2F363D', pressed='#39414A', border='#39414A', muted='#9DA7B0', rule='#2B3238',
            err='#FF8A80', ok='#81C784', groove='#38414A', tip='#262C32', disabled='#6B747C',
            light='#3B444D', midlight='#323A42', mid='#22282E', dark='#0A0C0E', shadow='#000000',
            header='#1F2429', zebra='#161A1E', sel='#1C3346', scroll='#39414A', scroll_hover='#4C5660',
            fig_bg='#1A1E22', fig_fg='#E4E8EB', fig_grid='0.25', ink='#F2F2F4', link='#7FB8E6',
            accent_fg='#7FB8E6', accent2='#5CC28A', grad='#66BB6A', btn_off='#2F363D', btn_off_text='#6B747C')

T = dict(LIGHT)
listeners = []                                             # callables run after a theme change
_windows = []                                              # top-level windows that carry the stylesheet (weak refs)
current_accent = '#1B6CA8'                                 # the accent colour (blue), used by the stylesheet


def set(name):
    T.clear()
    T.update(DARK if name == 'dark' else LIGHT)


def is_dark():
    return T['name'] == 'dark'


def mpl_rc():
    """Matplotlib rc settings that turn plots dark (an empty-effect set for the light theme)."""
    return {'figure.facecolor': T['fig_bg'], 'axes.facecolor': T['fig_bg'], 'savefig.facecolor': T['fig_bg'],
            'text.color': T['fig_fg'], 'axes.edgecolor': T['fig_fg'], 'axes.labelcolor': T['fig_fg'],
            'xtick.color': T['fig_fg'], 'ytick.color': T['fig_fg'], 'grid.color': T['fig_grid'],
            'legend.labelcolor': T['fig_fg']}


def lighten_ink(ink):
    """A very dark 'ink' (live trace colour) is replaced by the theme's light ink in the dark theme."""
    if not is_dark():
        return ink
    from matplotlib import colors
    r, g, b = colors.to_rgb(ink)
    return T['ink'] if 0.299 * r + 0.587 * g + 0.114 * b < 0.35 else ink


# ------------------------------------------------------------------------------------------------ Qt side
def settings():
    """The per-user settings of the program (theme, window positions, recent files): the registry on Windows.
    SED_SETTINGS_FILE=path.ini keeps them in that file instead (the tests use it to stay off the real settings)."""
    from PyQt6.QtCore import QSettings
    ini = os.environ.get('SED_SETTINGS_FILE')
    return QSettings(ini, QSettings.Format.IniFormat) if ini else QSettings('StructureEditor', 'StructureEditor')


def load_pref():
    v = settings().value('theme', 'system')
    return v if v in ('light', 'dark', 'system') else 'system'


def save_pref(v):
    settings().setValue('theme', v)


def restore_geometry(window, key):
    """Put a top-level window back where it was last closed (size, position, maximised); False the first time."""
    g = settings().value(f'geometry/{key}')
    return bool(g) and window.restoreGeometry(g)


def save_geometry(window, key):
    settings().setValue(f'geometry/{key}', window.saveGeometry())


def resolve(pref, app):
    if pref in ('light', 'dark'):
        return pref
    try:                                                   # Qt 6.5+: the operating system's setting
        from PyQt6.QtCore import Qt
        return 'dark' if app.styleHints().colorScheme() == Qt.ColorScheme.Dark else 'light'
    except Exception:
        return 'light'


def palette(accent=None):
    """A complete palette for the current theme: every role Fusion paints with, in every colour group."""
    from PyQt6.QtGui import QColor, QPalette
    accent = accent or current_accent
    p = QPalette()
    R, G = QPalette.ColorRole, QPalette.ColorGroup
    for role, key in ((R.Window, 'window'), (R.WindowText, 'text'), (R.Base, 'base'), (R.AlternateBase, 'zebra'),
                      (R.Text, 'text'), (R.Button, 'button'), (R.ButtonText, 'text'), (R.ToolTipBase, 'tip'),
                      (R.ToolTipText, 'text'), (R.PlaceholderText, 'disabled'), (R.Link, 'link'), (R.LinkVisited, 'link'),
                      (R.Light, 'light'), (R.Midlight, 'midlight'), (R.Mid, 'mid'), (R.Dark, 'dark'), (R.Shadow, 'shadow'),
                      (R.BrightText, 'err')):
        p.setColor(role, QColor(T[key]))
    for role in (R.Highlight, getattr(R, 'Accent', None)):     # Accent: Qt 6.6+
        if role is not None:
            p.setColor(role, QColor(accent))
    p.setColor(R.HighlightedText, QColor('#FFFFFF'))
    for role in (R.Text, R.ButtonText, R.WindowText):
        p.setColor(G.Disabled, role, QColor(T['disabled']))
    p.setColor(G.Disabled, R.Base, QColor(T['window']))
    p.setColor(G.Disabled, R.Button, QColor(T['window']))
    p.setColor(G.Disabled, R.Highlight, QColor(T['mid']))
    return p


# ------------------------------------------------------------------------------------------- stylesheet images
def _assets(accent):
    """Small images the stylesheet refers to (drop-down / spin arrows, check mark, radio dot), drawn for the current
    theme into a cache folder at 1x and 2x (Qt picks name@2x.png on high-DPI screens). Returns {name: path}."""
    from PyQt6.QtCore import QPointF, Qt
    from PyQt6.QtGui import QColor, QImage, QPainter, QPainterPath, QPen
    folder = os.path.join(tempfile.gettempdir(), 'structureeditor_style', f"{T['name']}_{accent.lstrip('#')}")
    os.makedirs(folder, exist_ok=True)

    def chevron(p, s, up):                                 # a 10 x 10 chevron, pointing down (or up)
        y0, y1 = (6.5, 3.5) if up else (3.5, 6.5)
        path = QPainterPath(QPointF(2.5 * s, y0 * s)); path.lineTo(5 * s, y1 * s); path.lineTo(7.5 * s, y0 * s)
        p.drawPath(path)

    def check(p, s):
        path = QPainterPath(QPointF(2.3 * s, 5.2 * s)); path.lineTo(4.2 * s, 7.1 * s); path.lineTo(7.8 * s, 3.0 * s)
        p.drawPath(path)

    def dot(p, s):
        p.setBrush(p.pen().color()); p.setPen(Qt.PenStyle.NoPen); p.drawEllipse(QPointF(5 * s, 5 * s), 2.2 * s, 2.2 * s)

    jobs = {'down': (chevron, False, T['muted']), 'up': (chevron, True, T['muted']),
            'down_off': (chevron, False, T['disabled']), 'up_off': (chevron, True, T['disabled']),
            'check': (check, None, '#FFFFFF'), 'check_off': (check, None, T['btn_off_text']), 'dot': (dot, None, '#FFFFFF')}
    out = {}
    for name, (fn, arg, col) in jobs.items():
        for scale, suffix in ((1, ''), (2, '@2x')):
            path = os.path.join(folder, f'{name}{suffix}.png')
            if not os.path.exists(path):
                img = QImage(10 * scale, 10 * scale, QImage.Format.Format_ARGB32_Premultiplied)
                img.fill(Qt.GlobalColor.transparent)
                p = QPainter(img); p.setRenderHint(QPainter.RenderHint.Antialiasing)
                pen = QPen(QColor(col), 1.5 * scale); pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin); p.setPen(pen)
                fn(p, scale) if arg is None else fn(p, scale, arg)
                p.end(); img.save(path)
        out[name] = os.path.join(folder, f'{name}.png').replace('\\', '/')
    return out


def stylesheet(accent=None):
    """The application stylesheet for the current theme."""
    from PyQt6.QtGui import QColor
    a = accent or current_accent
    A = _assets(a)
    a_hover, a_dark = QColor(a).lighter(118).name(), QColor(a).darker(115).name()
    t = dict(T, a=a, a_hover=a_hover, a_dark=a_dark, **{f'img_{k}': v for k, v in A.items()})
    return """
QMainWindow, QDialog {{ background: {window}; }}
QToolTip {{ background: {tip}; color: {text}; border: 1px solid {border}; padding: 5px 7px; border-radius: 4px; }}

/* buttons */
QPushButton {{ background: {button}; color: {text}; border: 1px solid {border}; border-radius: 5px; padding: 3px 10px;
              min-height: 18px; }}
QPushButton:hover {{ background: {hover}; }}
QPushButton:pressed {{ background: {pressed}; }}
QPushButton:checked {{ background: {sel}; border-color: {a}; color: {accent_fg}; }}
QPushButton:disabled {{ background: {window}; color: {disabled}; border-color: {rule}; }}
QPushButton:default {{ border-color: {a}; }}

/* text fields, spin boxes, drop-downs */
QLineEdit, QAbstractSpinBox, QComboBox {{ background: {base}; color: {text}; border: 1px solid {border}; border-radius: 4px;
              padding: 2px 6px; min-height: 18px; selection-background-color: {a}; selection-color: white; }}
QLineEdit:hover, QAbstractSpinBox:hover, QComboBox:hover {{ border-color: {scroll_hover}; }}
QLineEdit:focus, QAbstractSpinBox:focus, QComboBox:focus, QComboBox:on {{ border-color: {a}; }}
QLineEdit:disabled, QAbstractSpinBox:disabled, QComboBox:disabled {{ background: {window}; color: {disabled}; border-color: {rule}; }}
QAbstractSpinBox {{ padding-right: 16px; }}
QAbstractSpinBox::up-button, QAbstractSpinBox::down-button {{ subcontrol-origin: border; width: 16px; border: none;
              background: transparent; }}
QAbstractSpinBox::up-button {{ subcontrol-position: top right; }}
QAbstractSpinBox::down-button {{ subcontrol-position: bottom right; }}
QAbstractSpinBox::up-button:hover, QAbstractSpinBox::down-button:hover {{ background: {hover}; }}
QAbstractSpinBox::up-arrow {{ image: url({img_up}); width: 8px; height: 8px; }}
QAbstractSpinBox::down-arrow {{ image: url({img_down}); width: 8px; height: 8px; }}
QAbstractSpinBox::up-arrow:disabled, QAbstractSpinBox::up-arrow:off {{ image: url({img_up_off}); }}
QAbstractSpinBox::down-arrow:disabled, QAbstractSpinBox::down-arrow:off {{ image: url({img_down_off}); }}
QComboBox {{ padding-right: 20px; }}
QComboBox::drop-down {{ subcontrol-origin: padding; subcontrol-position: center right; width: 20px; border: none; }}
QComboBox::down-arrow {{ image: url({img_down}); width: 10px; height: 10px; }}
QComboBox::down-arrow:disabled {{ image: url({img_down_off}); }}
QComboBox QAbstractItemView {{ background: {card}; border: 1px solid {border}; padding: 3px; outline: none;
              selection-background-color: {sel}; selection-color: {text}; }}

/* check boxes, radio buttons */
QCheckBox, QRadioButton {{ spacing: 7px; }}
QCheckBox::indicator, QRadioButton::indicator {{ width: 14px; height: 14px; border: 1px solid {border}; background: {base}; }}
QCheckBox::indicator {{ border-radius: 4px; }}
QRadioButton::indicator {{ border-radius: 8px; }}
QCheckBox::indicator:hover, QRadioButton::indicator:hover {{ border-color: {a}; }}
QCheckBox::indicator:checked {{ background: {a}; border-color: {a}; image: url({img_check}); }}
QRadioButton::indicator:checked {{ background: {a}; border-color: {a}; image: url({img_dot}); }}
QCheckBox::indicator:disabled, QRadioButton::indicator:disabled {{ background: {window}; border-color: {rule}; }}
QCheckBox::indicator:checked:disabled {{ background: {btn_off}; border-color: {btn_off}; image: url({img_check_off}); }}
QCheckBox:disabled, QRadioButton:disabled {{ color: {disabled}; }}

/* tabs */
QTabWidget::pane {{ background: {card}; border: 1px solid {rule}; border-radius: 6px; top: -1px; }}
QTabBar {{ qproperty-drawBase: 0; }}
QTabBar::tab {{ background: transparent; color: {muted}; border: none; border-bottom: 2px solid transparent;
              padding: 6px 11px 5px 11px; margin-right: 1px; }}
QTabBar::tab:hover {{ color: {text}; background: {hover}; border-top-left-radius: 4px; border-top-right-radius: 4px; }}
QTabBar::tab:selected {{ color: {text}; border-bottom: 2px solid {accent2}; background: transparent; }}
QTabBar::tab:disabled {{ color: {disabled}; }}
QTabBar QToolButton {{ background: {card}; border: 1px solid {rule}; border-radius: 3px; }}

/* scroll bars */
QScrollBar:vertical {{ background: transparent; width: 12px; margin: 0; }}
QScrollBar:horizontal {{ background: transparent; height: 12px; margin: 0; }}
QScrollBar::handle:vertical {{ background: {scroll}; border-radius: 4px; min-height: 32px; margin: 2px 3px; }}
QScrollBar::handle:horizontal {{ background: {scroll}; border-radius: 4px; min-width: 32px; margin: 3px 2px; }}
QScrollBar::handle:hover {{ background: {scroll_hover}; }}
QScrollBar::add-line, QScrollBar::sub-line {{ width: 0; height: 0; border: none; background: none; }}
QScrollBar::add-page, QScrollBar::sub-page {{ background: none; }}
QAbstractScrollArea::corner {{ background: transparent; }}

/* menus, tool bar, status bar */
QMenuBar {{ background: {window}; color: {text}; padding: 2px 4px; }}
QMenuBar::item {{ background: transparent; padding: 4px 9px; border-radius: 4px; }}
QMenuBar::item:selected, QMenuBar::item:pressed {{ background: {hover}; }}
QMenu {{ background: {card}; color: {text}; border: 1px solid {border}; padding: 4px; }}
QMenu::item {{ padding: 5px 26px 5px 24px; border-radius: 4px; }}
QMenu::item:selected {{ background: {sel}; }}
QMenu::item:disabled {{ color: {disabled}; }}
QMenu::separator {{ height: 1px; background: {rule}; margin: 4px 8px; }}
QMenu::indicator {{ width: 13px; height: 13px; left: 6px; border-radius: 3px; }}
QMenu::indicator:checked {{ background: {a}; image: url({img_check}); }}
QToolBar {{ background: {window}; border: none; border-bottom: 1px solid {rule}; padding: 3px 8px; spacing: 2px; }}
QToolBar::separator {{ background: {rule}; width: 1px; margin: 5px 6px; }}
QToolButton {{ background: transparent; color: {text}; border: 1px solid transparent; border-radius: 5px; padding: 4px 7px; }}
QToolButton:hover {{ background: {hover}; border-color: {rule}; }}
QToolButton:pressed, QToolButton:checked {{ background: {pressed}; }}
QToolButton:disabled {{ color: {disabled}; }}
QStatusBar {{ background: {window}; color: {muted}; border-top: 1px solid {rule}; }}
QStatusBar::item {{ border: none; }}
QStatusBar QLabel {{ color: {muted}; padding: 0 6px; }}

/* lists, tables, text views */
QAbstractItemView {{ background: {base}; color: {text}; border: 1px solid {rule}; border-radius: 6px;
              alternate-background-color: {zebra}; selection-background-color: {sel}; selection-color: {text};
              gridline-color: {rule}; outline: none; }}
QHeaderView::section {{ background: {header}; color: {muted}; border: none; border-bottom: 1px solid {rule};
              border-right: 1px solid {rule}; padding: 4px 6px; font-weight: 600; }}
QTableCornerButton::section {{ background: {header}; border: none; border-bottom: 1px solid {rule}; }}
QPlainTextEdit, QTextEdit {{ background: {base}; color: {text}; border: 1px solid {rule}; border-radius: 6px;
              selection-background-color: {a}; }}
QGroupBox {{ border: 1px solid {rule}; border-radius: 6px; margin-top: 16px; padding: 8px 6px 6px 6px; }}
QGroupBox::title {{ subcontrol-origin: margin; subcontrol-position: top left; left: 10px; padding: 0 4px; color: {muted};
              font-weight: 600; }}
/* the viewer's own pieces (qtwidgets gives them these object names / roles) */
QWidget#pdaCard {{ background: {card}; border: 1px solid {rule}; border-radius: 8px; }}
QWidget#pdaNav {{ background: {alt}; border: 1px solid {rule}; border-radius: 6px; }}
QScrollArea#pdaScroll {{ background: {base}; border: 1px solid {rule}; border-radius: 6px; }}
QWidget#pdaScrollInner {{ background: {base}; }}
QPlainTextEdit#pdaOutput {{ color: {err}; }}
QLabel#pdaH {{ font-size: 10.5pt; font-weight: 600; color: {text}; border: 0px solid transparent;
              border-left: 3px solid {accent2}; padding: 1px 0 1px 8px; }}
QPushButton[pdaPrimary="true"] {{ background: {a}; color: white; border: 1px solid {a}; }}
QPushButton[pdaPrimary="true"]:hover {{ background: {a_hover}; border-color: {a_hover}; }}
QPushButton[pdaPrimary="true"]:pressed {{ background: {a_dark}; border-color: {a_dark}; }}
QPushButton[pdaPrimary="true"]:disabled {{ background: {btn_off}; border-color: {btn_off}; color: {btn_off_text}; }}
QPushButton[pdaStep="true"] {{ padding: 0px; font-size: 14pt; }}
QSlider {{ min-height: 20px; }}
QSlider::groove:horizontal {{ height: 4px; background: {groove}; border-radius: 2px; }}
QSlider::sub-page:horizontal {{ background: {accent2}; border-radius: 2px; }}
QSlider::handle:horizontal {{ background: {a}; width: 12px; height: 12px; margin: -6px 0; border-radius: 8px;
              border: 2px solid {card}; }}
QSlider::handle:horizontal:hover {{ background: {a_hover}; }}
QSlider::sub-page:horizontal:disabled, QSlider::handle:horizontal:disabled {{ background: {disabled}; }}
QSlider:vertical {{ min-width: 20px; }}
QSlider::groove:vertical {{ width: 4px; background: {groove}; border-radius: 2px; }}
QSlider::add-page:vertical {{ background: {accent2}; border-radius: 2px; }}
QSlider::handle:vertical {{ background: {a}; width: 12px; height: 12px; margin: 0 -6px; border-radius: 8px;
              border: 2px solid {card}; }}
QSlider::handle:vertical:hover {{ background: {a_hover}; }}
QSlider::add-page:vertical:disabled, QSlider::handle:vertical:disabled {{ background: {disabled}; }}
QSplitter::handle {{ background: {window}; }}
QSplitter::handle:hover {{ background: {a}; }}
QSplitter::handle:horizontal {{ width: 5px; }}
QSplitter::handle:vertical {{ height: 5px; }}
QProgressBar {{ background: {base}; border: 1px solid {rule}; border-radius: 4px; text-align: center; color: {text}; }}
QProgressBar::chunk {{ background: {accent2}; border-radius: 3px; }}
QScrollArea {{ background: transparent; border: none; }}
""".format(**t)


_applying = False
current_pref = 'system'                                    # what apply() last applied (not necessarily what is saved)


def apply(app, pref=None, accent=None, notify=True):
    """Apply 'light' / 'dark' / 'system' to the application palette and stylesheet (and tell the open windows)."""
    global _applying, current_pref, current_accent
    if _applying:                                          # setColorScheme below re-enters through colorSchemeChanged
        return
    _applying = True
    try:
        pref = current_pref = pref or load_pref()
        current_accent = accent or current_accent
        _set_scheme(app, pref)
        set(resolve(pref, app))
        app.setPalette(palette())
        _restyle_windows()
        if notify:
            for fn in list(listeners):
                try:
                    fn()
                except RuntimeError:                       # a window that no longer exists
                    pass
    finally:
        _applying = False


def register_window(window):
    """Style a top-level window (and everything in it: its menus, drop-down lists and dialogs) with the stylesheet.
    The stylesheet sits on each window, not on the application: replacing an application-wide stylesheet makes Qt
    re-style every widget several times over (~1.5 s for the viewer), a window's own takes a fraction of that."""
    import weakref
    _windows.append(weakref.ref(window))
    window.setStyleSheet(stylesheet())


def _restyle_windows():
    css = stylesheet()
    for ref in list(_windows):
        win = ref()
        if win is None:
            _windows.remove(ref); continue
        try:
            win.setStyleSheet(css)
        except RuntimeError:                               # deleted on the Qt side
            _windows.remove(ref)


def set_accent(color):
    """The viewer's accent colour changed (Display tab): restyle the controls that use it (focus rings, check boxes,
    selected tabs). Before the QApplication exists it is only remembered."""
    global current_accent
    if color == current_accent:
        return
    current_accent = color
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance()
    if app is not None and not _applying:
        app.setPalette(palette())
        _restyle_windows()


def _set_scheme(app, pref):
    """Tell Qt which scheme the program uses (Qt 6.8+): on Windows this makes the title bar and the native menus
    and dialogs dark or light along with the window; 'system' hands the choice back to the operating system."""
    try:
        from PyQt6.QtCore import Qt
        hints = app.styleHints()
        want = {'light': Qt.ColorScheme.Light, 'dark': Qt.ColorScheme.Dark}.get(pref, Qt.ColorScheme.Unknown)
        if hasattr(hints, 'setColorScheme'):
            hints.setColorScheme(want)
    except (ImportError, AttributeError):
        pass


def resource(name):
    """Path of a data file shipped with the program (next to the sources, or inside the frozen bundle)."""
    for base in (getattr(sys, '_MEIPASS', None), os.path.dirname(os.path.abspath(sys.executable if getattr(sys, 'frozen', False) else __file__))):
        if base and os.path.isfile(os.path.join(base, name)):
            return os.path.join(base, name)
    return None


def ui_font():
    """The interface font: Segoe UI on Windows (the system's own), else the platform default."""
    from PyQt6.QtGui import QFont, QFontDatabase
    f = QFontDatabase.systemFont(QFontDatabase.SystemFont.GeneralFont)
    if sys.platform == 'win32' and 'Segoe UI' in QFontDatabase.families():
        f = QFont('Segoe UI')
    f.setPointSizeF(9.0)
    f.setHintingPreference(QFont.HintingPreference.PreferNoHinting)
    return f


def make_app(argv, name, accent=None):
    from PyQt6 import QtGui, QtWidgets
    app = QtWidgets.QApplication(argv)
    app.setApplicationName(name)
    app.setStyle('Fusion')                                 # same look everywhere
    app.setFont(ui_font())
    apply(app, accent=accent, notify=False)
    try:                                                   # follow the operating system when 'system' is chosen
        from PyQt6.QtCore import QTimer
        # Only when following the system: the signal also fires when apply() itself picks light or dark (the saved
        # choice may lag behind). Deferred: Qt resets the palette for the new system scheme after this signal.
        app.styleHints().colorSchemeChanged.connect(
            lambda _s: current_pref == 'system' and QTimer.singleShot(
                0, lambda: current_pref == 'system' and apply(app, 'system')))
    except AttributeError:                                 # Qt before 6.5: no system scheme to follow
        pass
    ico = resource('structureeditor.ico')
    if ico:
        app.setWindowIcon(QtGui.QIcon(ico))
    return app


def add_menu(window, menubar_menu, on_change=None):
    """Adds Theme > Light / Dark / Follow system to a menu; returns the actions."""
    from PyQt6.QtGui import QAction, QActionGroup
    sub = menubar_menu.addMenu('&Theme')
    grp = QActionGroup(window); grp.setExclusive(True)
    cur = load_pref()
    for label, key in (('&Light', 'light'), ('&Dark', 'dark'), ('&Follow system', 'system')):
        a = QAction(label, window, checkable=True); a.setChecked(key == cur); grp.addAction(a); sub.addAction(a)
        a.triggered.connect(lambda _c=False, k=key: (save_pref(k), apply(window_app(window), k), on_change and on_change()))
    sub.theme_actions = {k: a for a, (_l, k) in zip(grp.actions(), (('', 'light'), ('', 'dark'), ('', 'system')))}
    return sub


def window_app(window):
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance()
