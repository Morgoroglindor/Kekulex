Structure editor - draw chemical structures for figures and slides, and read off their formula, exact mass and the
m/z of the usual LC-MS adducts. A sibling of LC viewer and Spectra viewer.

Run from source:   pip install -r requirements.txt  ;  python structure_editor.py [drawing.sed | file.mol | file.sdf | file.smi]
Windows installer: double-click BUILD_INSTALLER.bat; it creates installer_output\StructureEditor_Setup.exe
Tests:             pip install -r requirements-dev.txt  ;  python -m pytest

Drawing (Help > How to draw in the program says the same)
  Bond tools       click empty space, click an atom to grow a bond (zigzag), drag to point it (15 degree steps); drop the
                   end on an atom to close a ring. Click a bond: single -> double -> triple. Double-bond tool on a double
                   bond: moves the second line (automatic -> other side -> centred).
  Wedge / hash     start at the stereocentre; click again to reverse. The SMILES and InChIKey carry the stereo.
  Chain            drag; the number of bonds follows the mouse.
  Rings            empty space: new ring; a bond: fused ring (on the outside); an atom: attached ring.
  Keys over atoms  C N O S P F H I, L = Cl, B = Br, + / - charge, Delete. Over bonds: 1 2 3.
  Select           click / Shift+click / drag a box; drag to move; drop one atom on another to merge; arrows nudge.
  SMILES box       paste a SMILES, Enter. Ctrl+V pastes SMILES or MOL text. File > Open reads .mol, .sdf, .smi.
  Clean up         Ctrl+Shift+K lays the selected structures out again (RDKit), keeping stereo.

Figures
  Ctrl+Shift+C copies the selection (or everything) as a picture (PNG + SVG on the clipboard) for PowerPoint / Word;
  Format > Picture resolution and Transparent background apply. File > Export picture writes PNG, SVG or PDF.
  The default drawing style is ACS 1996 (14.4 pt bonds, 0.6 pt lines, Arial 10 pt); Format > Drawing style changes it
  per drawing. Text is drawn as outlines, so pictures look the same everywhere.

Masses (panel on the right)
  Formula (Hill order), monoisotopic exact mass, average MW, charge, and m/z of [M+H]+, [M+NH4]+, [M+Na]+, [M+K]+,
  [M+H-H2O]+, [M+2H]2+, [2M+H]+, [2M+Na]+ and [M-H]-, [M+Cl]-, [M+HCOO]-, [M+CH3COO]-, [M-H-H2O]-, [M-2H]2-, [2M-H]-.
  m/z = (neutral atoms + adduct - z electrons) / |z|. A structure drawn with a net charge is treated as the ion [M]z.
  With several structures, pick one in the list, or select atoms. Ctrl+C in the table copies rows.
  Atoms over their usual valence get a red box.

Files
  structure_editor.py   window + entry point: tools, menus, files, clipboard, export, style dialog
  canvas.py        the drawing area: tools, mouse, keys, selection, undo / redo, clean up
  model.py         the drawing (atoms, bonds, arrows, text) and its geometry; MOL writing; RDKit molecules in
  render.py        drawing with QPainter in points (screen and every export the same); PNG / SVG / PDF
  chem.py          formula, masses, adducts, SMILES, InChIKey (RDKit); reading SMILES / MOL / SDF
  info_panel.py    the formula and masses panel
  theme.py, applog.py   look (light / dark) and log file + crash dialog - copied from Spectra viewer
  icons.py         tool icons drawn at start-up
  version.py       the version number (About box and installer)
  tools\make_icon.py   draws structureeditor.ico

The .sed file is JSON: atoms (x, y in bond lengths, y down; element; charge), bonds (atoms, order 1-3, stereo
wedge / hash with the narrow end at atom a, double-bond side), arrows, texts and the drawing style.
Log: %APPDATA%\StructureEditor\structureeditor.log
