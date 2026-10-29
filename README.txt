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
  Hotkeys          ChemDraw's layout (taken from ChemDraw's hotkeys.xml; case matters, capital = Shift). Point at:
                   an atom - c n o s p f h i l(Cl) b(Br), Shift+C Cl, Shift+S Si, Shift+O OMe, Shift+N NO2,
                     Shift+P Ph, y Boc, e Et, m Me ... (groups are drawn out in full); 1 bond, 2 C=O, 4 / 5 wedge /
                     hash, 6 / 7 / v / u rings, a phenyl, z ethynyl, Shift+K tBu, 9 two methyls; + / - charge
                   a bond - 1 2 3 order, l / r / c double-bond side, w wedge, h hash, 4-8 / v fused ring,
                     a fused benzene, z fused cyclopentadiene
                   nothing - j benzene tool, x bond, Shift+X chain, e arrow, t text, Space select
                   Help > Keyboard shortcuts lists all of them.
  Select           click / Shift+click / drag a box; drag to move; drop one atom on another to merge; arrows nudge.
  SMILES box       paste a SMILES, Enter. Ctrl+V pastes SMILES or MOL text. File > Open reads .mol, .sdf, .smi.
  Clean up         Ctrl+Shift+K lays the selected structures out again (RDKit), keeping stereo.

ChemDraw-like features (0.3.0)
  Page view            the drawing sits on an A4 / Letter page (File > Document settings); View > Whole page.
  Nicknames            atom labels OMe, CO2H, Boc, Ph, NO2, OTs, TMS ... stay short on the drawing (turned round,
                       MeO2C, when bonded on the right) and count in full for formula / SMILES / MOL;
                       Structure > Expand label draws them out. Labels and captions are typed in place.
  Bond styles          wedge, hashed wedge, bold, hashed, dashed (dashed second line on a double bond), wavy.
  Arrows               reaction, equilibrium, retrosynthetic, resonance, dashed, curved, fish-hook (click a curved
                       arrow again with its tool to bend it the other way).
  Shapes               brackets, parentheses, braces, rectangle, rounded rectangle, oval, line.
  Captions             size, bold, italic, formula style (H2SO4 with subscripts), left / centre / right, colour.
  Object menu          align (edges, centres), distribute, Flip Horizontal / Vertical (the mirror image: stereo
                       inverts), Rotate 180 degrees Horizontal / Vertical (Ctrl+Shift+H / V: turned over, same
                       molecule), Rotate... (Ctrl+R), 15 degree steps (Ctrl+] / Ctrl+[), Scale... (Ctrl+K), colour.
  Lasso                the tool beside Select: drag a free-hand outline around what to select.
  Right-click menu     right-click an object (a right drag still pans): cut / copy / paste / clear, display
                       warnings, align, distribute, add frame (rectangle, brackets ...), flip, rotate, scale,
                       colour, the two colour palettes, then Atom (element, charge, label, highlight), Bond (order,
                       style, double-bond position, reverse, highlight), Text, Bracket / shape, Table (stoichiometry
                       grid), Molecule (select structure, clean up, expand labels, check, ring fill) and Analysis.
  Highlight / ring fill  View > Show Highlight Colors Toolbar and Show Ring Fill Colors Toolbar: a colour band behind
                       the selected atoms and bonds, a colour inside the selected structure's rings (Custom... any
                       colour, Remove takes it off). Kept in .sed files and in pictures; not written to .cdxml yet.
  Structure menu       clean up, expand label, analysis caption (formula, exact mass, m/z with isotopes, MW,
                       elemental analysis), analyze stoichiometry, check structure, name <-> structure through
                       PubChem (online, only when asked: sends the typed name or the structure's InChIKey),
                       R/S and E/Z labels, atom numbers, reaction interpretation.
  Stoichiometry grid   Structure > Analyze stoichiometry puts a Reactants / Products grid under the reaction
                       (structures left of the arrow and above / below it are reactants, right of it products).
                       Double-click it to type values, with or without units (250 mg, 1.5 mmol, 2 mL, 0.5 M,
                       1.2 g/mL, 95 %): limiting reactant, equivalents, sample mass, %weight, molarity, density,
                       volume, moles; for products equivalents, %completion, measured mass, purity. Moles, masses,
                       volumes, expected mass and %yield follow and are redrawn whenever the reaction changes.
                       Typed values are bold, calculated ones normal; rows can be hidden (right-click a row name).
  View menu            show reaction interpretation (reactants red, products magenta, above / below the arrow
                       cyan / blue), rulers in cm, whole page.
  Templates            rings and ring systems, heterocycles, functional groups, reagents, solvents.
  ChemDraw files       .cdxml opens and saves (structures, nicknames, bond styles, captions, arrows, shapes, colours,
                       document style); the stoichiometry grid and its typed values open too (not written back to
                       .cdxml yet: keep .sed for that). .cdx opens the structures (through RDKit). Not read: TLC
                       plates, orbitals, spectra.
  Analysis panel       Analysis tab: formula, exact mass, m/z with isotopes, MW, elemental analysis (each with
                       Paste, which puts it under the structure as a caption), adduct m/z, isotope pattern of the
                       chosen ion. Properties tab: logP, MR, tPSA, H-bond donors / acceptors, rotatable bonds,
                       rings, heavy atoms, fraction sp3, stereocentres, rule of five (RDKit).

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
  model.py         the drawing (atoms, bonds, arrows, text, shapes, nicknames) and its geometry; MOL; RDKit in
  cdxml.py         ChemDraw .cdxml reading and writing
  templates.py     the Templates menu (SMILES)
  stoich.py        stoichiometry grid: reaction columns, unit parsing, calculations; stoich_dialog.py its editor
  pubchem.py       name <-> structure through PubChem (online, on request only)
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
