Kekulex for Android phones and tablets

What it is
  A touch version of the drawing side of Kekulex, written in JavaScript around RDKit's WebAssembly build (RDKit.js), and
  packaged as an Android app by Capacitor. It works offline, except the name lookup. Python and Pyodide are not used:
  Pyodide has no RDKit. Drawings are the desktop's .sed files, both ways.

  www/core/doc.js      the drawing (as model.py)
  www/core/chem.js     formula, masses, ions, SMILES, InChIKey (as chem.py)
  www/core/layout.js   usual drawings: skeleton templates, even bonds, ring-fusion H, checked wedges (as layout2d.py)
  www/core/draw.js     painting on a canvas (as render.py)
  www/core/condensed.js typed labels read as condensed formulas: SO3K, N(CH3)2, COCH3 ... (as condensed.py)
  www/core/data.js     made from the desktop's Python by tools/export_data.py: do not edit
  www/app.js           tools, gestures, sheets, files
  tests/core.test.mjs  against the desktop's results (tests/reference.json, also from tools/export_data.py)

On the phone
  Tools (bottom): select, bonds, rings, atoms, chain, charge, eraser, arrow, text. The row above holds the tool's choices.
  Bonds: tap an atom for a new bond, drag for its direction, tap a bond to change it (single: cycles 1, 2, 3).
  Rings: tap empty space, an atom or a bond (fused).   Atoms: tap an atom. Long-press any atom: type a label or pick a group.
  Select: tap a molecule, or drag a box; drag it, also from inside the dashed box, to move it. The bar above has Delete, Copy,
  Clean up, SMILES, Flip and rotate.   Two fingers: move and zoom.
  Top: the menu; the formula (tap it: masses, adducts, SMILES, InChIKey; it covers the selection, if any); undo, redo.
  Menu: structure from a name (PubChem, only the name is sent), paste SMILES / MOL, templates, clean up, open and send files.
  The drawing is kept between sessions.
  Labels: elements, nicknames and condensed formulas (SO3K, CH2CH3, N(CH3)2, PO3H2 ...), also in lower case (so3k);
  CH2 keeps two H, CH2• / C•• / O•+ make radicals, carbenes and radical cations, as on the desktop. Groups made on the
  desktop select and move as one here, and captions pasted there (formula, masses, MS fragments) follow edits made here
  (names, properties and isotope patterns keep their text until the desktop updates them).

Desktop only
  NMR prediction and the 3D tab need 3D conformers (ETKDG + MMFF), which RDKit.js does not have; so do MD exports. Also .cdxml /
  .cdx, stoichiometry, TLC plates, bio objects and sequences: these are kept when a .sed is opened and saved here, not shown.

Build (Windows; needs Node.js, JDK 21, Android SDK with platform 36 + build-tools 36)
  cd mobile
  npm install
  npm test                 chemistry against the desktop
  npm run apk              -> dist\Kekulex-<version>-debug.apk
  npm run serve            tries it in a PC browser on http://localhost:8766 (run npm run prepare-web first)
After changing the desktop's model.py, chem.py, templates.py or layout2d.py:
  ..\.venv\Scripts\python tools\export_data.py
Icons: ..\.venv\Scripts\python tools\make_android_icons.py (run after a fresh "npx cap add android").
The debug APK is for your own devices. Google Play needs a signed release build (assembleRelease and your own key).
