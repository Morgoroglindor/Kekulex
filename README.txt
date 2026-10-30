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
                       elemental analysis), analyze stoichiometry, check structure, Convert Structure to Name
                       (Ctrl+Alt+N: systematic IUPAC name made offline by Openclatura, MIT licence, rule-based on
                       the IUPAC 2013 recommendations, stereo included; also shown as Name in the Analysis tab -
                       check it before publishing), name -> structure and common names through PubChem (online,
                       only when asked: sends the typed name or the structure's InChIKey), R/S and E/Z labels,
                       atom numbers, reaction interpretation.
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
  NMR prediction       Structure > Predict 1H / 13C NMR (also right-click > Analysis): CSP5 (MIT, machine learning
                       on measured spectra, solvent aware: CDCl3, DMSO-d6, CD3OD, D2O ...). A window with the
                       simulated spectrum (mouse as in Spectra viewer: drag a box to zoom, wheel zooms ppm, Ctrl+wheel
                       the height, right drag pans, double-click / F whole spectrum, click picks a signal; the splitter
                       gives it more room; saved as PNG / SVG / PDF), the signals (shift +/- uncertainty, first-order
                       multiplicity with expected J: vicinal 7 Hz, ring CH-CH by Karplus on an RDKit 3D conformer,
                       aromatic ortho 8 / meta 2, alkene cis 10-11 / trans 16 ...), spectrometer 300-800 MHz, the
                       experimental-section line with J to copy, Label atoms, Paste text. OH / NH shifts are flagged as
                       unreliable, and left out in D2O / CD3OD. Optional: install as in requirements-nmr.txt (CSP5 is patched to install on
                       Windows by tools\patch_csp5.py). The first prediction loads the model (a few seconds).
  Charges              as in ChemDraw: + / - change the hydrogens by the octet (NH2 -> NH3+ -> NH2 2+, C+ = CH+ with
                       its C shown); the sign goes on the atom's symbol, at the first spot clear of bonds and H
                       (top right, bottom right, top left, bottom left, sides, else the widest gap). Drag the sign
                       with Select / Lasso to put it elsewhere (it follows its atom); right-click > Atom > Reset
                       charge position.
  Bond-count warnings  a red box on atoms with more bonds than they take (C 4, neutral N 3, N+ 4, O 2, O+ 3 ...;
                       a nickname such as OMe only one); pointing at one says what is wrong, Check structure lists
                       them all. View > Display warnings turns the boxes off.
  Deleting             Delete (or the eraser) over N, Cl, OMe, a charged atom ... turns it back into a carbon and
                       keeps its bonds; over a plain carbon it removes the atom.
  Style sheets         File > New from Style Sheet and Format > Apply Style Sheet: ACS 1996, New Document, New Slide,
                       RSC, Wiley, Nature Research, Synthesis / Synlett, Science of Synthesis, Adv. Synth. Catal.,
                       Helv. Chim. Acta, J. Chin. Chem. Soc., J. Het. Chem., J. Mol. Mod., Phytomedicine (bond length,
                       line widths, margins, label and caption fonts as ChemDraw's style sheets set them).
  Page and printing    File > Document settings: page A4, A3, A5, B5, Letter, Legal, Tabloid, portrait or landscape;
                       caption font and size, bold labels, terminal CH3 labels. File > Print (Ctrl+P), Print Preview.

Biology (like ChemDraw's HELM / biopolymer tools and BioDraw; this program's own artwork and colours)
  Sequences            the sequence tool (three beads, under Text) or Biopolymer > Insert sequence (Ctrl+Shift+Q):
                       type or paste a peptide, DNA or RNA as FASTA one-letter codes, three-letter text
                       (H-Ala-Gly-Pro-NH2, Ac-...-OH) or HELM (PEPTIDE1{A.C.D}$$$$, RNA1{R(A)P.R(C)}$$$$). Drawn as
                       monomers coloured by residue class (peptides circles, DNA / RNA rounded squares), wrapped every
                       N residues, with H- / -OH (Ac, NH2) or 5' / 3' (5'-p, p-3') ends. Peptides: cyclic (head to
                       tail), disulfide (Cys-Cys, also Sec) and lactam (Lys-Asp / Lys-Glu) bridges drawn as arcs; the
                       22 amino acids incl. Sec (U) and Pyl (O); B Z J X allowed but have no structure.
                       Double-click to edit; point at a residue for its name and number.
  Biopolymer menu      Paste Special (FASTA peptide / DNA / RNA, biopolymer text, HELM - a pasted HELM string becomes
                       a sequence anyway; U pasted as DNA becomes T and back, as in ChemDraw), Copy As HELM / HELM
                       natural analog (FASTA) / three-letter text, toggle A / Ala labels, residue numbers, residues
                       per line, Expand to structure (the full molecule, one undo step), Show analysis (formula,
                       exact mass, MW, [M+zH]z+ of peptides / [M-zH]z- of oligonucleotides as a caption).
  Masses of sequences  pick "Sequence n" in the masses panel: formula, exact mass, adducts and the multiply charged
                       ions (z up to the basic sites + 1 for peptides), isotope pattern. Peptides are built residue
                       by residue (checked against RDKit's MolFromSequence); DNA / RNA through RDKit.
  BioDraw objects      the bio drop-down under the sequence tool (also Biopolymer > BioDraw objects): membrane
                       bilayer, membrane arc, micelle, DNA double helix, helix protein, receptor, ion channel,
                       enzyme, G protein (α β γ), antibody, ribosome, tRNA, mRNA, mitochondrion, nucleus, cell,
                       bacterium, virus, cloud, plasmid map. Click for the default size or drag; membranes, DNA,
                       helices and mRNA run along the drag and add subunits as they get longer. Right-click > Bio
                       object: fill colour, thicker / thinner, plasmid map (size in bp, name, regions as arrows,
                       markers such as restriction sites; also by double-click).
  Templates            amino acids (22), sugars (D-glucose, galactose, mannose α / β pyranoses, fructose, ribose,
                       deoxyribose, open-chain glucose, sucrose), nucleobases, nucleosides, deoxynucleosides, NMPs,
                       ATP, NAD+. (Betalains stay out, as asked earlier.)
  Files                sequences and bio objects are kept in .sed files; .cdxml leaves them out (a message says so) -
                       expand sequences first to keep them as structures there.

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
  biopolymer.py    sequences: monomers, FASTA / three-letter / HELM reading and writing, the molecule, analysis
  biodraw.py       drawing of sequences and BioDraw objects; bio_dialog.py the sequence and plasmid editors
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
