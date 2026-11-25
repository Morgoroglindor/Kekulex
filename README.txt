Kekulex (formerly Structure editor) - draw chemical structures for figures and slides, and read off their formula, exact mass and the
m/z of the usual LC-MS adducts. A sibling of LC viewer and Spectra viewer.

Run from source:   pip install -r requirements.txt  ;  python kekulex.py [drawing.sed | file.mol | file.sdf | file.smi]
Windows installer: double-click BUILD_INSTALLER.bat; it creates installer_output\Kekulex_Setup.exe
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
                     hash, 6 / 7 / v / u rings, a phenyl, z ethynyl, Shift+K tBu, 9 two methyls; + / - charge;
                     t (or Enter) types its label; with the Text tool, just point at it and type;
                     the element's key again takes one H away each time, as in ChemDraw (NH2 -> NH -> N -> NH2,
                     OH -> O -> OH, c on a vertex: CH3 -> CH2 -> CH -> C -> CH3); an atom left short of H gets the
                     red box
                   a bond - 1 2 3 order, l / r / c double-bond side, w wedge, h hash, 4-8 / v fused ring,
                     a fused benzene, z fused cyclopentadiene
                   nothing - j benzene tool, x bond, Shift+X chain, e arrow, t text, Space select
                   Help > Keyboard shortcuts lists all of them.
  Select           click / Shift+click / drag a box; drag to move; drop one atom on another to merge; arrows nudge.
                   Drag the round handle at the selection's top right corner to turn it freely (Shift: 15° steps).
  SMILES box       paste a SMILES, Enter. Ctrl+V pastes SMILES or MOL text. File > Open reads .mol, .sdf, .smi.
  Clean up         Ctrl+Shift+K or Ctrl+Alt+K lays the selected structures out again (RDKit), keeping stereo.
                   Ctrl+Shift+K stays Kekulex's while its window is in front, even if another program (Notion)
                   registered it system-wide (winkeys.py).
  Usual drawings   a structure from a name or SMILES, or cleaned up, follows the usual orientation of common skeletons
                   (steroids, triterpenes, betalains, purines, indoles, quinolines, flavonoids, penicillins, sugars
                   ...: layout2d.py). A betalain by name comes in its usual form (ring NH, iminium N+).
  3D tab           select a molecule, Structure > Send to 3D (Ctrl+3, also in the right-click Molecule menu): an
                   RDKit + MMFF94 model with its conformers. Left drag rotates, right drag pans, the wheel zooms;
                   click 2 / 3 / 4 atoms for a distance / angle / dihedral. It is a copy made when sent: it does
                   not follow later edits and is not saved in the file. Choose MMFF94, MMFF94s or UFF, search
                   conformers again, minimize one or all. Export for MD programs: PDB, MOL2 (Gasteiger charges),
                   GROMACS .gro, AMBER .inpcrd, LAMMPS data, SDF, XYZ. They hold coordinates and bonds; force-field
                   parameters come from ACPYPE / antechamber, CGenFF or LigParGen, which read the MOL2.
                   Surfaces: solvent accessible (1.4 Å probe) or van der Waals, solid, translucent, mesh or dots,
                   coloured by element, Gasteiger partial charge or Crippen hydrophobicity. Also H-bonds inside
                   the molecule, depth fading, perspective, axes, background colour, and the conformer's
                   properties: surface areas, volume, ovality, principal moments, radius of gyration, H-bond
                   donors / acceptors.
  Autosave         unsaved work is copied every minute to %APPDATA%\Kekulex\autosave (autosave.py); after a crash
                   or power cut, Kekulex offers the drawing back when it next starts. Saving or closing on purpose
                   removes the copy.
  Optional parts   the installer asks whether to include NMR prediction (about 70 MB) and the 3D
                   tab; run it again to add or remove them. A part left out disappears from the menus.

ChemDraw-like features (0.3.0)
  Page view            the drawing sits on an A4 / Letter page (File > Document settings); View > Whole page.
  Nicknames            atom labels OMe, CO2H, Boc, Ph, NO2, OTs, TMS ... stay short on the drawing (turned round,
                       MeO2C, when bonded on the right) and count in full for formula / SMILES / MOL; so do condensed
                       formulas typed on an atom (condensed.py): SO3K, SO3Na, CO2Na (as salts), CH2CH3, COCH3, C(O)CH3,
                       N(CH3)2, CH(CH3)2, (CH2)3CH3, PO3H2, P(O)(OH)2, SO2NH2, SO2Cl, OCF3, CH=CH2, NMe3+ ...;
                       all in lower case works too (so3k, ch2cl, nh2, boc: two letters are one symbol only for
                       Cl Br Na Li Mg Si and nicknames, so co is C + O, never cobalt).
                       Structure > Expand label draws them out. Labels and captions are typed in place.
  Bond styles          wedge, hashed wedge, bold, hashed, dashed (dashed second line on a double bond), wavy.
  Arrows               reaction, equilibrium, retrosynthetic, resonance, dashed, curved, fish-hook (click a curved
                       arrow again with its tool to bend it the other way).
  Radicals and ions     the charge palette (the ⊕ button's arrow) has +, −, radical •, diradical ••, radical cation
                       •+ and radical anion •−: click an atom (again: taken away). Or point at an atom and press . for
                       a radical (•), again for a carbene / diradical (••), again for none; or right-click > Atom >
                       Radical, Charge > Radical cation •+ / anion •−, Hydrogens (0-3 or automatic). Typed labels say it too: CH2• , CH3+, O•+, C••, and CH2 on a carbon keeps
                       exactly two H. Each unpaired electron takes one H away; carbocations (CH3+), acylium ions and
                       carbenes are allowed, so MS fragments and reactive intermediates can be drawn. Formula, masses
                       and SMILES follow; odd-electron ions are listed as [M]•+. Kept in .sed, MOL (RAD) and .cdxml.
  Isotopes             type 13C, 13CH3, 15NH2, [18O], ¹³C or D / T on an atom (or right-click > Atom > Isotope: the
                       element's stable isotopes and the usual radiolabels 14C, 3H, 18F, 125I ...). The mass number is
                       drawn up before the symbol (H3¹³C), D and T as letters. Formula (C¹³CH6O, CH3D), exact mass,
                       isotope pattern, SMILES ([13CH3]) and InChIKey follow; kept in .sed, MOL (ISO) and .cdxml.
  Missing charges      Structure > Add missing charges (also right-click > Molecule): the formal charge the bonds call
                       for - N or P with four bonds +, O or S with three +, B with four −, O typed with no H on one bond
                       −, NH3 typed on one bond +. Carbon is left alone (a cation or an anion is your call); Check
                       structure names the atoms that need one.
  Partial charges      View > Show partial charges: Gasteiger charges (with the atom's H) as discs behind the atoms,
                       red for δ− (electron rich, nucleophilic), blue for δ+ (electrophilic), deeper for stronger;
                       point at an atom for the value. On screen only (not printed or exported).
  MS fragments         right-click a bond > Bond > Cleave (MS fragments): each way of breaking it, named by its two
                       pieces and their m/z - a radical cation (M•+) gives A+ + B• or A• + B+ (alpha cleavage:
                       acetone M•+ -> CH3C≡O+ m/z 43.0178 + CH3•), an even-electron ion ([M+H]+) the inductive
                       cleavage A+ + B (neutral), a neutral molecule A• + B• or A+ + B−. The pieces are drawn under an
                       arrow with live captions (formula, m/z or neutral mass); break them again for a fragmentation
                       tree. Ring bonds are not offered (one piece). ms.py.
  Fragments tab        Structure > MS fragments, or the Fragments tab of the formula panel: the fragments the structure
                       can logically give, from the textbook rules, for ESI+ ([M+H]+), ESI− ([M−H]−) or EI (M•+) -
                       EI: alpha cleavage next to N, O, S and C=O (iminium, oxonium, acylium ions), benzylic / allylic
                       cleavage (C7H7+ 91), inductive cleavage at C-halogen, McLafferty rearrangement, retro-Diels-
                       Alder; ESI: C-O / C-N / C-S heterolysis (acylium and b-type ions, glycoside Y ions, y-type
                       amines, carboxylates, quinate ...), retro-Diels-Alder of flavonoids (1,3A ions); both: neutral
                       losses (H2O, NH3, CO, CO2, CH3OH, HCOOH, CH3•, ketene, C4H8, SO3, H3PO4, NO2, HX, C-glycoside
                       90 / 120). Steps 2 or 3: the fragments lose small molecules again. Each row: m/z, ion, what is
                       lost, the rule and a likelihood (likely / possible / minor - a plausibility, not an intensity).
                       Picking a row shows the atoms kept (green), lost (red) and the bonds that break; Draw (or
                       double-click) puts the ion under the structure after an arrow with the loss and a live caption;
                       Copy gives the list as text. A structure drawn as an ion is used as it is. The same fragments
                       label the peaks in Annotate spectrum. fragments.py, fragments_tab.py.
  Match m/z           Structure > Match m/z (Ctrl+M): a measured m/z within so many ppm against every drawn structure's
                       ions and adducts (charged fragments as themselves; double-click selects it) and against
                       elemental compositions of the ion (elements C H N O S by default, any of P F Cl Br I Si B Na K
                       and caps such as N4 can be added), with ppm, RDBE and even / odd electron; the usual ratio rules
                       (H/C, N/C, O/C ...) leave out impossible formulas.
  Annotate spectrum    Structure > Annotate mass spectrum: open a measured spectrum (two-column CSV / TXT, JCAMP-DX,
                       NIST MSP, MGF, mzML - its scan with the most ions first, any scan from the list) or paste one;
                       profile data is reduced to its peak tops. Each peak gets the best of the drawn structures' ions
                       within the tolerance (ppm, or Da for unit-mass EI spectra): adducts, M•+, a drawn fragment ion
                       as itself, the rule-based fragments of the Fragments tab and, with Bond cleavages, each piece of breaking one bond ([PH+H]+, [PH]•+, P+;
                       negative [PH−H]−). Stick plot with the matched peaks labelled (C2H3O+), table to copy
                       (double-click selects the structure or piece), Insert into drawing puts the annotated spectrum
                       on the page as a picture. msfile.py, spectrum_dialog.py.
  Mechanisms           curved and fish-hook arrows drawn from an atom or bond to an atom or bond stay on them when the
                       structure moves (drag the arrow alone to move its ends). Right-click > Atom > Lone pairs: none,
                       all it has (counted from its electrons) or 1-3, drawn as dot pairs in the free directions.
                       Structure > Reaction > Run the curved arrows (Ctrl+Shift+R, or right-click a curved arrow):
                       the structure the arrows make is drawn to the right of the step with a reaction arrow - a full
                       arrow moves a pair, a fish-hook one electron, from a lone pair or a bond to an atom (a lone pair
                       there, or a new bond) or to a bond (one order more); charges, radicals and bond orders follow
                       by counting electrons and every atom keeps its H. The selected curved arrows run, else those on
                       the selected structure, else the rightmost step's. Wrong arrows are named instead (from an atom
                       with no lone pair, a fish-hook with no partner, carbon with five bonds after the step: one more
                       arrow should take a pair away), and Check structure lists them too.
  Find                 Edit > Find (Ctrl+F): selects captions holding the text, atoms with that label or element, or a
                       substructure written as SMILES / SMARTS (C(=O)O, c1ccccc1, [OX2H]), nicknames included.
  Reactions            Structure > Reaction > Map atoms: atoms of the reactants and products get matching numbers (:1,
                       :2 ..., from the largest common substructures); Clear atom map; Check balance: what the products
                       have more or fewer than the reactants, and a by-product that would balance it (H2O, HCl, CO2 ...).
                       Only the arrow's own row and step count; reagents over / under the arrow are left out.
  Forms                Structure > Forms > Tautomers / Resonance forms: the other forms of the selected structure drawn
                       to its right, laid out like it, with equilibrium or resonance arrows (RDKit).
  Check structure      also lists stereocentres with no wedge, wedges that start at an atom that is not a stereocentre,
                       and atoms with too many wedges.
  Snap guides          dragging a selection lines its edges or centre up with another object's (pink guides); hold
                       Alt to move freely.
  Reaction schemes     right-click a reaction arrow > Reaction conditions: reagents over the arrow, solvent / temperature
                       / time under it (formula style, grouped with the arrow, which grows to fit). Object > Arrange as
                       reaction scheme: the selected structures and arrows in one row, centred on one line, evenly spaced.
  Batch export         File > Batch export: an SDF or a SMILES list (SMILES, then a name) to one grid sheet or one
                       picture per structure (PNG / SVG / PDF), named, in this drawing's style. batch.py.
  Shapes               brackets, parentheses, braces, rectangle, rounded rectangle, oval, line.
  Pictures             Ctrl+V pastes a copied picture (from a browser, a screenshot, Word, a copied PNG / JPEG file) at
                       its own size, File > Insert picture reads one, and image files dropped on the window land where
                       they fall. A picture moves, scales (Object > Scale) and flips with the drawing, sits behind the
                       structures, and goes into exports (PNG, SVG, PDF) and .sed files (not .cdxml). Copied SMILES text
                       still pastes as a structure even when a picture came with it. pictures.py.
  Readable drawings    nothing is drawn on anything else: an atom's H goes to the side that is free (OH beside HO),
                       charges and radical dots keep off other bonds and labels, atom numbers, R/S and map numbers look
                       round their atom for a clear spot, lone pairs sit in free directions just off the letter, labels
                       read the other way when that is clear (TBSO, BocHN, Me3N+), analysis captions and MS fragments
                       move down past what is already there, and the layout turns crowded end groups apart (two OH of a
                       diphosphate). overlaps.py finds what still overlaps; the tests keep it at none for common
                       molecules. Very crowded cages with every number switched on can still touch.
  Polymer brackets     right-click a bracket pair > Bracket / shape > Polymer repeat count: n (or 3, x ...) at the
                       foot of the right bracket. The atoms inside count that many times in the formula and masses;
                       a letter gives C2H6(C2H4)n, with the masses for n = 1. Written to .cdxml as ChemDraw's
                       bracketed group (SRU, with its label) and read back.
  Captions             size, bold, italic, formula style (H2SO4 with subscripts), left / centre / right, colour.
  Groups               Object > Group (Ctrl+G) / Ungroup (Ctrl+Shift+G): grouped structures, captions, arrows and shapes
                       select, move, turn and are analysed as one (a salt or mixture: one formula, mass and properties;
                       the Structure list shows "Group"). Captions pasted under a structure join its group and move with
                       it. Kept in .sed files, copies and .cdxml (as ChemDraw groups).
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
                       (structures left of the arrow and above / below it are reactants, right of it products;
                       grouped structures - a salt, a reagent with its acid - are one compound, also when grouped
                       after the grid was made). A card per side: the formulas on top (the limiting reactant
                       marked), one row per quantity, typed values in blue, calculated ones in ink, units in grey.
                       Double-click a white cell to type into it on the drawing, with or without units (250 mg,
                       1.5 mmol, 2 mL, 0.5 M, 1.2 g/mL, 95 %; a bare number is g, mol, mL ...): the grid follows each
                       keystroke, Enter keeps it, Esc puts the old value back, Tab goes to the next cell, an emptied
                       cell is calculated again. Double-click a reactant's formula to make it the limiting one;
                       double-click a title or row name (or right-click > Table > Edit values) for the grid's window,
                       where rows can be hidden. Values: equivalents, sample mass, %weight, molarity, density, volume,
                       moles; for products equivalents, %completion, measured mass, purity. Moles, masses, volumes,
                       expected mass and %yield follow and are redrawn whenever the reaction changes.
  View menu            show reaction interpretation (reactants red, products magenta, above / below the arrow
                       cyan / blue), rulers in cm, whole page.
  Templates            rings and ring systems, heterocycles, functional groups, reagents, solvents.
  ChemDraw files       .cdxml opens and saves (structures, nicknames, bond styles, captions, arrows, shapes, colours,
                       document style); the stoichiometry grid and its typed values open too (not written back to
                       .cdxml yet: keep .sed for that). .cdx opens the structures (through RDKit). Not read: TLC
                       plates, orbitals, spectra.
  Analysis panel       Analysis tab: formula, exact mass, m/z with isotopes, MW, elemental analysis (each with
                       Paste, which puts it under the structure as a caption; "Paste all under the structure" puts them
                       all in one caption), adduct m/z, isotope pattern of the chosen ion. Properties tab: logP, MR,
                       tPSA, H-bond donors / acceptors, rotatable bonds, rings, heavy atoms, fraction sp3,
                       stereocentres, rule of five (RDKit); tick rows and press Paste (none ticked: all of them), or
                       double-click one row. Pasted values (and Show analysis captions) are live: they are rewritten
                       whenever the structure changes, until you type over the caption or ungroup it.
  NMR prediction       Structure > Predict 1H / 13C NMR (also right-click > Analysis): CSP5 (MIT, machine learning
                       on measured spectra, solvent aware: CDCl3, DMSO-d6, CD3OD, D2O ...). A window with the
                       simulated spectrum (mouse as in Spectra viewer: drag a box to zoom, wheel zooms ppm, Ctrl+wheel
                       the height, right drag pans, double-click / F whole spectrum, click picks a signal; the splitter
                       gives it more room; saved as PNG / SVG / PDF), the signals (shift +/- uncertainty, first-order
                       multiplicity with expected J: vicinal 7 Hz, ring CH-CH by Karplus on an RDKit 3D conformer,
                       aromatic ortho 8 / meta 2, alkene cis 10-11 / trans 16 ...), spectrometer 300-800 MHz, the
                       experimental-section line with J to copy, Label atoms, Paste text. OH / NH shifts are flagged as
                       unreliable, and left out in D2O / CD3OD. An optional part of the installer; the CSP5 models run on numpy (nmr_engine.py), without PyTorch.
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

Drawing elements, orbitals, fragmentation, TLC (as in ChemDraw's manual; the palettes open as grids)
  Drawing elements     circle (drag out from the centre), oval, rounded rectangle, rectangle - each plain, dashed,
                       shaded, filled or with a shadow; lines plain, dashed, wavy, bold; arcs 270 / 180 / 120 / 90
                       degrees, solid or dashed (drag left to right: convex, right to left: concave). Right-click >
                       Style changes it later. Kept in .cdxml (RectangleType / OvalType / LineType words); arcs stay
                       in .sed.
  Orbitals (Shift+G)   s, σ, single lobe, p, hybrid sp3, dxy, dz2 - plain, shaded or filled (the + phase). Click an
                       atom to put one on it, or drag from the atom to set size and direction.
  Fragmentation        Mass fragmentation: drag across bonds - the formula and exact mass of each side appear under
                       the structure, as homolytic cleavage (radicals; Alt keeps the cut line; a ring bond alone
                       gives one formula). Retrosynthesis: the product, ⇒, then the pieces (filled up with H).
                       Synthesis: the pieces, →, then the product. Not predictions: you choose the bonds.
  TLC plate (g)        drag a plate (the width sets the lanes) or click. Drag a spot to set its Rf (shown while
                       dragging), Shift+drag to resize it, Ctrl+drag to copy it to another lane; Eraser on a spot
                       deletes it. Right-click: add / duplicate / delete lane, show origin, solvent front, border,
                       side ticks, Rf of all or one spot, spot style (filled, shaded, hollow) and colour, Set Rf.
                       Double-click: type the Rf values, one lane per line. Kept in .sed files.

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
  kekulex.py   window + entry point: tools, menus, files, clipboard, export, style dialog
  canvas.py        the drawing area: tools, mouse, keys, selection, undo / redo, clean up
  model.py         the drawing (atoms, bonds, arrows, text, shapes, nicknames) and its geometry; MOL; RDKit in
  cdxml.py         ChemDraw .cdxml reading and writing
  biopolymer.py    sequences: monomers, FASTA / three-letter / HELM reading and writing, the molecule, analysis
  biodraw.py       drawing of sequences and BioDraw objects; bio_dialog.py the sequence, plasmid and TLC editors
  elements.py      drawing elements (styles, arcs) and orbitals; fragment.py mass fragmentation and (retro)synthesis
  tlc.py           TLC plates: layout, drawing, spots and Rf
  templates.py     the Templates menu (SMILES)
  stoich.py        stoichiometry grid: reaction columns, unit parsing, calculations; stoich_dialog.py its editor
  pubchem.py       name <-> structure through PubChem (online, on request only)
  render.py        drawing with QPainter in points (screen and every export the same); PNG / SVG / PDF
  chem.py          formula, masses, adducts, SMILES, InChIKey (RDKit); reading SMILES / MOL / SDF
  info_panel.py    the formula and masses panel
  theme.py, applog.py   look (light / dark) and log file + crash dialog - copied from Spectra viewer
  icons.py         tool icons drawn at start-up
  version.py       the version number (About box and installer)
  tools\make_icon.py   draws kekulex.ico

The .sed file is JSON: atoms (x, y in bond lengths, y down; element; charge), bonds (atoms, order 1-3, stereo
wedge / hash with the narrow end at atom a, double-bond side), arrows, texts and the drawing style.
Log: %APPDATA%\Kekulex\kekulex.log
