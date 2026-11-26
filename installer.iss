; Inno Setup script - builds Kekulex_Setup.exe (a normal Windows installer wizard).
; Compiled by BUILD_INSTALLER.bat after PyInstaller has produced dist\Kekulex\.
#define AppName "Kekulex"
; the version comes from version.py (BUILD_INSTALLER.bat passes /DAppVersion=...); this is only the fallback
#ifndef AppVersion
#define AppVersion "0.11.0"
#endif
#define AppExe "Kekulex.exe"

[Setup]
AppId={{58EABA69-73C3-45F8-B7F6-A2061787A8E8}
AppName={#AppName}
AppVersion={#AppVersion}
DefaultDirName={autopf}\Kekulex
DefaultGroupName=Kekulex
UninstallDisplayIcon={app}\{#AppExe}
SetupIconFile=kekulex.ico
OutputDir=installer_output
OutputBaseFilename=Kekulex_Setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64compatible
; installs for the current user without asking for an administrator password (the user can still choose "all users")
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
DisableProgramGroupPage=yes
ChangesAssociations=yes

; The folders only NMR prediction uses (checked by starting the program without them)
#define NmrDirs "nmr_models,matplotlib,PIL"
#define NmrBuilt DirExists(SourcePath + "\dist\Kekulex\_internal\nmr_models")

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
#ifndef NOPT
Name: "ptbr"; MessagesFile: "compiler:Languages\BrazilianPortuguese.isl"
#endif

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"
Name: "assoc"; Description: "Open Kekulex drawings (.sed) when double-clicked"; GroupDescription: "File types:"

[Types]
Name: "full"; Description: "Everything"
Name: "compact"; Description: "Drawing only (no NMR prediction, no 3D)"
Name: "custom"; Description: "Choose"; Flags: iscustom

[Components]
Name: "main"; Description: "Kekulex"; Types: full compact custom; Flags: fixed
#if NmrBuilt
Name: "nmr"; Description: "NMR prediction (the CSP5 models and the spectrum plot)"; Types: full custom
#endif
Name: "view3d"; Description: "3D tab (very small)"; Types: full custom

; Choosing fewer parts when updating an installation removes what is no longer wanted
[InstallDelete]
; the program's earlier name
Type: files; Name: "{app}\StructureEditor.exe"
Type: files; Name: "{autoprograms}\Structure editor.lnk"
Type: files; Name: "{autodesktop}\Structure editor.lnk"
; earlier versions carried PyTorch and the rest of CSP5's dependencies
Type: filesandordirs; Name: "{app}\_internal\torch"
Type: filesandordirs; Name: "{app}\_internal\csp5"
Type: filesandordirs; Name: "{app}\_internal\pyarrow"
Type: filesandordirs; Name: "{app}\_internal\pyarrow.libs"
Type: filesandordirs; Name: "{app}\_internal\scipy"
Type: filesandordirs; Name: "{app}\_internal\scipy.libs"
Type: filesandordirs; Name: "{app}\_internal\sklearn"
Type: filesandordirs; Name: "{app}\_internal\pandas"
Type: filesandordirs; Name: "{app}\_internal\pandas.libs"
#if NmrBuilt
#define NmrRest NmrDirs + ","
#define NmrDir ""
#sub NmrDelete
Type: filesandordirs; Name: "{app}\_internal\{#NmrDir}"; Components: not nmr
#endsub
#define NmrCut() Local[0] = Pos(",", NmrRest), NmrDir = Copy(NmrRest, 1, Local[0] - 1), NmrRest = Copy(NmrRest, Local[0] + 1), 0
#for {0; Len(NmrRest) > 0; 0} (NmrCut(), NmrDelete)
Type: files; Name: "{app}\_internal\features\nmr.txt"; Components: not nmr
#endif
Type: files; Name: "{app}\_internal\features\3d.txt"; Components: not view3d

[Files]
#define NmrRest NmrDirs + ","
#define NmrExcludes ""
#define ExCut() Local[0] = Pos(",", NmrRest), NmrExcludes = NmrExcludes + (NmrExcludes == "" ? "" : ",") + "\_internal\" + Copy(NmrRest, 1, Local[0] - 1), NmrRest = Copy(NmrRest, Local[0] + 1), 0
#for {0; Len(NmrRest) > 0; 0} ExCut()
Source: "dist\Kekulex\*"; Excludes: "{#NmrExcludes}"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion
#if NmrBuilt
#define NmrRest NmrDirs + ","
#sub NmrFiles
Source: "dist\Kekulex\_internal\{#NmrDir}\*"; DestDir: "{app}\_internal\{#NmrDir}"; Components: nmr; Flags: recursesubdirs createallsubdirs ignoreversion skipifsourcedoesntexist
#endsub
#for {0; Len(NmrRest) > 0; 0} (NmrCut(), NmrFiles)
Source: "features_markers\nmr.txt"; DestDir: "{app}\_internal\features"; Components: nmr; Flags: ignoreversion
#endif
Source: "features_markers\3d.txt"; DestDir: "{app}\_internal\features"; Components: view3d; Flags: ignoreversion
Source: "README.txt"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\Kekulex"; Filename: "{app}\{#AppExe}"
Name: "{autodesktop}\Kekulex"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Registry]
; the file type under the program's earlier name
Root: HKA; Subkey: "Software\Classes\.sed\OpenWithProgids"; ValueType: none; ValueName: "StructureEditor.sed"; Flags: deletevalue
Root: HKA; Subkey: "Software\Classes\StructureEditor.sed"; ValueType: none; Flags: deletekey
Root: HKA; Subkey: "Software\Classes\.sed\OpenWithProgids"; ValueType: string; ValueName: "Kekulex.sed"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\Kekulex.sed"; ValueType: string; ValueData: "Kekulex drawing"; Flags: uninsdeletekey; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\Kekulex.sed\DefaultIcon"; ValueType: string; ValueData: "{app}\{#AppExe},0"; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\Kekulex.sed\shell\open\command"; ValueType: string; ValueData: """{app}\{#AppExe}"" ""%1"""; Tasks: assoc

[Run]
Filename: "{app}\{#AppExe}"; Description: "Start Kekulex now"; Flags: nowait postinstall skipifsilent
