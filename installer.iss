; Inno Setup script - builds StructureEditor_Setup.exe (a normal Windows installer wizard).
; Compiled by BUILD_INSTALLER.bat after PyInstaller has produced dist\StructureEditor\.
#define AppName "Structure editor"
; the version comes from version.py (BUILD_INSTALLER.bat passes /DAppVersion=...); this is only the fallback
#ifndef AppVersion
#define AppVersion "0.5.0"
#endif
#define AppExe "StructureEditor.exe"

[Setup]
AppId={{58EABA69-73C3-45F8-B7F6-A2061787A8E8}
AppName={#AppName}
AppVersion={#AppVersion}
DefaultDirName={autopf}\StructureEditor
DefaultGroupName=Structure editor
UninstallDisplayIcon={app}\{#AppExe}
SetupIconFile=structureeditor.ico
OutputDir=installer_output
OutputBaseFilename=StructureEditor_Setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64compatible
; installs for the current user without asking for an administrator password (the user can still choose "all users")
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
DisableProgramGroupPage=yes
ChangesAssociations=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
#ifndef NOPT
Name: "ptbr"; MessagesFile: "compiler:Languages\BrazilianPortuguese.isl"
#endif

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"
Name: "assoc"; Description: "Open Structure editor drawings (.sed) when double-clicked"; GroupDescription: "File types:"

[Files]
Source: "dist\StructureEditor\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion
Source: "README.txt"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\Structure editor"; Filename: "{app}\{#AppExe}"
Name: "{autodesktop}\Structure editor"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Registry]
Root: HKA; Subkey: "Software\Classes\.sed\OpenWithProgids"; ValueType: string; ValueName: "StructureEditor.sed"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\StructureEditor.sed"; ValueType: string; ValueData: "Structure editor drawing"; Flags: uninsdeletekey; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\StructureEditor.sed\DefaultIcon"; ValueType: string; ValueData: "{app}\{#AppExe},0"; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\StructureEditor.sed\shell\open\command"; ValueType: string; ValueData: """{app}\{#AppExe}"" ""%1"""; Tasks: assoc

[Run]
Filename: "{app}\{#AppExe}"; Description: "Start Structure editor now"; Flags: nowait postinstall skipifsilent
