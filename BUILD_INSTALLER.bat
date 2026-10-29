@echo off
rem ==========================================================================================
rem  Builds StructureEditor_Setup.exe, a normal Windows installer for the Structure editor.
rem  Double-click this file. It checks for Python and Inno Setup and installs them if missing.
rem  Needs an internet connection the first time; takes about 10-20 minutes.
rem  Result: installer_output\StructureEditor_Setup.exe  (the people who use the programs need neither
rem  Python nor this folder: they just run StructureEditor_Setup.exe)
rem ==========================================================================================
setlocal EnableExtensions
title Building the Structure editor installer
cd /d "%~dp0"
echo.
echo  ============================================================
echo   Structure editor : installer builder
echo  ============================================================
echo.

rem ---------------------------------------------------------------- 1. Python
call :findpython
if not defined PYEXE (
    echo [1/5] Python 3.10 - 3.13 was not found on this computer. Installing Python 3.13 ...
    call :installpython
    call :findpython
)
if not defined PYEXE (
    echo.
    echo  Python could not be installed automatically.
    echo  Please install Python 3.13 from https://www.python.org/downloads/
    echo  ^(tick "Add python.exe to PATH" in the first screen^) and run this file again.
    goto :fail
)
echo [1/5] Python found: %PYEXE%

rem ---------------------------------------------------------------- 2. libraries
echo [2/5] Installing the libraries the programs need ...
if not exist ".venv\Scripts\python.exe" ( "%PYEXE%" -m venv .venv || goto :fail )
call ".venv\Scripts\activate.bat"
python -m pip install --upgrade pip >nul
python -m pip install -r requirements-dev.txt || goto :fail
for /f "delims=" %%V in ('python -c "from version import __version__;print(__version__)"') do set "APPVER=%%V"
if not defined APPVER goto :fail
echo       version %APPVER%

rem ---------------------------------------------------------------- 3. programs
echo [3/5] Building the programs (this is the slow step) ...
python -m PyInstaller --noconfirm --clean --windowed --onedir --name StructureEditor --icon structureeditor.ico --paths . ^
  --add-data "structureeditor.ico;." --hidden-import theme --hidden-import icons --hidden-import stoich_dialog --hidden-import pubchem --collect-all rdkit --collect-all openclatura ^
  --exclude-module tkinter --exclude-module matplotlib --exclude-module pandas --exclude-module PyQt5 --exclude-module PySide2 --exclude-module PySide6 ^
  --exclude-module IPython --exclude-module ipywidgets --exclude-module notebook --exclude-module pytest ^
  structure_editor.py || goto :fail
if not exist "dist\StructureEditor\StructureEditor.exe" goto :fail

rem ---------------------------------------------------------------- 4. Inno Setup (makes the installer)
call :findiscc
if not defined ISCC (
    echo [4/5] Inno Setup was not found. Installing it ...
    call :installinno
    call :findiscc
)
if not defined ISCC (
    echo.
    echo  Inno Setup could not be installed automatically.
    echo  Please install it from https://jrsoftware.org/isdl.php and run this file again.
    goto :fail
)
echo [4/5] Inno Setup found: %ISCC%

rem ---------------------------------------------------------------- 5. installer
echo [5/5] Creating the installer ...
"%ISCC%" /DAppVersion=%APPVER% installer.iss
if errorlevel 1 (
    echo Retrying without the Portuguese language file ...
    "%ISCC%" /DAppVersion=%APPVER% /DNOPT installer.iss || goto :fail
)
if not exist "installer_output\StructureEditor_Setup.exe" goto :fail
echo.
echo  ============================================================
echo   Done.  installer_output\StructureEditor_Setup.exe
echo   Double-click it to install, or send it to colleagues.
echo  ============================================================
start "" explorer "%~dp0installer_output"
pause
exit /b 0

:fail
echo.
echo  The build did not finish. Scroll up to read the last messages.
pause
exit /b 1

rem ---------------------------------------------------------------- helpers
:findpython
set "PYEXE="
for %%V in (3.13 3.12 3.11 3.10) do (
    if not defined PYEXE (
        for /f "delims=" %%P in ('py -%%V -c "import sys;print(sys.executable)" 2^>nul') do set "PYEXE=%%P"
    )
)
if not defined PYEXE (
    python -c "import sys;sys.exit(0 if (3,10)<=sys.version_info[:2]<=(3,13) else 1)" >nul 2>&1
    if not errorlevel 1 (
        for /f "delims=" %%P in ('python -c "import sys;print(sys.executable)"') do set "PYEXE=%%P"
    )
)
if not defined PYEXE (
    for %%V in (313 312 311 310) do (
        if not defined PYEXE if exist "%LOCALAPPDATA%\Programs\Python\Python%%V\python.exe" set "PYEXE=%LOCALAPPDATA%\Programs\Python\Python%%V\python.exe"
    )
)
exit /b 0

:installpython
where winget >nul 2>&1
if not errorlevel 1 (
    winget install -e --id Python.Python.3.13 --scope user --silent --accept-package-agreements --accept-source-agreements
    exit /b 0
)
echo winget is not available: downloading the Python installer from python.org ...
powershell -NoProfile -ExecutionPolicy Bypass -Command "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12; Invoke-WebRequest 'https://www.python.org/ftp/python/3.13.1/python-3.13.1-amd64.exe' -OutFile '%TEMP%\python-installer.exe'"
if exist "%TEMP%\python-installer.exe" "%TEMP%\python-installer.exe" /quiet InstallAllUsers=0 PrependPath=1 Include_launcher=1 Include_test=0
exit /b 0

:findiscc
set "ISCC="
for %%P in ("%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe" "%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe" "%ProgramFiles%\Inno Setup 6\ISCC.exe") do (
    if exist %%P if not defined ISCC set "ISCC=%%~P"
)
exit /b 0

:installinno
where winget >nul 2>&1
if not errorlevel 1 (
    winget install -e --id JRSoftware.InnoSetup --silent --accept-package-agreements --accept-source-agreements
    exit /b 0
)
echo winget is not available: downloading Inno Setup from jrsoftware.org ...
powershell -NoProfile -ExecutionPolicy Bypass -Command "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12; Invoke-WebRequest 'https://jrsoftware.org/download.php/is.exe' -OutFile '%TEMP%\innosetup.exe'"
if exist "%TEMP%\innosetup.exe" "%TEMP%\innosetup.exe" /VERYSILENT /SUPPRESSMSGBOXES /NORESTART /CURRENTUSER
exit /b 0
