@echo off
rem Double-click to open a folder picker, or drag a results folder onto this file.
rem Extra options can be passed from a terminal:  run_viewer.bat D:\data\brain --no-meshes
rem
rem Needs nothing installed beforehand except an internet connection on the first run:
rem uv (the Python launcher) is installed for the current user if missing, then uv
rem downloads its own Python 3.12 and the exact package versions from
rem view_in_neuroglancer.py.lock. No Visual C++ is needed.
setlocal
set "UV=uv"
where uv >nul 2>nul
if not errorlevel 1 goto run
if exist "%USERPROFILE%\.local\bin\uv.exe" (
    set "UV=%USERPROFILE%\.local\bin\uv.exe"
    goto run
)

echo uv is not installed. Installing it for this user (no administrator rights needed)...
powershell -NoProfile -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
if not exist "%USERPROFILE%\.local\bin\uv.exe" (
    echo.
    echo Could not install uv automatically ^(no internet, or blocked by a firewall/antivirus^).
    echo Install it by hand from https://docs.astral.sh/uv/getting-started/installation/
    echo then close this window and try again.
    pause
    exit /b 1
)
set "UV=%USERPROFILE%\.local\bin\uv.exe"

:run
rem --python 3.12: ignore whatever Python is installed (too new or too old) and use
rem   the version every compiled dependency has prebuilt wheels for (no compiler needed;
rem   the only source-only package, atomicwrites, is pure Python).
"%UV%" run --python 3.12 "%~dp0view_in_neuroglancer.py" %*
if errorlevel 1 pause
