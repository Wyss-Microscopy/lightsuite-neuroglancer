@echo off
rem Double-click to open a folder picker, or drag a results folder onto this file.
rem Extra options can be passed from a terminal:  run_viewer.bat D:\data\brain --no-meshes
setlocal
where uv >nul 2>nul
if errorlevel 1 (
    echo uv was not found. Install it once with:
    echo   powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
    echo then close this window and try again.
    pause
    exit /b 1
)
uv run "%~dp0view_in_neuroglancer.py" %*
if errorlevel 1 pause
