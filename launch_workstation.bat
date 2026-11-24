@echo off
setlocal
set "ROOT=%~dp0"
if not exist "%ROOT%.tmp" mkdir "%ROOT%.tmp"
set "TEMP=%ROOT%.tmp"
set "TMP=%ROOT%.tmp"
echo Starting Blind Signal Discovery (BSD) Local Measurement Workstation...
"%ROOT%.venv\Scripts\streamlit.exe" run "%ROOT%app\bsd_workstation.py" --server.headless false
endlocal
