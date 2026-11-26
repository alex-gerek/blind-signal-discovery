@echo off
setlocal EnableExtensions
title BSD Local Measurement Workstation

set "ROOT=%~dp0"
if "%ROOT:~-1%"=="\" set "ROOT=%ROOT:~0,-1%"
cd /d "%ROOT%"

set "VENV_PY=%ROOT%\.venv\Scripts\python.exe"
set "VENV_ST=%ROOT%\.venv\Scripts\streamlit.exe"
set "APP_FILE=%ROOT%\app\bsd_workstation.py"
set "PORT=8501"
set "URL=http://127.0.0.1:%PORT%"

if not exist "%ROOT%\.tmp" mkdir "%ROOT%\.tmp"
set "TEMP=%ROOT%\.tmp"
set "TMP=%ROOT%\.tmp"
set "PYTHONPATH=%ROOT%\src;%PYTHONPATH%"
set "STREAMLIT_BROWSER_GATHER_USAGE_STATS=false"

echo ========================================================================
echo   Blind Signal Discovery -- Local Measurement Workstation
echo   Working Directory: %ROOT%
echo ========================================================================
echo.

if not exist "%VENV_PY%" goto :no_venv
if not exist "%VENV_ST%" goto :no_streamlit
if not exist "%APP_FILE%" goto :no_app

"%VENV_PY%" -c "import streamlit, numpy, scipy, matplotlib, bsd" >nul 2>&1
if errorlevel 1 goto :bad_imports

REM Check whether a BSD Workstation Streamlit instance is already running on PORT
"%VENV_PY%" -c "import urllib.request, sys; r = urllib.request.urlopen('%URL%/_stcore/health', timeout=1.5); sys.exit(0 if r.status == 200 else 1)" >nul 2>&1
if not errorlevel 1 goto :already_running

echo [INFO] Starting BSD Workstation on %URL% ...
echo [INFO] Your default browser will open automatically once ready.
echo [INFO] Press Ctrl+C or close this window to shut down the local server.
echo.

REM Open the default browser automatically as soon as the health endpoint responds
start "" /b powershell -NoProfile -WindowStyle Hidden -Command "for ($i=0; $i -lt 40; $i++) { Start-Sleep -Milliseconds 500; try { $r = Invoke-WebRequest -Uri '%URL%/_stcore/health' -TimeoutSec 2 -UseBasicParsing; if ($r.StatusCode -eq 200) { Start-Process '%URL%'; break } } catch {} }"

"%VENV_ST%" run "%APP_FILE%" --server.address 127.0.0.1 --server.port %PORT% --server.headless true --browser.gatherUsageStats false
set "EXIT_CODE=%ERRORLEVEL%"

if not "%EXIT_CODE%"=="0" (
    echo.
    echo [ERROR] BSD Workstation exited unexpectedly with code %EXIT_CODE%.
    pause
    exit /b %EXIT_CODE%
)
goto :eof

:already_running
echo [INFO] BSD Workstation is already running at %URL%.
echo        Opening existing session in your default browser -- no duplicate server started.
start "" "%URL%"
exit /b 0

:no_venv
echo [ERROR] Local Python virtual environment not found:
echo         %VENV_PY%
echo.
echo Please create the virtual environment in the project folder and install dependencies:
echo   python -m venv .venv
echo   .venv\Scripts\python.exe -m pip install -e ".[gui,ligo,dev]"
echo.
pause
exit /b 1

:no_streamlit
echo [ERROR] Streamlit executable not found in the local virtual environment:
echo         %VENV_ST%
echo.
echo Install the GUI dependency with:
echo   .venv\Scripts\python.exe -m pip install "streamlit>=1.35"
echo.
pause
exit /b 1

:no_app
echo [ERROR] Workstation entry point not found:
echo         %APP_FILE%
echo.
pause
exit /b 1

:bad_imports
echo [ERROR] Required Python packages failed to import: streamlit, numpy, scipy, matplotlib, or bsd.
echo Run the following command to inspect or repair dependencies:
echo   "%VENV_PY%" -c "import streamlit, numpy, scipy, matplotlib, bsd"
echo.
pause
exit /b 1
