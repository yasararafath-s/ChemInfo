@echo off
echo ============================================
echo    ChemInfo Dashboard - Installation
echo ============================================
echo.

echo [1/3] Creating Python virtual environment...
python -m venv venv
if errorlevel 1 (
    echo ERROR: Python not found. Please install Python 3.10+ from https://www.python.org/downloads/
    pause
    exit /b 1
)

echo.
echo [2/3] Activating virtual environment...
call venv\Scripts\activate.bat

echo.
echo [3/3] Installing dependencies (this may take a few minutes)...
pip install -r requirements.txt

echo.
echo ============================================
echo    Installation Complete!
echo    Run START.bat to launch the dashboard
echo ============================================
echo.
pause
