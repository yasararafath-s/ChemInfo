@echo off
echo ============================================
echo    ChemInfo Dashboard - Starting...
echo ============================================
echo.

call venv\Scripts\activate.bat

echo Starting Streamlit server...
echo.
echo    Dashboard URL:  http://localhost:8501
echo    Press Ctrl+C to stop
echo.
echo ============================================
echo.

streamlit run app.py --server.port 8501 --server.headless true --browser.gatherUsageStats false
