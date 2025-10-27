@echo off
echo ====================================
echo Setting up Mol-to-HELM Environment
echo ====================================
echo.

echo Installing required packages...
py -m pip install --upgrade pip
py -m pip install rdkit>=2023.3.1 pandas>=2.0.0 numpy>=1.24.0

echo.
echo ====================================
echo Installation Complete!
echo ====================================
echo.
echo To run the tests, execute:
echo   cd logics
echo   py main.py
echo.
pause

