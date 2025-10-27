#!/bin/bash

echo "===================================="
echo "Setting up Mol-to-HELM Environment"
echo "===================================="
echo ""

echo "Installing required packages..."
pip install --upgrade pip
pip install -r requirements.txt

echo ""
echo "===================================="
echo "Installation Complete!"
echo "===================================="
echo ""
echo "To run the tests, execute:"
echo "  cd logics"
echo "  python main.py"
echo ""

