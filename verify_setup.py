"""
Verification script to check if the environment is properly set up.
"""
import os
import sys

def check_dependencies():
    """Check if all required dependencies are installed."""
    print("Checking dependencies...")
    
    missing = []
    
    try:
        import pandas
        print("[OK] pandas installed")
    except ImportError:
        missing.append("pandas")
        print("[MISSING] pandas not installed")
    
    try:
        import numpy
        print("[OK] numpy installed")
    except ImportError:
        missing.append("numpy")
        print("[MISSING] numpy not installed")
    
    try:
        from rdkit import Chem
        print("[OK] rdkit installed")
    except ImportError:
        missing.append("rdkit")
        print("[MISSING] rdkit not installed")
    
    return missing

def check_files():
    """Check if all required files exist."""
    print("\nChecking required files...")
    
    missing = []
    
    # Check library file
    library_path = os.path.join("libraries", "HELMCoreLibrary.json")
    if os.path.exists(library_path):
        print(f"[OK] {library_path} exists")
    else:
        missing.append(library_path)
        print(f"[MISSING] {library_path} not found")
    
    # Check test files
    test_files = [
        os.path.join("test-sets", "HELM_LINEAR.csv"),
        os.path.join("test-sets", "HELM_cyclic.csv")
    ]
    
    for test_file in test_files:
        if os.path.exists(test_file):
            print(f"[OK] {test_file} exists")
        else:
            missing.append(test_file)
            print(f"[MISSING] {test_file} not found")
    
    # Check logic files
    logic_files = [
        os.path.join("logics", "main.py"),
        os.path.join("logics", "pipeline.py"),
        os.path.join("logics", "monomer_library.py"),
        os.path.join("logics", "target_processor.py")
    ]
    
    for logic_file in logic_files:
        if os.path.exists(logic_file):
            print(f"[OK] {logic_file} exists")
        else:
            missing.append(logic_file)
            print(f"[MISSING] {logic_file} not found")
    
    return missing

def main():
    print("=" * 50)
    print("Mol-to-HELM Environment Verification")
    print("=" * 50)
    print()
    
    missing_deps = check_dependencies()
    missing_files = check_files()
    
    print("\n" + "=" * 50)
    if not missing_deps and not missing_files:
        print("[SUCCESS] Setup verification passed!")
        print("=" * 50)
        print("\nYou can now run the tests:")
        print("  cd logics")
        if sys.platform == "win32":
            print("  py main.py")
        else:
            print("  python main.py")
        return 0
    else:
        print("[FAILED] Setup verification failed!")
        print("=" * 50)
        
        if missing_deps:
            print("\nMissing dependencies:")
            for dep in missing_deps:
                print(f"  - {dep}")
            print("\nRun setup script to install dependencies:")
            if sys.platform == "win32":
                print("  setup.bat")
            else:
                print("  ./setup.sh")
        
        if missing_files:
            print("\nMissing files:")
            for file in missing_files:
                print(f"  - {file}")
        
        return 1

if __name__ == "__main__":
    exit(main())

