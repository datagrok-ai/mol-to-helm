#!/usr/bin/env python3
"""
Aggregation script that combines all Python files from the logics folder
into a single Python file, removing internal imports and keeping only external library imports.
"""

import os
import ast
import re
from typing import List, Set, Dict


class ImportAnalyzer(ast.NodeVisitor):
    """Analyzes Python AST to extract import information"""
    
    def __init__(self):
        self.imports = []
        self.from_imports = []
        
    def visit_Import(self, node):
        """Handle 'import module' statements"""
        for alias in node.names:
            self.imports.append({
                'type': 'import',
                'module': alias.name,
                'alias': alias.asname,
                'lineno': node.lineno
            })
        self.generic_visit(node)
    
    def visit_ImportFrom(self, node):
        """Handle 'from module import ...' statements"""
        module = node.module if node.module else ''
        for alias in node.names:
            self.from_imports.append({
                'type': 'from_import',
                'module': module,
                'name': alias.name,
                'alias': alias.asname,
                'lineno': node.lineno
            })
        self.generic_visit(node)


def get_python_files(directory: str) -> List[str]:
    """Get all Python files in the directory"""
    python_files = []
    skipped_files = {'main.py', 'debug_aze_3abz.py', 'debug_hArg.py'}
    for file in os.listdir(directory):
        if file.endswith('.py') and not file.startswith('__') and not file in skipped_files:
            python_files.append(os.path.join(directory, file))
    return sorted(python_files)


def analyze_file_imports(file_path: str) -> tuple:
    """Analyze imports in a Python file"""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    try:
        tree = ast.parse(content)
        analyzer = ImportAnalyzer()
        analyzer.visit(tree)
        return analyzer.imports, analyzer.from_imports, content
    except SyntaxError as e:
        print(f"Warning: Could not parse {file_path}: {e}")
        return [], [], content


def get_local_modules(logics_dir: str) -> Set[str]:
    """Get set of local module names (without .py extension)"""
    local_modules = set()
    for file in os.listdir(logics_dir):
        if file.endswith('.py') and not file.startswith('__'):
            local_modules.add(file[:-3])  # Remove .py extension
    return local_modules


def is_external_import(module_name: str, local_modules: Set[str]) -> bool:
    """Check if an import is external (not a local module)"""
    if not module_name:
        return True
    
    # Split module name to handle nested imports like 'package.module'
    base_module = module_name.split('.')[0]
    return base_module not in local_modules


def remove_internal_imports(content: str, local_modules: Set[str]) -> str:
    """Remove internal imports from file content"""
    lines = content.split('\n')
    filtered_lines = []
    
    for line in lines:
        stripped_line = line.strip()
        
        # Skip empty lines and comments at the start
        if not stripped_line or stripped_line.startswith('#'):
            if not any(filtered_lines):  # Only skip at the beginning
                continue
            else:
                filtered_lines.append(line)
                continue
        
        # Check for import statements
        is_internal_import = False
        
        # Handle 'from module import ...' statements
        if stripped_line.startswith('from ') and ' import ' in stripped_line:
            match = re.match(r'from\s+(\S+)\s+import', stripped_line)
            if match:
                module_name = match.group(1)
                if not is_external_import(module_name, local_modules):
                    is_internal_import = True
        
        # Handle 'import module' statements
        elif stripped_line.startswith('import '):
            match = re.match(r'import\s+(\S+)', stripped_line)
            if match:
                module_name = match.group(1)
                if not is_external_import(module_name, local_modules):
                    is_internal_import = True
        
        # Keep the line if it's not an internal import
        if not is_internal_import:
            filtered_lines.append(line)
    
    return '\n'.join(filtered_lines)


def collect_external_imports(logics_dir: str, local_modules: Set[str]) -> Set[str]:
    """Collect all unique external imports from all files"""
    external_imports = set()
    
    python_files = get_python_files(logics_dir)
    
    for file_path in python_files:
        imports, from_imports, _ = analyze_file_imports(file_path)
        
        # Process regular imports
        for imp in imports:
            if is_external_import(imp['module'], local_modules):
                if imp['alias']:
                    external_imports.add(f"import {imp['module']} as {imp['alias']}")
                else:
                    external_imports.add(f"import {imp['module']}")
        
        # Process from imports
        for imp in from_imports:
            if is_external_import(imp['module'], local_modules):
                if imp['alias']:
                    external_imports.add(f"from {imp['module']} import {imp['name']} as {imp['alias']}")
                else:
                    external_imports.add(f"from {imp['module']} import {imp['name']}")
    
    return external_imports


def aggregate_files(logics_dir: str, output_file: str):
    """Main function to aggregate all files"""
    if not os.path.exists(logics_dir):
        print(f"Error: Directory {logics_dir} does not exist")
        return
    
    # Get local modules
    local_modules = get_local_modules(logics_dir)
    print(f"Found local modules: {sorted(local_modules)}")
    
    # Collect all external imports
    external_imports = collect_external_imports(logics_dir, local_modules)
    
    # Get all Python files
    python_files = get_python_files(logics_dir)
    
    if not python_files:
        print(f"No Python files found in {logics_dir}")
        return
    
    # Start building the aggregated content
    aggregated_content = []
    
    
    # Add Datagrok part annotations ***
    aggregated_content.append('#language: python')
    aggregated_content.append('#name: molToHelmConverterPy')
    aggregated_content.append('#description: Converts molecules to HELM notation based on monomer library')
    aggregated_content.append('#input: dataframe moleculesDataframe')
    aggregated_content.append('#input: column moleculesColumn {semType: Molecule}')
    aggregated_content.append('#input: string libraryJSON')
    aggregated_content.append('#output: dataframe result_helm {action:join(moleculesDataframe)} [Sequences, in HELM format]')
    aggregated_content.append('molListToProcess = moleculesDataframe[moleculesColumn].tolist()')
    aggregated_content.append('import pandas as pd')
    aggregated_content.append('import numpy as np')
    # END Datagrok part annotations
    
    
    # Add header comment
    aggregated_content.append('"""')
    aggregated_content.append('Aggregated file combining all modules from the logics folder.')
    aggregated_content.append('Generated automatically - do not edit manually.')
    aggregated_content.append('"""')
    aggregated_content.append('')
    
    # Add all external imports
    if external_imports:
        aggregated_content.append('# External library imports')
        for imp in sorted(external_imports):
            aggregated_content.append(imp)
        aggregated_content.append('')
    
    # Process each file
    for file_path in python_files:
        file_name = os.path.basename(file_path)
        print(f"Processing {file_name}...")
        
        _, _, content = analyze_file_imports(file_path)
        
        # Remove internal imports and clean up
        cleaned_content = remove_internal_imports(content, local_modules)
        
        # Remove leading/trailing whitespace but preserve internal structure
        cleaned_content = cleaned_content.strip()
        
        if cleaned_content:
            # Add file separator comment
            aggregated_content.append(f'# ============================================================================')
            aggregated_content.append(f'# Content from: {file_name}')
            aggregated_content.append(f'# ============================================================================')
            aggregated_content.append('')
            aggregated_content.append(cleaned_content)
            aggregated_content.append('')
    
    
    # add datagrok func call result_helm
    aggregated_content.append('res_helm_list = convert_molecules_batch(molListToProcess, library_json=libraryJSON)')
    # create dataframe result_helm
    aggregated_content.append('result_helm = pd.DataFrame(map(lambda x: x[1], res_helm_list), columns=["regenerated sequences"])')
    
    
    # Write the aggregated file
    final_content = '\n'.join(aggregated_content)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(final_content)
    
    print(f"Aggregated file created: {output_file}")
    print(f"Combined {len(python_files)} files")
    print(f"External imports: {len(external_imports)}")


if __name__ == "__main__":
    # Define paths
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    logics_dir = os.path.join(project_root, 'logics')
    output_file = os.path.join(current_dir, 'mol_to_helm_aggregated.py')
    
    print("Starting aggregation process...")
    print(f"Source directory: {logics_dir}")
    print(f"Output file: {output_file}")
    print()
    
    aggregate_files(logics_dir, output_file)
    print("\nAggregation complete!")