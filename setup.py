#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import subprocess
import os
import shutil
from pathlib import Path

class CustomBuildExt(build_ext):
    def run(self):
        jess_dir = Path('src/isopeptor/resources/jess')
        
        # Ensure the directory exists
        if not jess_dir.exists():
            raise RuntimeError(f"Source directory {jess_dir} does not exist")

        # Clean any existing object files
        try:
            subprocess.check_call(['make', 'clean'], cwd=str(jess_dir))
        except subprocess.CalledProcessError:
            print("Warning: Clean failed, continuing anyway")

        # Compile Jess
        try:
            print("Building Jess")
            subprocess.check_call(['make'], cwd=str(jess_dir))
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Error building Jess: {e}")

        # Create destination directory
        dest_dir = Path(self.build_lib) / 'isopeptor' / 'resources' / 'jess'
        dest_dir.mkdir(parents=True, exist_ok=True)

        # Copy the executable
        exec_path = jess_dir / 'jess'
        if not exec_path.exists():
            raise RuntimeError("Jess executable not found after compilation")
        print("Jess built succesfully")
        shutil.copy2(exec_path, dest_dir)
        os.chmod(dest_dir / 'jess', 0o755)
        
        build_ext.run(self)

setup(
    name='isopeptor',
    version='0.0.1',
    description='Package for isopeptide bond prediction and analysis based on Jess',
    author='Francesco Costa',
    author_email='fcosta@ebi.ac.uk',
    
    # Package structure
    package_dir={'': 'src'},
    packages=find_packages(where='src'),

    package_data={
    "isopeptor": ["resources/data/*"],
    },
    
    # Include non-Python files
    include_package_data=True,
    
    # Package dependencies
    install_requires=[
        'numpy',
        'biotite',
        'joblib',
        'scikit-learn'
    ],

    entry_points={
        'console_scripts': [
            'isopeptor=isopeptor.cli:main',
        ],
    },
    
    # Custom build command
    cmdclass={
        'build_ext': CustomBuildExt,
    },
    
    python_requires='>=3.8',
    
    # Additional metadata
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    
    project_urls={
        'Source': 'https://github.com/FranceCosta/isopeptor',
        'Bug Reports': 'https://github.com/FranceCosta/isopeptor/issues',
    },
)