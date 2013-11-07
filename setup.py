import os
import pkgutil
import subprocess

from setuptools import setup
from setuptools.command.test import test as TestCommand
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
from sphinx.setup_command import BuildDoc

import dcaf

# If we are in a git repository, write the current version (based on
# git tag) to version.py, so that it can be loaded from a source
# package.

with open("dcaf/version.py", "w") as h:
    h.write("__version__ = '%s'\n" % dcaf.__version__)
    
# Import all submodules so that the entry points will be properly registered
# for wrapper script autogeneration.

def import_all():
    for loader, module_name, is_pkg in pkgutil.walk_packages(dcaf.__path__):
        module_name = "dcaf." + module_name
        loader.find_module(module_name).load_module(module_name)

import_all()

# Also install any scripts that are in the top level of the script/ directory

scripts = [os.path.abspath("script/" + p) \
           for p in os.listdir("script/") if os.path.isfile(p)]

# Run py.test tests from setup.py

class Test(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True
    
    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        raise SystemExit(errno)

# Declare extensions

extensions = [
    Extension(
        "dcaf.bbi",
        include_dirs=["include"],
        sources=["dcaf/bbi.pyx", "src/bbi/bbi.cpp"],
        extra_compile_args=["-std=c++0x"],
        libraries=["z"],
        language="c++")
]

setup(
    # Metadata
    name="dcaf",
    version=dcaf.__version__,
    author=dcaf.__author__,
    author_email=dcaf.__author_email__,
    description="Utilities for genome analysis, expression analysis, and text mining used by the Wren Lab at the Oklahoma Medical Research Foundation.",
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Programming Language :: Unix Shell",
        "Programming Language :: SQL",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    license="AGPLv3+",

    # Modules, data, and extensions to be installed 
    packages=["dcaf"],
    package_data={"dcaf": "data"},
    tests_require=["pytest"],
    extras_require={"doc": "sphinx"},
    ext_modules = extensions,

    # Executable scripts
    scripts=scripts,
    entry_points={
        "console_scripts": 
        ["%s = %s:%s" % (script_name, fn.__module__, fn.__name__)
         for script_name, fn in dcaf._entry_points.items()]
    },
    
    # setup.py entry points
    cmdclass= {
        "build_ext": build_ext,
        "test": Test,
        "doc": BuildDoc
    }
)
