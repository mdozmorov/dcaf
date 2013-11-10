import sys
import os
import pkgutil
import subprocess

from setuptools import setup
from setuptools.command.test import test as TestCommand
from distutils.extension import Extension
from pip.req import parse_requirements

import dcaf

NAME = "dcaf"
VERSION = dcaf.__version__
AUTHOR = dcaf.__author__
AUTHOR_EMAIL = dcaf.__author_email__
DESCRIPTION = "Utilities for genome analysis, expression analysis, and text mining used by the Wren Lab at the Oklahoma Medical Research Foundation."
CLASSIFIERS = [
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
]
LICENSE = "AGPLv3+"
REQUIREMENTS = [str(item.req) for item in parse_requirements("requirements.txt")]

# If we are in a git repository, write the current version (based on
# git tag) to version.py, so that it can be loaded from a source
# package.

with open("dcaf/version.py", "w") as h:
    h.write("__version__ = '%s'\n" % dcaf.__version__)
 
cmdclass = {}
extensions = []

try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize

    extensions.extend([
        Extension(
            "dcaf.io.bbi",
            include_dirs=["include"],
            sources=["dcaf/io/bbi.pyx", "src/bbi/bbi.cpp"],
            extra_compile_args=["-std=c++0x"],
            libraries=["z"],
            language="c++")
    ])
    cmdclass["build_ext"] = build_ext
except ImportError:
    pass

try:
    from sphinx.setup_command import BuildDoc
    cmdclass["doc"] = BuildDoc
except ImportError:
    pass

# Find registered entry points

entry_points = {}

try:
    # Import all submodules so that the entry points will be properly registered
    # for wrapper script autogeneration.

    for loader, module_name, is_pkg in pkgutil.walk_packages(dcaf.__path__):
        module_name = "dcaf." + module_name
        loader.find_module(module_name).load_module(module_name)

    entry_points={
        "console_scripts": 
        ["%s = %s:%s" % (script_name, fn.__module__, fn.__name__)
         for script_name, fn in dcaf._entry_points.items()]
    }
except ImportError:
    pass

   
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

cmdclass["test"] = Test

# Create list of data files with paths relative to the base dcaf directory
package_data = []
for root, dirs, files in os.walk("data"):
    print(files)
    for file in files:
        package_data.append(os.path.relpath(os.path.join(root, file), "dcaf"))

setup(
    # Metadata
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    classifiers=CLASSIFIERS,
    license=LICENSE,

    # Modules, data, and extensions to be installed 
    packages=["dcaf", "dcaf.io"],
    package_data={"dcaf": package_data},
    install_requires=REQUIREMENTS,
    tests_require=REQUIREMENTS + ["pytest"],
    extras_require={"doc": REQUIREMENTS},
    ext_modules=extensions,

    # Executable scripts
    scripts=scripts,
    entry_points=entry_points,
    
    # setup.py entry points
    cmdclass=cmdclass
)
