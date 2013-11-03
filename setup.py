import glob
import os
import subprocess

from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import dcaf

with open("dcaf/version.py", "w") as h:
    h.write("__version__ = '%s'\n" % dcaf.__version__)
    
scripts = [os.path.abspath("script/" + p) for p in os.listdir("script/") if os.path.isfile(p)]

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

    packages=["dcaf"],

    scripts=scripts,
    entry_points={
        "console_scripts": [
            "genome-runner = dcaf.genome:main"
        ]
    },
    
    cmdclass= {"build_ext": build_ext},
    ext_modules = extensions
)
