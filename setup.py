import glob
import subprocess

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

scripts = ["script/genome_runner"]

def version():
    return subprocess.check_output(["git", "describe", "--tags"]).decode("ascii").strip()

extensions = [
        Extension(
            "dcaf.bbi",
            include_dirs=["src"],
            sources=["dcaf/bbi.pyx", "src/bbi.cpp", "src/mm.cpp"],
            extra_compile_args=["-std=c++0x"],
            libraries=["z"],
            language="c++")
        ]

setup(
        name="dcaf",
        version=version(),
        author="Cory Giles, Mikhail Dozmorov",
        author_email="mail@corygil.es; dozmorovm@omrf.org",
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
        ]
        license="AGPLv3+",

        packages=["dcaf"],
        scripts=scripts,
        cmdclass= {"build_ext": build_ext},
        ext_modules = extensions
)
