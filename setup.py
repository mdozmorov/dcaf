from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import glob

scripts = ["script/genome_runner"]

setup(
        name="dcaf",
        version="0.0.1-HEAD",
        author="Mikhail Dozmorov; Cory Giles",
        author_email="dozmorovm@omrf.org; mail@corygil.es",
        description="General Wren Lab scripts and utilities.",

        packages=["dcaf"],
        scripts=scripts,
        cmdclass= {"build_ext": build_ext},
        ext_modules = [
            Extension(
                "dcaf.bbi",
                include_dirs=["src"],
                sources=["dcaf/bbi.pyx", "src/bbi.cpp", "src/mm.cpp"],
                extra_compile_args=["-std=c++0x"],
                libraries=["z"],
                language="c++")
            ]
)
