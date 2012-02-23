from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(ext_modules=[Extension(
                   "pyvcflib.cvcflib",                 # name of extension
                   ["cvcflib.pyx", "Variant.cpp"], #  our Cython source
                   language="c++")],  # causes Cython to create C++ source
                   include_dirs=["../", "."],
                   library_dirs=["../", "."],
                   cmdclass={'build_ext': build_ext})