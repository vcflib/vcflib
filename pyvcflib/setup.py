from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    ext_modules=[
      Extension("pyvcflib.cvcflib",                 # name of extension
                 sources=["../Variant.cpp", "cvcflib.pyx"], #  our Cython source
                 libraries=["stdc++", 'z'],
                 language="c++",  # causes Cython to create C++ source
                 include_dirs=["../", "./"],
                 library_dirs=["../"]),
#                 extra_objects=["../Variant.o"]),
      ],
    cmdclass = {'build_ext': build_ext},
)                   
