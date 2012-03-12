from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import glob

dependencies = ['smithwaterman/SmithWatermanGotoh.h',
                'tabixpp/bgzf.h',
                'tabixpp/tabix.h',
                'tabixpp/tabix.hpp',
                'multichoose/multichoose.h']
print dependencies
setup(
    ext_modules=[
      Extension("cvcflib",
                 sources=["split.cpp", 
                          "smithwaterman/SmithWatermanGotoh.cpp",
                          "multichoose/multichoose.cpp",
                          "tabixpp/tabix.cpp",
                          "tabixpp/bgzf.c",
                          "Variant.cpp",
                          "pyvcflib/cvcflib.pyx"],
                 libraries=["stdc++", 'z', 'm', 'tabix'],
                 language="c++",
                 include_dirs=["./", "smithwaterman/", "tabixpp/", "multichoose/"],
                 depends = dependencies,
                 library_dirs=["tabixpp/"]),
      ],
    cmdclass = {'build_ext': build_ext},
)                   

setup(
        name="pyvcflib",
        version="0.1.0",
        packages=['pyvcflib'],
        author="arq5x",
        description='Wrapper around vcflib',
        url="none",
        package_data = {'pyvcflib':["*.pyx",
                                      "*.pxi",
                                      "*.pxd",
                                      "*.cpp"]
                       },
        py_modules=["pyvcflib/__init__"],
        author_email="arq5x@virginia.edu",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']

    )