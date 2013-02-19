from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os
import glob


vcflib_dir = os.getcwd()
smithwaterman_dir = os.path.join(vcflib_dir, 'smithwaterman')
tabixpp_dir = os.path.join(vcflib_dir, 'tabixpp')


def get_vcflib_sources():
    sources = list()
    sources += [os.path.join(vcflib_dir, s) for s in ('Variant.cpp', 'ssw.c', 'ssw_cpp.cpp', 'split.cpp')]
    sources += [os.path.join(smithwaterman_dir, s) for s in ('BandedSmithWaterman.cpp', 'SmithWatermanGotoh.cpp', 'Repeats.cpp', 'disorder.c', 'LeftAlign.cpp', 'IndelAllele.cpp')]    
    sources += [os.path.join(tabixpp_dir, s) for s in ('bedidx.c', 'bgzf.c', 'index.c', 'knetfile.c', 'kstring.c', 'tabix.cpp')]
    return sources


vcflib_extension = Extension('vcflib',
                             sources=['vcflib.pyx'] + get_vcflib_sources(),
                             language='c++',
                             include_dirs=[vcflib_dir, smithwaterman_dir, tabixpp_dir],
                             libraries=['m', 'z'],)
    
    
setup(
    name = 'vcflib',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [vcflib_extension]
)
