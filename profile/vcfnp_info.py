#!/usr/bin/env python

import sys
import time
import pstats
import cProfile
import timeit


sys.path.append('.')
import vcfnp


def profile():
    a = vcfnp.info(sys.argv[1], count=int(sys.argv[2]))


prof_fn = 'profile/tmp.prof'
cmd = 'profile()'
cProfile.runctx(cmd, globals(), locals(), prof_fn)
s = pstats.Stats(prof_fn)
s.strip_dirs().sort_stats('time').print_stats()
print timeit.repeat(cmd, number=1, repeat=3, setup='from __main__ import profile')


