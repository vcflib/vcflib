#!/usr/bin/env python

import sys
import time
from vcflib import count_variants


if __name__ == '__main__':
    before = time.time()
    print count_variants(sys.argv[1])
    after = time.time()
    print after-before

