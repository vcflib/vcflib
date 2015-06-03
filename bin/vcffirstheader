#!/usr/bin/env python

import sys

header=True
for line in sys.stdin:
    if line.startswith('##'):
        if header:
            print line.strip()
        continue
    elif line.startswith('#'):
        if header:
           print line.strip()
           header=False
        continue
    print line.strip()
