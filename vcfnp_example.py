#!/usr/bin/env python
"""
Some simple examples of working with data loaded into numpy arrays from a VCF 
file.
 
"""

import sys
from vcfnp import variants, info, samples, view2d
try:
    import numpy as np
except:
    print 'please install numpy'
try:
    import matplotlib.pyplot as plt
except:
    print 'please install matplotlib'


if __name__ == '__main__':
    filename = sys.argv[1]
    
    print 'loading variants from %s ...' % filename
    v = variants(filename).view(np.recarray)
    print 'found %s variants (%s SNPs)' % (v.size, np.count_nonzero(v.is_snp))
    print 'QUAL mean (std): %s (%s)' % (np.mean(v.QUAL), np.std(v.QUAL))
    
    print 'loading INFO ...'
    i = info(filename).view(np.recarray)
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.hist(i.DP)
    ax.set_title('DP histogram')
    ax.set_xlabel('DP')
    plt.show()
    
    print 'loading samples ...'
    s = samples(filename).view(np.recarray)
    c = view2d(s)
    missing = np.count_nonzero(~c.is_called)
    hom_ref = np.count_nonzero(np.all(c.genotype == 0, axis=2))
    hom_alt = np.count_nonzero(np.all(c.genotype == 1, axis=2))
    het = np.count_nonzero(np.any(c.genotype == 0, axis=2) & np.any(c.genotype > 0, axis=2))
    print 'calls (missing, hom ref, hom alt, het): %s (%s, %s, %s, %s)' % (c.flatten().size, missing, hom_ref, hom_alt, het)
    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    ax.hist(c.GQ.flatten())
    ax.set_title('GQ histogram')
    ax.set_xlabel('GQ')
    plt.show()
    
    
    
    

