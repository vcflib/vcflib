#!/usr/bin/env python
"""
Some simple examples of working with data loaded into numpy arrays from a VCF 
file.
 
"""

import sys

try:
    from vcfnp import variants, info, calldata, view2d
except:
    print 'please build the vcflib Cython extension: python setup.py build_ext --inplace'
    raise

try:
    import numpy as np
except:
    print 'please install numpy'
    raise

try:
    import matplotlib.pyplot as plt
except:
    print 'please install matplotlib'
    raise


if __name__ == '__main__':
    
    if len(sys.argv) < 2:
        print 'usage: python vcfnp_example.py [VCFFILE]'
        sys.exit(1)
        
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
    c = calldata(filename).view(np.recarray)
    c = view2d(c)
    count_phased = np.count_nonzero(c.is_phased)
    count_variant = np.count_nonzero(np.any(c.genotype > 0, axis=2)) 
    count_missing = np.count_nonzero(~c.is_called)
    print 'calls (phased, variant, missing): %s (%s, %s, %s)' % (c.flatten().size, count_phased, count_variant, count_missing)
    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    ax.hist(c.GQ.flatten())
    ax.set_title('GQ histogram')
    ax.set_xlabel('GQ')
    plt.show()
    
    
    
    

