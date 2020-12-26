% PERMUTESMOOTH(1) permuteSmooth (vcflib) | permuteSmooth (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

permuteSmooth

# SYNOPSIS

usage: permuteSmooth -s wcFst.smooth.txt -f wcFst.txt -n 5 -s 1

# DESCRIPTION

./build/permuteSmooth: invalid option -- 'h' FATAL: no file was provided

# OPTIONS

```





     permuteSmooth is a method for adding empirical p-values  smoothed wcFst scores.

Required:
      file:     f   -- argument: original wcFst data     
      smoothed: s   -- argument: smoothed wcFst data     
      format:   y   -- argument: [swcFst, segwcFst]      
Optional:
      number:   n   -- argument: the number of permutations to run for each value [1000]
      success:  u   -- argument: stop permutations after 's' successes [1]
      success:  x   -- argument: number of threads [1]

OUTPUT: permuteSmooth will append three additional columns:
        1. The number of successes                            
        2. The number of trials                               
        3. The empirical p-value                              


```

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[permuteSmooth.cpp](https://github.com/vcflib/vcflib/blob/master/src/permuteSmooth.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
