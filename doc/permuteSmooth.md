% PERMUTESMOOTH(1) permuteSmooth (vcflib) | permuteSmooth (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**permuteSmooth**

# SYNOPSIS

**permuteSmooth** -s wcFst.smooth.txt -f wcFst.txt -n 5 -s 1

# DESCRIPTION

**permuteSmooth** is a method for adding empirical p-values smoothed wcFst scores.



# OPTIONS

```


Required:
      file:     f   -- argument: original wcFst data     
      smoothed: s   -- argument: smoothed wcFst data     
      format:   y   -- argument: [swcFst, segwcFst]      
Optional:
      number:   n   -- argument: the number of permutations to run for each value [1000]
      success:  u   -- argument: stop permutations after 's' successes [1]
      success:  x   -- argument: number of threads [1]

OUTPUT: **permuteSmooth** will append three additional columns:
        1. The number of successes                            
        2. The number of trials                               
        3. The empirical p-value                              


Type: statistics


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[permuteSmooth.cpp](https://github.com/vcflib/vcflib/blob/master/src/permuteSmooth.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
