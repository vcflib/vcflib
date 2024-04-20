% PERMUTEGPAT++(1) permuteGPAT++ (vcflib) | permuteGPAT++ (VCF phenotype)
% Erik Garrison and vcflib contributors

# NAME

**permuteGPAT++**

# SYNOPSIS

**permuteGPAT++** -f gpat.txt -n 5 -s 1

# DESCRIPTION

**permuteGPAT++** is a method for adding empirical p-values to a GPAT++ score.



# OPTIONS

```


     Currently **permuteGPAT++** only supports wcFst, but will be extended.    

OUTPUT: **permuteGPAT++** will append three additional columns:
        1. The number of successes                         
        2. The number of trials                            
        3. The empirical p-value                           

file:    f   -- argument: the input file     
number:  n   -- argument: the number of permutations to run for each value [1000]
success: s   -- argument: stop permutations after 's' successes [1]

Type: phenotype

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

[permuteGPAT++.cpp](https://github.com/vcflib/vcflib/blob/master/src/permuteGPAT++.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
