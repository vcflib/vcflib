% VCFLD(1) vcfld (vcflib) | vcfld (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfld**

# SYNOPSIS

**vcfld** --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf -e -d -r

# DESCRIPTION

Compute LD



# OPTIONS

```


required: t,target     -- argument: a zero base comma separated list of target individuals corresponding to VCF columns        
required: b,background -- argument: a zero base comma separated list of background individuals corresponding to VCF columns    
required: f,file       -- argument: a properly formatted phased VCF file                                                       
required: y,type       -- argument: type of genotype likelihood: PL, GL or GP                                                  
optional: w,window     -- argument: window size to average LD; default is 1000                                                 
optional: e,external   -- switch: population to calculate LD expectation; default is target                                    
optional: d,derived    -- switch: which haplotype to count "00" vs "11"; default "00",                                   


Type: transformation


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

[vcfld.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfld.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
