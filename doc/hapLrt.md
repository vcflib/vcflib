% HAPLRT(1) hapLrt (vcflib) | hapLrt (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

hapLrt

# SYNOPSIS



# DESCRIPTION

 HapLRT is a likelihood ratio test for haplotype lengths. The lengths are modeled with an exponential distribution. The sign denotes if the target has longer haplotypes (1) or the background (-1). 

# OPTIONS

```


Output : 4 columns :                             
     1. seqid                                    
     2. position                                 
     3. mean target haplotype length             
     4. mean background haplotype length         
     5. p-value from LRT                         
     6. sign                                     

hapLRT  --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --type GP --file my.vcf                                     

required: t,target     -- argument: a zero base comma separated list of target individuals corrisponding to VCF columns        
required: b,background -- argument: a zero base comma separated list of background individuals corrisponding to VCF columns    
required: f,file       -- argument: a properly formatted phased VCF file                                                       
required: y,type       -- argument: type of genotype likelihood: PL, GL or GP                                                  
optional: r,region     -- argument: a genomice range to calculate hapLrt on in the format : "seqid:start-end" or "seqid" 

```

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[hapLrt.cpp](https://github.com/vcflib/vcflib/blob/master/src/hapLrt.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
