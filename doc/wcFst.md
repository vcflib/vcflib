% WCFST(1) wcFst (vcflib) | wcFst (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**wcFst**

# SYNOPSIS

**wcFst** --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1 --type PL

# DESCRIPTION

**wcFst** is Weir & Cockerham's Fst for two populations. Negative values are VALID, they are sites which can be treated as zero Fst. For more information see Evolution, Vol. 38 N. 6 Nov 1984. Specifically **wcFst** uses equations 1,2,3,4.



# OPTIONS

```


Output : 3 columns :                 
     1. seqid                        
     2. position                     
     3. target allele frequency      
     4. background allele frequency  
     5. **wcFst**                        

required: t,target     -- argument: a zero based comma separated list of target individuals corrisponding to VCF columns        
required: b,background -- argument: a zero based comma separated list of background individuals corrisponding to VCF columns    
required: f,file       -- argument: proper formatted VCF                                                                        
required, y,type       -- argument: genotype likelihood format; genotype : GT,GL,PL,GP                                             
optional: r,region     -- argument: a tabix compliant genomic range: seqid or seqid:start-end                                   
optional: d,deltaaf    -- argument: skip sites where the difference in allele frequencies is less than deltaaf, default is zero 

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

[wcFst.cpp](https://github.com/vcflib/vcflib/blob/master/src/wcFst.cpp)

# LICENSE

Copyright 2011-2025 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2025 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
