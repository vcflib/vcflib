% PVST(1) pVst (vcflib) | pVst (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**pVst**

# SYNOPSIS

**pVst** --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --type CN

# DESCRIPTION

**pVst** calculates vst, a measure of CNV stratification.



# OPTIONS

```




The statistic Vst is used to test the difference in copy numbers at
each SV between two groups: Vst = (Vt-Vs)/Vt, where Vt is the overall
variance of copy number and Vs the average variance within
populations.

Output : 4 columns :     
     1. seqid            
     2. position         
     3. end              
     3. vst              
     4. probability      

required: t,target     -- argument: a zero based comma separated list of target individuals corresponding to VCF columns       
required: b,background -- argument: a zero based comma separated list of background individuals corresponding to VCF columns   
required: f,file       -- argument: a properly formatted VCF.                                                                  
required: y,type       -- argument: the genotype field with the copy number: e.g. CN|CNF                           
optional: r,region     -- argument: a tabix compliant genomic range : seqid or seqid:start-end                                 
optional: x,cpu        -- argument: number of CPUs [1] 
optional: n,per        -- argument: number of permutations [1000] 

Type: statistics



```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[pVst.cpp](https://github.com/vcflib/vcflib/blob/master/src/pVst.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
