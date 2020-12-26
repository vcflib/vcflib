% PFST(1) pFst (vcflib) | pFst (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

pFst

# SYNOPSIS

usage: pFst --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1 --type PL

# DESCRIPTION

pFst is a probabilistic approach for detecting differences in allele frequencies between two populations.



# OPTIONS

```


Output : 3 columns :     
     1. seqid            
     2. position         
     3. pFst probability 

required: t,target     -- argument: a zero based comma separated list of target individuals corresponding to VCF columns       
required: b,background -- argument: a zero based comma separated list of background individuals corresponding to VCF columns   
required: f,file       -- argument: a properly formatted VCF.                                                                  
required: y,type       -- argument: genotype likelihood format ; genotypes: GP, GL or PL; pooled: PO                           
optional: d,deltaaf    -- argument: skip sites where the difference in allele frequencies is less than deltaaf, default is zero
optional: r,region     -- argument: a tabix compliant genomic range : seqid or seqid:start-end                                 
optional: c,counts     -- switch  : use genotype counts rather than genotype likelihoods to estimate parameters, default false 

Type: statistics

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[pFst.cpp](https://github.com/vcflib/vcflib/blob/master/src/pFst.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
