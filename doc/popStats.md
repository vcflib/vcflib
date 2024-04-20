% POPSTATS(1) popStats (vcflib) | popStats (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**popStats**

# SYNOPSIS

popStat --type PL --target 0,1,2,3,4,5,6,7 --file my.vcf

# DESCRIPTION

General population genetic statistics for each SNP



# OPTIONS

```



  Calculates basic population statistics at bi-allelic sites. The allele frequency is the number of non-reference alleles divided by the total number of alleles.  The expected hetrozygosity is 2*p*q, where p is the non-reference allele frequency and q is 1-p.  The observed heterozgosity is the fraction of 0/1 genotypes out of all genotypes.  The inbreeding coefficient, Fis, is the relative heterozygosity of each individual vs. compared to the target group. 

Output : 9 columns :                 
     1. seqid                        
     2. position                     
     3. target allele frequency      
     4. expected heterozygosity      
     5. observed heterozygosity      
     6. number of hets               
     7. number of homozygous ref     
     8. number of homozygous alt     
     9. target Fis                   
required: t,target     -- a zero based comma separated list of target individuals corresponding to VCF columns        
required: f,file       -- proper formatted VCF                                                                        
required, y,type       -- genotype likelihood format; genotype : GL,PL,GP                                             
optional, r,region     -- a tabix compliant region : chr1:1-1000 or chr1                                              

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

[popStats.cpp](https://github.com/vcflib/vcflib/blob/master/src/popStats.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
