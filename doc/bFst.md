% BFST(1) bFst (vcflib) | bFst (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**bFst**

# SYNOPSIS

**bFst** --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1

# DESCRIPTION

**bFst** is a Bayesian approach to Fst. Importantly **bFst** accounts for genotype uncertainty in the model using genotype likelihoods. For a more detailed description see: `A Bayesian approach to inferring population structure from dominant markers' by Holsinger et al. Molecular Ecology Vol 11, issue 7 2002. The likelihood function has been modified to use genotype likelihoods provided by variant callers. There are five free parameters estimated in the model: each subpopulation's allele frequency and Fis (fixation index, within each subpopulation), a free parameter for the total population's allele frequency, and Fst.



# OPTIONS

```


Output : 11 columns :                          
     1.  Seqid                                     
     2.  Position				     
     3.  Observed allele frequency in target.	     
     4.  Estimated allele frequency in target.     
     5.  Observed allele frequency in background.  
     6.  Estimated allele frequency in background. 
     7.  Observed allele frequency combined. 	     
     8.  Estimated allele frequency in combined.   
     9.  ML estimate of Fst (mean)		     
     10. Lower bound of the 95% credible interval  
     11. Upper bound of the 95% credible interval  

required: t,target     -- a zero bases comma separated list of target individuals corrisponding to VCF columns
required: b,background -- a zero bases comma separated list of background individuals corrisponding to VCF columns
required: f,file a     -- a proper formatted VCF file.  the FORMAT field MUST contain "PL"
required: d,deltaaf    -- skip sites were the difference in allele frequency is less than deltaaf


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

[bFst.cpp](https://github.com/vcflib/vcflib/blob/master/src/bFst.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
