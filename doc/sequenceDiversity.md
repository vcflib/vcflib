% SEQUENCEDIVERSITY(1) sequenceDiversity (vcflib) | sequenceDiversity (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

**sequenceDiversity**

# SYNOPSIS

**sequenceDiversity** --target 0,1,2,3,4,5,6,7 --file my.vcf

# DESCRIPTION

The **sequenceDiversity** program calculates two popular metrics of haplotype diversity: pi and extended haplotype homozygoisty (eHH). Pi is calculated using the Nei and Li 1979 formulation. eHH a convenient way to think about haplotype diversity. When eHH = 0 all haplotypes in the window are unique and when eHH = 1 all haplotypes in the window are identical.



# OPTIONS

```


Output : 5 columns:
         1.  seqid
         2.  start of window
         3.  end of window  
         4.  pi             
         5.  eHH            


required: t,target     -- argument: a zero base comma separated list of target individuals corresponding to VCF columns        
required: f,file       -- argument: a properly formatted phased VCF file                                                       
required: y,type       -- argument: type of genotype likelihood: PL, GL or GP                                                  
optional: a,af         -- sites less than af  are filtered out; default is 0                                          
optional: r,region     -- argument: a tabix compliant region : "seqid:0-100" or "seqid"                                    
optional: w,window     -- argument: the number of SNPs per window; default is 20                                               


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[sequenceDiversity.cpp](https://github.com/vcflib/vcflib/blob/master/src/sequenceDiversity.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
