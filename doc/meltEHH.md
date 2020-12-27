% MELTEHH(1) meltEHH (vcflib) | meltEHH (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**meltEHH**

# SYNOPSIS

**meltEHH** --target 0,1,2,3,4,5,6,7 --pos 10 --file my.phased.vcf \ --region chr1:1-1000 > STDOUT 2> STDERR

# DESCRIPTION

**meltEHH** provides the data to plot extended haplotype homozygosity (EHH) curves.



# OPTIONS

```


Output : 4 columns :                  
     1. seqid                         
     2. position                      
     3. EHH                           
     4. ref or alt [0 == ref]         
Params:
       required: t,target   <STRING>  A zero base comma separated list of target
                                     individuals corresponding  to VCF columns  
       required: r,region   <STRING>  A tabix compliant genomic range           
                                     format: "seqid:start-end" or "seqid"  
       required: f,file     <STRING>  Proper formatted and phased VCF.          
       required: y,type     <STRING>  Genotype likelihood format: GT,PL,GL,GP   
       required: p,position <INT>     Variant position to melt.                 
       optional: a,af       <DOUBLE>  Alternative  alleles with frequencies less   
                                     than [0.05] are skipped.                  

Type: statistics



```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[meltEHH.cpp](https://github.com/vcflib/vcflib/blob/master/src/meltEHH.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
