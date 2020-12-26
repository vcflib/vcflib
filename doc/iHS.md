% IHS(1) iHS (vcflib) | iHS (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

iHS

# SYNOPSIS

Usage: iHS --target 0,1,2,3,4,5,6,7 --file my.phased.vcf \ --region chr1:1-1000 > STDOUT 2> STDERR

# DESCRIPTION

iHS calculates the integrated ratio of haplotype decay between the reference and non-reference allele. Output : 4 columns : 1. seqid 2. position 3. target allele frequency 4. integrated EHH (alternative) 5. integrated EHH (reference) 6. iHS ln(iEHHalt/iEHHref)



# OPTIONS

```


     7. != 0 integration failure                    

     8. != 0 integration failure                    

Params:
       required: t,target  <STRING>  A zero base comma separated list of target
                                     individuals corresponding  to VCF columns  
       required: r,region  <STRING>  A tabix compliant genomic range           
                                     format: "seqid:start-end" or "seqid"  
       required: f,file    <STRING>  Proper formatted and phased VCF.          
       required: y,type    <STRING>  Genotype likelihood format: GT,PL,GL,GP   
       optional: a,af      <DOUBLE>  Alternative  alleles with frequencies less   
                                     than [0.05] are skipped.                  
       optional: x,threads <INT>     Number of CPUS [1].                       
       recommended: g,gen <STRING>   A PLINK formatted map file.               


```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[iHS.cpp](https://github.com/vcflib/vcflib/blob/master/src/iHS.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
