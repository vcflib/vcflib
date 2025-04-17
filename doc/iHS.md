% IHS(1) iHS (vcflib) | iHS (VCF statistics)
% Erik Garrison and vcflib contributors

# NAME

**iHS**

# SYNOPSIS

**iHS** --target 0,1,2,3,4,5,6,7 --file my.phased.vcf \ --region chr1:1-1000 > STDOUT 2> STDERR

# DESCRIPTION

**iHS** calculates the integrated haplotype score which measures the relative decay of extended haplotype homozygosity (EHH) for the reference and alternative alleles at a site (see: voight et al. 2006, Spiech & Hernandez 2014).



# OPTIONS

```


Our code is highly concordant with both implementations mentioned. However, we do not set an upper limit to the allele frequency.  **iHS** can be run without a genetic map, in which case the change in EHH is integrated over a constant.  Human genetic maps for GRCh36 and GRCh37 (hg18 & hg19) can be found at: http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ . **iHS** by default interpolates SNV positions to genetic position (you don't need a genetic position for every VCF entry in the map file).

**iHS** analyses requires normalization by allele frequency.  It is important that **iHS** is calculated over large regions so that the normalization does not down weight real signals.  For genome-wide runs it is recommended to run slightly overlapping windows and throwing out values that fail integration (columns 7 & 8 in the output) and then removing duplicates by using the 'sort' and 'uniq' linux commands.  Normalization of the output is as simple as running 'normalize-**iHS**'.



     **iHS** calculates the integrated ratio of haplotype decay between the reference and non-reference allele.
Output : 4 columns :
     1. seqid
     2. position
     3. target allele frequency
     4. integrated EHH (alternative)
     5. integrated EHH (reference)
     6. **iHS** ln(iEHHalt/iEHHref)
     7. != 0 integration failure
     8. != 0 integration failure

Params:
       required: t,target  <STRING>  A zero base comma separated list of target
                                     individuals corresponding to VCF columns
       required: r,region  <STRING>  A tabix compliant genomic range
                                     format: "seqid:start-end" or "seqid"
       required: f,file    <STRING>  Proper formatted and phased VCF.
       required: y,type    <STRING>  Genotype likelihood format: GT,PL,GL,GP
       optional: a,af      <DOUBLE>  Alternative alleles with frquences less
                                     than [0.05] are skipped.
       optional: x,threads <INT>     Number of CPUS [1].
       recommended: g,gen <STRING>   A PLINK formatted map file.



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

[iHS.cpp](https://github.com/vcflib/vcflib/blob/master/src/iHS.cpp)

# LICENSE

Copyright 2011-2025 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2025 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
