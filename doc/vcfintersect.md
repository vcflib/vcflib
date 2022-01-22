% VCFINTERSECT(1) vcfintersect (vcflib) | vcfintersect (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfintersect**

# SYNOPSIS

**vcfintersect** [options] [<vcf file>]

# DESCRIPTION

VCF set analysis



# OPTIONS

```


options:
    -b, --bed FILE            use intervals provided by this BED file
    -R, --region REGION       use 1-based tabix-style region (e.g. chrZ:10-20), multiples allowed
    -S, --start-only          don't use the reference length information in the record to determine
                              overlap status, just use the start posiion
    -v, --invert              invert the selection, printing only records which would
                                not have been printed out
    -i, --intersect-vcf FILE  use this VCF for set intersection generation
    -u, --union-vcf FILE      use this VCF for set union generation
    -w, --window-size N       compare records up to this many bp away (default 30)
    -r, --reference FILE      FASTA reference file, required with -i and -u
    -l, --loci                output whole loci when one alternate allele matches
    -m, --ref-match           intersect on the basis of record REF string
    -t, --tag TAG             attach TAG to each record's info field if it would intersect
    -V, --tag-value VAL       use this value to indicate that the allele is passing
                              '.' will be used otherwise.  default: 'PASS'
    -M, --merge-from FROM-TAG
    -T, --merge-to   TO-TAG   merge from FROM-TAG used in the -i file, setting TO-TAG
                              in the current file.

For bed-vcf intersection, alleles which fall into the targets are retained.

Haplotype aware intersection, union and complement. Use for intersection and union of VCF files: unify on equivalent alleles within window-size bp
as determined by haplotype comparison alleles.

type: transformation

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

[vcfintersect.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfintersect.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
