% VCFALLELICPRIMITIVES(1) vcfallelicprimitives (vcflib) | vcfallelicprimitives (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcfallelicprimitives**

# SYNOPSIS

**vcfallelicprimitives** [options] [file]

# DESCRIPTION

Realign reference and alternate alleles with WFA, parsing out the primitive alleles
into multiple VCF records. New records have IDs that reference the source record ID.
Genotypes are handled. Deletion alleles will result in haploid (missing allele) genotypes.

# OPTIONS

```

options:
    -m, --use-mnps          Retain MNPs as separate events (default: false).
    -t, --tag-parsed FLAG   Annotate decomposed records with the source record position
                            (default: ORIGIN).
    -L, --max-length LEN    Do not manipulate records in which either the ALT or
                            REF is longer than LEN (default: unlimited).
    -k, --keep-info         Maintain site and allele-level annotations when decomposing.
                            Note that in many cases, such as multisample VCFs, these won't
                            be valid post-decomposition.  For biallelic loci in single-sample
                            VCFs, they should be usable with caution.
    -d, --debug             debug mode.

Type: transformation

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

[vcfallelicprimitives.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfallelicprimitives.cpp)

# LICENSE

Copyright 2011-2022 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
