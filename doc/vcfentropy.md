% VCFENTROPY(1) vcfentropy (vcflib) | vcfentropy (VCF metrics)
% Erik Garrison and vcflib contributors

# NAME

**vcfentropy**

# SYNOPSIS

**vcfentropy** [options] <vcf file>

# DESCRIPTION

Annotate VCF records with the Shannon entropy of flanking sequence. Anotates the output VCF file with, for each record, EntropyLeft, EntropyRight, EntropyCenter, which are the entropies of the sequence of the given window size to the left, right, and center of the record. Also adds EntropyRef and EntropyAlt for each alt.



# OPTIONS

```

options:
    -f, --fasta-reference  FASTA reference file to use to obtain flanking sequences
    -w, --window-size      Size of the window over which to calculate entropy



Type: metrics

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

[vcfentropy.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfentropy.cpp)

# LICENSE

Copyright 2011-2023 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
