% VCF2FASTA(1) vcf2fasta (vcflib) | vcf2fasta (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcf2fasta

# SYNOPSIS

usage: ./build/vcf2fasta [options] [file]

# DESCRIPTION

options: -f, --reference REF Use this reference when decomposing samples. -p, --prefix PREFIX Affix this output prefix to each file, none by default -P, --default-ploidy N Set a default ploidy for samples which do not have information in the first record (2).

# OPTIONS

```


Outputs sample_seq:N.fa for each sample, reference sequence, and chromosomal copy N in [0,1... ploidy].
Each sequence in the fasta file is named using the same pattern used for the file name, allowing them to be combined.

```

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcf2fasta.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcf2fasta.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
