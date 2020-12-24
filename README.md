# vcflib

### A C++ library for parsing and manipulating VCF files.

![Github-CI](https://github.com/vcflib/vcflib/workflows/CI/badge.svg) [![Travis-CI](https://travis-ci.org/vcflib/vcflib.svg?branch=master)](https://travis-ci.org/vcflib/vcflib) [![AnacondaBadge](https://anaconda.org/bioconda/vcflib/badges/installer/conda.svg)](https://anaconda.org/bioconda/vcflib) [![DL](https://anaconda.org/bioconda/vcflib/badges/downloads.svg)](https://anaconda.org/bioconda/vcflib) [![BrewBadge](https://img.shields.io/badge/%F0%9F%8D%BAbrew-vcflib-brightgreen.svg)](https://github.com/brewsci/homebrew-bio) [![GuixBadge](https://img.shields.io/badge/gnuguix-vcflib-brightgreen.svg)](https://www.gnu.org/software/guix/packages/V/) [![DebianBadge](https://badges.debian.net/badges/debian/testing/libvcflib-dev/version.svg)](https://packages.debian.org/testing/libvcflib-dev)
[![C++0x](https://img.shields.io/badge/Language-C++0x-steelblue.svg)](https://www.cprogramming.com/c++11/what-is-c++0x.html)
[![Gitter](https://badges.gitter.im/ekg/vcflib.svg)](https://gitter.im/ekg/vcflib?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## overview

The [Variant Call Format (VCF)](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)
is a flat-file, tab-delimited textual format
intended to concisely describe reference-indexed variations between individuals.
VCF provides a common interchange format for the description of variation in individuals and populations of samples,
and has become the _defacto_ standard reporting format for a wide array of genomic variant detectors.

vcflib provides methods to manipulate and interpret sequence variation as it can be described by VCF.
It is both:

 * an API for parsing and operating on records of genomic variation as it can be described by the VCF format,
 * and a collection of command-line utilities for executing complex manipulations on VCF files.

The API itself provides a quick and extremely permissive method to read and write VCF files.
Extensions and applications of the library provided in the included utilities (*.cpp) comprise the vast bulk of the library's utility for most users.

Short index:

- [Install](#INSTALL)
- [Usage](#USAGE)
- [TOOLS](#TOOLS)
  * [Transformation](#tools-for-transformation)
  * [Metrics](#tools-for-metrics)
  * [Statistics](#tools-for-statistics)
  * [Validation](#tools-for-validation)
- [Link library](#link-library)
- [Build from source](#build-from-source)
- [Development](#Development)
- [LICENSE](#LICENSE)

## INSTALL

### [Bioconda](https://bioconda.github.io/user/install.html)

Conda installs in user land without root access

```sh
conda install -c bioconda vcflib
```

### [Homebrew](https://brew.sh)

Homebrew installs on Linux and Mac OSX

```sh
brew install brewsci/bio/vcflib
```

### [Debian](https://debian.org/)

For Debian and Ubuntu

```sh
apt-get install libvcflib-tools libvcflib-dev
```

### [GNU Guix](https://guix.gnu.org/)

We develop against guix

```sh
guix package -i vcflib
```

## USAGE

vcflib provides a variety of tools and functions for VCF manipulation.

Quick overview:

### transform

 * Generate **haplotype-aware intersections** ([vcfintersect](#vcfintersect) -i), **unions** (vcfintersect -u), and **complements** (vcfintersect -v -i).
 * **Overlay-merge** multiple VCF files together, using provided order as precedence ([vcfoverlay](#vcfoverlay)).
 * **Combine** multiple VCF files together, handling samples when alternate allele descriptions are identical ([vcfcombine](#vcfcombine)).
 * **Validate** the integrity and identity of the VCF by verifying that the VCF record's REF matches a given reference file ([vcfcheck](#vcfcheck)).

### format conversion

 * Convert a VCF file into a per-allele or per-genotype **tab-separated (.tsv)** file ([vcf2tsv](#vcf2tsv)).
 * Store a VCF file in an **SQLite3** database (vcf2sqlite.py).
 * Make a **BED file** from the intervals in a VCF file (vcf2bed.py).

### filtering and subsetting

 * **Filter** variants and genotypes using arbitrary expressions based on values in the INFO and sample fields ([vcffilter](#vcffilter)).
 * **Randomly sample** a subset of records from a VCF file, given a rate ([vcfrandomsample](#vcfrandomsample)).
 * **Select variants** of a certain type (vcfsnps, vcfbiallelic, vcfindels, vcfcomplex, etc.)

### annotation

 * **Annotate** one VCF file with fields from the INFO column of another, based on position ([vcfaddinfo](#vcfaddinfo), [vcfintersect](#vcfintersect)).
 * Incorporate annotations or targets provided by a *BED* file ([vcfannotate](#vcfannotate), [vcfintersect](#vcfintersect)).
 * Examine **genotype correspondence** between two VCF files by annotating samples in one file with genotypes from another ([vcfannotategenotypes](#vcfannotategenotypes)).
 * Annotate variants with the **distance** to the nearest variant ([vcfdistance](#vcfdistance)).
 * Count the number of alternate alleles represented in samples at each variant record ([vcfaltcount](#vcfaltcount)).
 * **Subset INFO fields** to decrease file size and processing time ([vcfkeepinfo](#vcfkeepinfo)).
 * Lighten up VCF files by keeping only a **subset of per-sample information** ([vcfkeepgeno](#vcfkeepgeno)).
 * **Numerically index** alleles in a VCF file ([vcfindex](#vcfindex)).

### samples

 * Quickly obtain the **list of samples** in a given VCF file ([vcfsamplenames](#vcfsamplenames)).
 * **Remove samples** from a VCF file ([vcfkeepsamples](#vcfkeepsamples), [vcfremovesamples](#vcfremovesamples)).

### ordering

 * **Sort variants** by genome coordinate ([vcfstreamsort](#vcfstreamsort)).
 * **Remove duplicate** variants in vcfstreamsort'ed files according to their REF and ALT fields ([vcfuniq](#vcfuniq)).

### variant representation

 * **Break multiallelic** records into multiple records ([vcfbreakmulti](#vcfbreakmulti)), retaining allele-specific INFO fields.
 * **Combine overlapping biallelic** records into a single record ([vcfcreatemulti](#vcfcreatemulti)).
 * **Decompose complex variants** into a canonical SNP and indel representation ([vcfallelicprimitives](#vcfallelicprimitives)), generating phased genotypes for available samples.
 * **Reconstitute complex variants** provided a phased VCF with samples ([vcfgeno2haplo](#vcfgeno2haplo)).
 * **Left-align indel and complex variants** ([vcfleftalign](#vcfleftalign)).

### genotype manipulation

 * **Set genotypes** in a VCF file provided genotype likelihoods in the GL field ([vcfglxgt](#vcfglxgt)).
 * Establish putative **somatic variants** using reported differences between germline and somatic samples ([vcfsamplediff](#vcfsamplediff)).
 * Remove samples for which the reported genotype (GT) and observation counts disagree (AO, RO) ([vcfremoveaberrantgenotypes](#vcfremoveaberrantgenotypes)).

### interpretation and classification of variants

 * Obtain aggregate **statistics** about VCF files ([vcfstats](#vcfstats)).
 * Print the **receiver-operating characteristic (ROC)** of one VCF given a truth set ([vcfroc](#vcfroc)).
 * Annotate VCF records with the **Shannon entropy** of flanking sequence ([vcfentropy](#vcfentropy)).
 * Calculate the heterozygosity rate ([vcfhetcount](#vcfhetcount)).
 * Generate potential **primers** from VCF records ([vcfprimers](#vcfprimers)), to check for genome uniqueness.
 * Convert the numerical represenation of genotypes provided by the GT field to a **human-readable genotype format** ([vcfgenotypes](#vcfgenotypes)).
 * Observe how different alignment parameters, including context and entropy-dependent ones, influence **variant classification and interpretation** ([vcfremap](#vcfremap)).
 * **Classify variants** by annotations in the INFO field using a self-organizing map ([vcfsom](#vcfsom)); **re-estimate their quality** given known variants.

A number of "helper" perl and python3 scripts (e.g. vcf2bed.py,
vcfbiallelic) further extend functionality.

Users are encouraged to drive the utilities in the library in a
streaming fashion, using pipes, to fully utilize resources on
multi-core systems during interactive work.  Piping provides a
convenient method to interface with other libraries (vcf-tools,
BedTools, GATK, htslib, [bio-vcf](https://github.com/vcflib/bio-vcf),
bcftools, [freebayes](https://github.com/freebayes)) which interface
via VCF files, allowing the composition of an immense variety of
processing functions.



# TOOLS

<!--

    >>> from pytest.rtest import run_stdout, head, cat

-->

## tools for transformation

### vcfintersect

```sh

>>> cat("vcfintersect -h")
vcflib set analysis
>
usage: vcfintersect [options] [<vcf file>]
>
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
>
For bed-vcf intersection, alleles which fall into the targets are retained.
>
For vcf-vcf intersection and union, unify on equivalent alleles within window-size bp
as determined by haplotype comparison alleles.

```

Docs [vcfintersect.md](./test/pytest/vcfintersect.md) | source code [vcfintersect.cpp](./src/vcfintersect.cpp)

### vcf2tsv


```sh

>>> cat("vcf2tsv -h")
usage: vcf2tsv [-n null_string] [-g] [vcf file]
>
Converts stdin or given VCF file to tab-delimited format, using null string to replace empty values in the table.
Specifying -g will output one line per sample with genotype information.
When there is more than one alt allele there will be multiple rows, one for each allele and, the info will match the 'A' index

```

vcf2tsv converts a CSV to a tabulated test file, e.g.

       vcf2tsv ../samples/sample.vcf


```python

>>> head("vcf2tsv ../samples/sample.vcf")
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  AA      AC      AF      AN      DB      DP      H2      NS
19      111     .       A       C       9.6     .       .       .       .       .       .       .       .       .
19      112     .       A       G       10      .       .       .       .       .       .       .       .       .
20      14370   rs6054257       G       A       29      PASS    .       .       0.5     .       .       14      .       3

```

Docs [vcf2tsv.md](./test/pytest/vcf2tsv.md) | source code [vcf2tsv.cpp](./src/vcf2tsv.cpp)

### vcfaddinfo

```sh

>>> cat("vcfaddinfo -h")
usage: vcfaddinfo <vcf file> <vcf file>
Adds info fields from the second file which are not present in the first vcf file.

```

Source code [vcfaddinfo.cpp](./src/vcfaddinfo.cpp)

### vcfannotate

```sh

>>> cat("vcfannotate -h")
usage: vcfannotate [options] [<vcf file>]
>
options:
    -b, --bed   use annotations provided by this BED file
    -k, --key   use this INFO field key for the annotations
    -d, --default  use this INFO field key for records without annotations
>
Intersect the records in the VCF file with targets provided in a BED file.
Intersections are done on the reference sequences in the VCF file.
If no VCF filename is specified on the command line (last argument) the VCF
read from stdin.

```

Source code [vcfannotate.cpp](./src/vcfannotate.cpp)


### vcfallelicprimitives


```sh

>>> cat("vcfallelicprimitives -h")
usage: vcfallelicprimitives [options] [file]
>
options:
    -m, --use-mnps          Retain MNPs as separate events (default: false).
    -t, --tag-parsed FLAG   Tag records which are split apart of a complex allele with this flag.
    -L, --max-length LEN    Do not manipulate records in which either the ALT or
                            REF is longer than LEN (default: 200).
    -k, --keep-info         Maintain site and allele-level annotations when decomposing.
                            Note that in many cases, such as multisample VCFs, these won't
                            be valid post-decomposition.  For biallelic loci in single-sample
                            VCFs, they should be usable with caution.
    -g, --keep-geno         Maintain genotype-level annotations when decomposing.  Similar
                            caution should be used for this as for --keep-info.
>
If multiple allelic primitives (gaps or mismatches) are specified in
a single VCF record, split the record into multiple lines, but drop all
INFO fields.  Does not handle genotypes (yet).  MNPs are split into
multiple SNPs unless the -m flag is provided.  Records generated by splits have th

```

Source code [vcfallelicprimitives.cpp](./src/vcfallelicprimitives.cpp)

### vcfannotategenotypes

```sh

>>> cat("vcfannotategenotypes -h")
usage: vcfannotategenotypes <annotation-tag> <vcf file> <vcf file>
annotates genotypes in the first file with genotypes in the second
adding the genotype as another flag to each sample filed in the first file.
annotation-tag is the name of the sample flag which is added to store the annotation.
also adds a 'has_variant' flag for sites where the second file has a variant.

```

Source code [vcfannotategenotypes.cpp](./src/vcfannotategenotypes.cpp)

### vcfbreakmulti


```sh

>>> cat("vcfbreakmulti -h")
usage: vcfbreakmulti [options] [file]
>
If multiple alleles are specified in a single record, break the record into
multiple lines, preserving allele-specific INFO fields.

```

Source code [vcfbreakmulti.cpp](./src/vcfbreakmulti.cpp)



### vcfcleancomplex

Removes reference-matching sequence from complex alleles and adjusts
records to reflect positional change.

```sh

>>> cat("vcfcleancomplex")
usage: vcfcleancomplex <vcf file>
outputs a VCF stream in which 'long' non-complexalleles have their position corrected.
assumes that VCF records can't overlap 5'->3'

```

Source code [vcfbreakmulti.cpp](./src/vcfbreakmulti.cpp)

### vcfcombine

```sh

>>> cat("vcfcombine -h")
usage: vcfcombine [vcf file] [vcf file] ...
>
Combines VCF files positionally, combining samples when sites and alleles are identical.
Any number of VCF files may be combined.  The INFO field and other columns are taken from
one of the files which are combined when records in multiple files match.  Alleles must
have identical ordering to be combined into one record.  If they do not, multiple records
will be emitted.
>
options:
    -h --help           This text.
    -v --version        Print version.
    -r --region REGION  A region specifier of the form chrN:x-y to bound the merge

```

Source code [vcfcombine.cpp](./src/vcfcombine.cpp)

### vcfdistance

Adds a value to each VCF record indicating the distance to the nearest
variant in the file.

```sh

>>> cat("vcfdistance -h")
usage: vcfdistance [customtagname] < [vcf file]
adds a tag to each variant record which indicates the distance
to the nearest variant.
(defaults to BasesToClosestVariant if no custom tag name is given.

```

Source code [vcfdistance.cpp](./src/vcfdistance.cpp)

## tools for metrics

### vcfaltcount

Counts the number of alternate alleles in the record.

```sh

>>> cat("vcfaltcount")
usage: vcfaltcount <vcf file>
count the number of alternate alleles in all records in the vcf file

```

Source code [vcfaltcount.cpp](./src/vcfaltcount.cpp)

### vcfcommonsamples

```sh

>>> cat("vcfcommonsamples -h")
usage: vcfcommonsamples <vcf file> <vcf file>
outputs each record in the first file, removing samples not present in the second

```

Source code [vcfcommonsamples.cpp](./src/vcfcommonsamples.cpp)



### vcfcountalleles

Counts the total number of alleles in the input.

```sh

>> cat("vcfcountalleles -h")

```

Source code [vcfcountalleles.cpp](./src/vcfcountalleles.cpp)



### vcfcreatemulti


```sh

>>> cat("vcfcreatemulti -h")
usage: vcfcreatemulti [options] [file]
>
If overlapping alleles are represented across multiple records, merge
them into a single record.  Currently only for indels.

```

Source code [vcfcreatemulti.cpp](./src/vcfcreatemulti.cpp)



## tools for statistics

### vcfafpath

Uses allele frequencies in the AF info column to estimate phylogeny at multiallelic sites.


## tools for validation

### vcfcheck

```sh

>>> cat("vcfcheck -h")
usage: vcfcheck [options] <vcf file>
>
options:
    -f, --fasta-reference  FASTA reference file to use to obtain primer sequences
    -x, --exclude-failures If a record fails, don't print it.  Otherwise do.
    -k, --keep-failures    Print if the record fails, otherwise not.
    -h, --help       Print this message.
    -v, --version    Print version.
>
Verifies that the VCF REF field matches the reference as described.
>

```

Source code [vcfcheck.cpp](./src/vcfcheck.cpp)

<!--
  I am here organising @pjotrp


```sh

>> cat("vcfxx -h")

```

Source code [vcfxx.cpp](./src/vcfxx.cpp)


-->



### vcfentropy

```sh

>>> cat("vcfentropy -h")
usage: vcfentropy [options] <vcf file>
>
options:
    -f, --fasta-reference  FASTA reference file to use to obtain flanking sequences
    -w, --window-size      Size of the window over which to calculate entropy
>
Anotates the output VCF file with, for each record, EntropyLeft, EntropyRight,
EntropyCenter, which are the entropies of the sequence of the given window size to the
left, right, and center  of the record.  Also adds EntropyRef and EntropyAlt for each alt.
>

```


### vcffilter

```sh

>>> cat("vcffilter -h")
vcflib filtering
>
usage: vcffilter [options] <vcf file>
>
options:
    -f, --info-filter     specifies a filter to apply to the info fields of records,
                          removes alleles which do not pass the filter
    -g, --genotype-filter specifies a filter to apply to the genotype fields of records
    -k, --keep-info       used in conjunction with '-g', keeps variant info, but removes genotype
    -s, --filter-sites    filter entire records, not just alleles
    -t, --tag-pass        tag vcf records as positively filtered with this tag, print all records
    -F, --tag-fail        tag vcf records as negatively filtered with this tag, print all records
    -A, --append-filter   append the existing filter tag, don't just replace it
    -a, --allele-tag      apply -t on a per-allele basis.  adds or sets the corresponding INFO field tag
    -v, --invert          inverts the filter, e.g. grep -v
    -o, --or              use logical OR instead of AND to combine filters
    -r, --region          specify a region on which to target the filtering, requires a BGZF
                          compressed file which has been indexed with tabix.  any number of
                          regions may be specified.
>
Filter the specified vcf file using the set of filters.
Filters are specified in the form "<ID> <operator> <value>:
 -f "DP > 10"  # for info fields
 -g "GT = 1|1" # for genotype fields
 -f "CpG"  # for 'flag' fields
>
Operators can be any of: =, !, <, >, |, &
>
Any number of filters may be specified.  They are combined via logical AND
unless --or is specified on the command line.  Obtain logical negation through
the use of parentheses, e.g. "! ( DP = 10 )"
>
For convenience, you can specify "QUAL" to refer to the quality of the site, even
though it does not appear in the INFO fields.
>

```

Docs [vcffilter.md](./test/pytest/vcffilter.md) | source code [vcffilter.cpp](./src/vcffilter.cpp)


### vcffixup

```sh

>>> cat("vcffixup -h")
usage: vcffixup <vcf file>
outputs a VCF stream where AC and NS have been generated for each record using sample genotypes

```

Count the allele frequencies across alleles present in each record in the VCF file. (Similar to vcftools --freq.)

Uses genotypes from the VCF file to correct AC (alternate allele count), AF
(alternate allele frequency), NS (number of called), in the VCF records.  For
example:

    % vcfkeepsamples file.vcf NA12878 | vcffixup - | vcffilter -f "AC > 0"

Would downsample file.vcf to only NA12878, removing sites for which the sample
was not called as polymorphic.


### vcfflatten

```sh

>>> cat("vcfflatten -h")
usage: vcfflatten [file]
>
Removes multi-allelic sites by picking the most common alternate.  Requires
allele frequency specification 'AF' and use of 'G' and 'A' to specify the
fields which vary according to the Allele or Genotype. VCF file may be
specified on the command line or piped as stdin.

```

### vcfgeno2haplo


```sh

>>> cat("vcfgeno2haplo -h")
usage: vcfgeno2haplo [options] [<vcf file>]
>
options:
    -h, --help              Print this message
    -v, --version           Print version
    -r, --reference FILE    FASTA reference file
    -w, --window-size N     Merge variants at most this many bp apart (default 30)
    -o, --only-variants     Don't output the entire haplotype, just concatenate
                            REF/ALT strings (delimited by ":")
>
Convert genotype-based phased alleles within --window-size into haplotype alleles.
Will break haplotype construction when encountering non-phased genotypes on input.
>

```

### vcfgenotypecompare

```sh

>>> cat("vcfgenotypecompare -h")
usage: vcfgenotypecompare <other-genotype-tag> <vcf file>
adds statistics to the INFO field of the vcf file describing the
amount of discrepancy between the genotypes (GT) in the vcf file and the
genotypes reported in the <other-genotype-tag>.  use this after
vcfannotategenotypes to get correspondence statistics for two vcfs.

```

### vcfgenotypes

Converts numerical representation of genotypes (standard in GT field) to the
alleles provided in the call's ALT/REF fields.


### vcfglxgt

    usage: vcfglxgt [options] <vcf file>

    options:
        -n, --fix-null-genotypes   only apply to null and partly-null genotypes

Set genotypes using the maximum genotype likelihood for each sample.



### vcfhetcount

Count the number of heterozygotes in the input VCF.


### vcfhethomratio

Provides the ratio between heterozygotes and homozygotes.

### vcfindex

Adds a field (id) which contains an allele-specific numerical index.


### vcfkeepgeno

    usage: vcfkeepgeno <vcf file> [FIELD1] [FIELD2] ...

Outputs each record in the vcf file, removing FORMAT fields not listed on the command line from sample specifications in the output.


### vcfkeepinfo

    usage: vcfkeepinfo <vcf file> [FIELD1] [FIELD2] ...

 Outputs each record in the vcf file, removing INFO fields not listed on the command line.


### vcfkeepsamples

    usage: vcfkeepsamples <vcf file> [SAMPLE1] [SAMPLE2] ...

Outputs each record in the vcf file, removing samples not listed on the command line.


### vcfleftalign

Left-align indels and complex variants in the input using a pairwise ref/alt
alignment followed by a heuristic, iterative left realignment process that
shifts indel representations to their absolute leftmost (5') extent.  This is
the same procedure used in the internal left alignment in freebayes, and can be
used when preparing VCF files for input to freebayes to decrease positional
representation differences between the input alleles and left-realigned
alignments.

    usage: vcfleftalign [options] [file]

    options:
        -r, --reference FILE  Use this reference as a basis for realignment.
        -w, --window N        Use a window of this many bp when left aligning (150).

Left-aligns variants in the specified input file or stdin.
Window size is determined dynamically according to the entropy of the regions flanking the indel.
These must have entropy > 1 bit/bp, or be shorter than ~5kb.


### vcflength

Adds the length of the variant record (in [-/+]) relative to the reference allele to each VCF record.


### vcfnumalt

Annotates the VCF stream on stdin with the number of alternate alleles at the site.


### vcfoverlay

VCF merging.

```sh

>>> cat("vcfoverlay -h")
usage: vcfoverlay [options] [<vcf file> ...]
>
options:
    -h, --help       this dialog
    -v, --version    prints version
>
Overlays records in the input vcf files in the order in which they appear.

```

### vcfparsealts

Demonstration of alternate allele parsing method.  This method uses pairwise
alignment of REF and ALTs to determine component allelic primitives for each
alternate allele.

Use `vcfallelicprimitives` to decompose records while preserving format.


### vcfprimers

    usage: vcfprimers [options] <vcf file>

    options:
        -f, --fasta-reference  FASTA reference file to use to obtain primer sequences
        -l, --primer-length    The length of the primer sequences on each side of the variant

For each VCF record, extract the flanking sequences, and write them to stdout as FASTA records suitable for alignment.
This tool is intended for use in designing validation experiments.
Primers extracted which would flank all of the alleles at multi-allelic sites.
The name of the FASTA "reads" indicates the VCF record which they apply to.
The form is >CHROM_POS_LEFT for the 3' primer and >CHROM_POS_RIGHT for the 5' primer, for example:

    >20_233255_LEFT
    CCATTGTATATATAGACCATAATTTCTTTATCCAATCATCTGTTGATGGA
    >20_233255_RIGHT
    ACTCAGTTGATTCCATACCTTTGCCATCATGAATCATGTTGTAATAAACA



### vcfrandomsample

    usage: vcfrandomsample [options] [<vcf file>]

    options:
        -r, --rate RATE      base sampling probability per locus
        -s, --scale-by KEY   scale sampling likelihood by this Float info field
        -p, --random-seed N  use this random seed

Randomly sample sites from an input VCF file, which may be provided as stdin.
Scale the sampling probability by the field specified in KEY.
This may be used to provide uniform sampling across allele frequencies, for instance.


### vcfremap

    usage: vcfremap [options] [<vcf file>]

    options:
        -w, --ref-window-size N      align using this many bases flanking each side of the reference allele
        -s, --alt-window-size N      align using this many flanking bases from the reference around each alternate allele
        -r, --reference FILE         FASTA reference file, required with -i and -u
        -m, --match-score N          match score for SW algorithm
        -x, --mismatch-score N       mismatch score for SW algorithm
        -o, --gap-open-penalty N     gap open penalty for SW algorithm
        -e, --gap-extend-penalty N   gap extension penalty for SW algorithm
        -z, --entropy-gap-open       use entropy scaling for the gap open penalty
        -R, --repeat-gap-extend N    penalize non-repeat-unit gaps in repeat sequence
        -a, --adjust-vcf TAG         supply a new cigar as TAG in the output VCF

For each alternate allele, attempt to realign against the reference with lowered gap open penalty.
If realignment is possible, adjust the cigar and reference/alternate alleles.


### vcfremoveaberrantgenotypes

Strips genotypes which are homozygous but have observations implying
heterozygosity.  Requires RA (reference allele observation) and AA (alternate
allele observation) for each genotype.


### vcfremovesamples

    usage: vcfremovesamples <vcf file> [SAMPLE1] [SAMPLE2] ...

Outputs each record in the vcf file, removing samples listed on the command line.


### vcfroc

    usage: vcfroc [options] [<vcf file>]

    options:
        -t, --truth-vcf FILE      use this VCF as ground truth for ROC generation
        -w, --window-size N       compare records up to this many bp away (default 30)
        -r, --reference FILE      FASTA reference file

Generates a pseudo-ROC curve using sensitivity and specificity estimated against a putative truth set.
Thresholding is provided by successive QUAL cutoffs.


### vcfsamplediff

    usage: vcfsamplediff <tag> <sample> <sample> [ <sample> ... ] <vcf file>

Tags each record where the listed sample genotypes differ with <tag>
The first sample is assumed to be germline, the second somatic.
Each record is tagged with <tag>={germline,somatic,loh} to specify the type of variant given the genotype difference between the two samples.


### vcfsamplenames

Prints the names of the samples in the VCF file.


### vcfsom

    usage: vcfsom [options] [vcf file]

    training:
        vcfsom -s output.som -f "AF DP ABP" training.vcf

    application:
        vcfsom -a output.som -f "AF DP ABP" test.vcf >results.vcf

vcfsom trains and/or applies a self-organizing map to the input VCF data on stdin, adding two columns for the x and y coordinates of the winning neuron in the network and an optional euclidean distance from a given node (--center).

If a map is provided via --apply,  map will be applied to input without training.
Automated filtering to an estimated FP rate is

    options:

        -h, --help             this dialog

    training:

        -f, --fields "FIELD ..."  INFO fields to provide to the SOM
        -a, --apply FILE       apply the saved map to input data to FILE
        -s, --save  FILE       train on input data and save the map to FILE
        -t, --print-training-results
                               print results of SOM on training input
                               (you can also just use --apply on the same input)
        -x, --width X          width in columns of the output array
        -y, --height Y         height in columns of the output array
        -i, --iterations N     number of training iterations or epochs
        -d, --debug            print timing information

    recalibration:

        -c, --center X,Y       annotate with euclidean distance from center
        -p, --paint-true VCF   use VCF file to annotate true variants (multiple)
        -f, --paint-false VCF  use VCF file to annotate false variants (multiple)
        -R, --paint-tag TAG    provide estimated FDR% in TAG in variant INFO
        -N, --false-negative   replace FDR% (false detection) with FNR% (false negative)


### vcfstats

    usage: vcfstats [options] <vcf file>

        -r, --region          specify a region on which to target the stats, requires a BGZF
                              compressed file which has been indexed with tabix.  any number of
                              regions may be specified.
        -a, --add-info        add the statistics intermediate information to the VCF file,
                              writing out VCF records instead of summary statistics
        -l, --no-length-frequency    don't out the indel and mnp length-frequency spectra
        -m, --match-score N          match score for SW algorithm
        -x, --mismatch-score N       mismatch score for SW algorithm
        -o, --gap-open-penalty N     gap open penalty for SW algorithm
        -e, --gap-extend-penalty N   gap extension penalty for SW algorithm

Prints statistics about variants in the input VCF file.


### vcfstreamsort

Reads VCF on stdin and guarantees that the positional order is correct provided out-of-order
variants are no more than 100 positions in the VCF file apart.


### vcfuniq

Like GNU uniq, but for VCF records.  Remove records which have the same positon, ref, and alt
as the previous record.


### vcfuniqalleles

For each record, remove any duplicate alternate alleles that may have resulted from merging
separate VCF files.

## GPAT++

The application of population genomics to non-model organisms is greatly facilitated by the low cost of next generation sequencing (NGS). Barriers, however, exist for using NGS data for population level analyses. Traditional population genetic metrics, such as Fst, are not robust to the genotyping errors inherent in noisy NGS data. Additionally, many older software tools were never designed to handle the volume of data produced by NGS pipelines. To overcome these limitations we have developed a flexible software library designed specifically for large and noisy NGS datasets. The Genotype Phenotype Association Toolkit (GPAT++) implements both traditional and novel population genetic methods in a single user-friendly framework. GPAT consists of a suite of command-line tools and a Perl API that programmers can use to develop new applications. To date GPAT++ has been used successfully to identity genotype-phenotype associations in several real-world datasets including: domestic pigeons, Pox virus and pine rust fungus. GPAT++ is open source and freely available for academic use.

### Functions
 - [X] Basic population stats (Af, Pi, eHH, oHet, genotypeCounts)
 - [X] Several flavors of Fst
 - [X] Linkage
 - [X] Association testing (genotypic and pooled data)
 - [X] Haplotype methods (hapLrt)
 - [X] Smoothing
 - [X] Permutation
 - [X] Plotting

### Documentation , basic usage, FAQ


1. Most GPAT++ tools write to both STDERR and STDOUT.
2. All GPAT++ tools group individuals using a zero-based comma separated index (e.g. 0,1,2 ; first three individuals in VCF)
3. Some GPAT++ tools (haplotype methods) require a region.
4. What is the genotype likelihood format?  When in doubt use GT! Only a few GPAT++ tools make use of the genotype likelihoods.
 GT: The genotype is correct
 GL: Genotype likelihood (Freebayes)
 GP: Genotype probability (Beagle)
 PL: Scaled genotype likelihood (GATK)
5. pFst is the only tool that will work on pooled data.

### wcFst

Calculates Weir and Cockerham's Fst estimator bi-allelic genotype data (Weir and Cockerham 1984).  Sites with less than five genotypes in the target and background are skipped because they provide unreliable estimates of Fst.  Fix sites are also ignored.

```
INFO: help
INFO: description:
      wcFst is Weir & Cockerham's Fst for two populations.  Negative values are VALID,
      they are sites which can be treated as zero Fst. For more information see Evolution, Vol. 38 N. 6 Nov 1984.
      Specifically wcFst uses equations 1,2,3,4.

Output : 3 columns :
     1. seqid
     2. position
     3. target allele frequency
     4. background allele frequency
     5. wcFst

INFO: usage:  wcFst --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1 --type PL

INFO: required: t,target     -- argument: a zero based comma separated list of target individuals corrisponding to VCF columns
INFO: required: b,background -- argument: a zero based comma separated list of background individuals corrisponding to VCF columns
INFO: required: f,file       -- argument: proper formatted VCF
INFO: required, y,type       -- argument: genotype likelihood format; genotype : GT,GL,PL,GP
INFO: optional: r,region     -- argument: a tabix compliant genomic range: seqid or seqid:start-end
INFO: optional: d,deltaaf    -- argument: skip sites where the difference in allele frequencies is less than deltaaf, default is zero
```

### segmentFst

This program provides a way to find continious regions with high Fst values.  It takes the output of wcFst and produces a BED file.  These high Fst region can be permutated with 'permuteGPATwindow'.

```
INFO: help
INFO: description:
      Creates genomic segments (bed file) for regions with high wcFst
Output : 8 columns :
     1. Seqid
     2. Start (zero based)
     3. End   (zero based)
     4. Average Fst
     5. Average high Fst (Fst > -s)
     6. N Fst values in segment
     7. N high fst values in segment
     8. Segment length
INFO: usage:  segmentFst -s 0.7 -f wcFst.output.txt

INFO: required: -f            -- Output from wcFst
INFO: optional: -s            -- High Fst cutoff [0.8]

```


### popStats
Calculates basic population statistics at bi-allelic sites. The allele frequency is the number of non-reference alleles divided by the total number of alleles.  The expected hetrozygosity is 2*p*q, where p is the non-reference allele frequency and q is 1-p.  The observed heterozgosity is the fraction of 0/1 genotypes out of all genotypes.  The inbreeding coefficent, Fis, is the relative heterozygosity of each individual vs. compared to the target group.

```
INFO: help
INFO: description:
      General population genetic statistics for each SNP

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
INFO: usage:  popStat --type PL --target 0,1,2,3,4,5,6,7 --file my.vcf

INFO: required: t,target     -- a zero based comma separated list of target individuals corresponding to VCF columns
INFO: required: f,file       -- proper formatted VCF
INFO: required, y,type       -- genotype likelihood format; genotype : GL,PL,GP
INFO: optional, r,region     -- a tabix compliant region : chr1:1-1000 or chr1
```
### genotypeSummary

Generates a table of genotype counts.

```
INFO: help
INFO: description:
      Summarizes genotype counts for bi-allelic SNVs and indel

INFO: usage:  genotypeSummmary --type PL --target 0,1,2,3,4,5,6,7 --file my.vcf --snp

INFO: required: t,target     -- a zero based comma separated list of target individuals corresponding to VCF columns
INFO: required: f,file       -- proper formatted VCF
INFO: required, y,type       -- genotype likelihood format; genotype : GL,PL,GP
INFO: optional, r,region     -- a tabix compliant region : chr1:1-1000 or chr1
INFO: optional, s,snp        -- Only count SNPs

```

### pFst

pFst is a likelihood ratio test (LRT) quantifying allele frequency differences between populations.  The LRT by default uses the binomial distribution.  If Genotype likelihoods are provided it uses a modified binomial that weights each allele count by its certainty.  If type is set to 'PO' the LRT uses a beta distribution to fit the allele frequency spectrum of the target and background.  PO requires the AD and DP genotype fields and requires at least two pools for the target and background.  The p-value calculated in pFst is based on the chi-squared distribution with one degree of freedom.

```
INFO: help
INFO: description:
     pFst is a probabilistic approach for detecting differences in allele frequencies between two populations.

Output : 3 columns :
     1. seqid
     2. position
     3. pFst probability

INFO: usage:  pFst --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1 --type PL

INFO: required: t,target     -- argument: a zero based comma separated list of target individuals corresponding to VCF columns
INFO: required: b,background -- argument: a zero based comma separated list of background individuals corresponding to VCF columns
INFO: required: f,file       -- argument: a properly formatted VCF.
INFO: required: y,type       -- argument: genotype likelihood format ; genotypes: GP, GL or PL; pooled: PO
INFO: optional: d,deltaaf    -- argument: skip sites where the difference in allele frequencies is less than deltaaf, default is zero
INFO: optional: r,region     -- argument: a tabix compliant genomic range : seqid or seqid:start-end
INFO: optional: c,counts     -- switch  : use genotype counts rather than genotype likelihoods to estimate parameters, default false
```

### EHH and PI

The 'sequenceDiversity' program calculates extended haplotype homozygosity and pi within a fixed-width sliding window.  This requires phased data.

```
INFO: help
INFO: description:
      The sequenceDiversity program calculates two popular metrics of  haplotype diversity: pi and
      extended haplotype homozygoisty (eHH).  Pi is calculated using the Nei and Li 1979 formulation.
      eHH a convenient way to think about haplotype diversity.  When eHH = 0 all haplotypes in the window
      are unique and when eHH = 1 all haplotypes in the window are identical.

Output : 5 columns:
         1.  seqid
         2.  start of window
         3.  end of window
         4.  pi
         5.  eHH


INFO: usage: sequenceDiversity --target 0,1,2,3,4,5,6,7 --file my.vcf

INFO: required: t,target     -- argument: a zero base comma separated list of target individuals corresponding to VCF columns
INFO: required: f,file       -- argument: a properly formatted phased VCF file
INFO: required: y,type       -- argument: type of genotype likelihood: PL, GL or GP
INFO: optional: a,af         -- sites less than af  are filtered out; default is 0
INFO: optional: r,region     -- argument: a tabix compliant region : "seqid:0-100" or "seqid"
INFO: optional: w,window     -- argument: the number of SNPs per window; default is 20

```

### meltEHH

The program 'meltEHH' produces the data to generate the following plot:

<img src="https://github.com/vcflib/vcflib/blob/master/examples/example-ehh.png?raw=true" alt="" width=400>

```
INFO: help
INFO: description:
     meltEHH provides the data to plot EHH curves.
Output : 4 columns :
     1. seqid
     2. position
     3. EHH
     4. ref or alt [0 == ref]
Usage:
      meltEHH --target 0,1,2,3,4,5,6,7 --pos 10 --file my.phased.vcf  \
           --region chr1:1-1000 > STDOUT 2> STDERR

Params:
       required: t,target   <STRING>  A zero base comma separated list of target
                                     individuals corresponding to VCF columns
       required: r,region   <STRING>  A tabix compliant genomic range
                                     format: "seqid:start-end" or "seqid"
       required: f,file     <STRING>  Proper formatted and phased VCF.
       required: y,type     <STRING>  Genotype likelihood format: GT,PL,GL,GP
       required: p,position <INT>     Variant position to melt.
       optional: a,af       <DOUBLE>  Alternative alleles with frequencies less
                                     than [0.05] are skipped.
```
### iHS

iHS calculates the integrated haplotype score which measures the relative decay of extended haplotype homozygosity (EHH) for the reference and alternative alleles at a site (see: voight et al. 2006, Spiech & Hernandez 2014).  Our code is highly concordant with both implementations mentioned. However, we do not set an upper limit to the allele frequency.  iHS can be run without a genetic map, in which case the change in EHH is integrated over a constant.  Human genetic maps for GRCh36 and GRCh37 (hg18 & hg19) can be found at: http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ . iHS by default interpolates SNV positions to genetic position (you don't need a genetic position for every VCF entry in the map file).

iHS analyses requires normalization by allele frequency.  It is important that iHS is calculated over large regions so that the normalization does not down weight real signals.  For genome-wide runs it is recommended to run slightly overlapping windows and throwing out values that fail integration (columns 7 & 8 in the output) and then removing duplicates by using the 'sort' and 'uniq' linux commands.  Normalization of the output is as simple as running 'normalize-iHS'.

```
INFO: help
INFO: description:
     iHS calculates the integrated ratio of haplotype decay between the reference and non-reference allele.
Output : 4 columns :
     1. seqid
     2. position
     3. target allele frequency
     4. integrated EHH (alternative)
     5. integrated EHH (reference)
     6. iHS ln(iEHHalt/iEHHref)
     7. != 0 integration failure
     8. != 0 integration failure

Usage:
      iHS  --target 0,1,2,3,4,5,6,7 --file my.phased.vcf  \
           --region chr1:1-1000 > STDOUT 2> STDERR

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

```
### smoother
```
A method for window smoothing many of the GPAT++ formats.

INFO: help
INFO: description:
      Smoother averages a set of scores over a sliding genomic window.
      Smoother slides over genomic positions not the SNP indices. In other words
      the number of scores within a window will not be constant. The last
      window for each seqid can be smaller than the defined window size.
      Smoother automatically analyses different seqids separately.
Output : 4 columns :
     1. seqid
     2. window start
     2. window end
     3. averaged score

INFO: usage: smoother --format pFst --file GPA.output.txt

INFO: required: f,file     -- argument: a file created by GPAT++
INFO: required: o,format   -- argument: format of input file, case sensitive
                              available format options:
                                wcFst, pFst, bFst, iHS, xpEHH, abba-baba
INFO: optional: w,window   -- argument: size of genomic window in base pairs (default 5000)
INFO: optional: s,step     -- argument: window step size in base pairs (default 1000)
INFO: optional: t,truncate -- flag    : end last window at last position (zero based) last window at last position (zero based)
```

# Development

## build from source

VCFLIB uses the cmake build system, after a recursive checkout
of the sources make the files in the ./build directory with:

```sh
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
mkdir -p build && cd build
cmake ..
cmake --build .
cmake --install .
```

and to run the tests

```sh
ctest --verbose
```

Executables are built into the `./build` directory in the repository.

Build dependencies can be viewed in the Travis-CI and github-CI
scripts (see badges above), as well as [guix.scm](./guix.scm) used by
us to create the build environment. Essentially:

- C++ compiler
- htslib
- tabixpp

For include files add

- libhts-dev
- libtabixpp-dev
- libtabixpp0

And for some of the VCF executables

- python
- perl

## link library

```make
# Put this in your Makefile

VCFLIB_DIR = /path/to/vcflib
VCFLIB_INC = -I $(VCFLIB_DIR)/include -I $(VCFLIB_DIR)/tabixpp/htslib
VCFLIB_LIB = -L $(VCFLIB_DIR)/tabixpp -L $(VCFLIB_DIR)/lib -lhts -lvcflib -lz -lm -llzma -lbz2
CXX="g++"
CXXFLAGS="-std=-std=c++0x -Ofast -D_FILE_OFFSET_BITS=64"

mytool: mytool.cpp
	$(CXX) $(CXXFLAGS) $(VCFLIB_INC) -o $@ $< $(VCFLIB_LIB)
```

## source code

See [vcfecho.cpp](./src/vcfecho.cpp) for basic usage.
[Variant.h](./src/Variant.h) and [Variant.cpp](./src/Variant.cpp)
describe methods available in the API.  vcflib is incorporated into
several projects, such as
[freebayes](https://github.com/freebayes/freebayes), which may provide
a point of reference for prospective developers.  Note vcflib contains
submodules (git repositories) comprising some dependencies. A full
Guix development environment we use is defined [here](./guix.scm).


# LICENSE

This software is distributed under the free software [MIT
LICENSE](./LICENSE).
