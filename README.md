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

See [vcflib.md](./doc/vcflib.md).

### smoother

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

The standard build creates `build/vcflib.a`. Take a hint from the
[cmake](./CMakeLists.txt) file that builds all the vcflib tools.

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
