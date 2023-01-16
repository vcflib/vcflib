# vcflib

### A C++ library for parsing and manipulating VCF files.

[![Github-CI](https://github.com/vcflib/vcflib/workflows/CI/badge.svg)](https://github.com/vcflib/vcflib/actions?query=workflow%3ACI) [![AnacondaBadge](https://anaconda.org/bioconda/vcflib/badges/version.svg)](https://anaconda.org/bioconda/vcflib) [![DL](https://anaconda.org/bioconda/vcflib/badges/downloads.svg)](https://anaconda.org/bioconda/vcflib) [![BrewBadge](https://img.shields.io/badge/%F0%9F%8D%BAbrew-vcflib-brightgreen.svg)](https://github.com/brewsci/homebrew-bio) [![GuixBadge](https://img.shields.io/badge/gnuguix-vcflib-brightgreen.svg)](https://www.gnu.org/software/guix/packages/V/) [![DebianBadge](https://badges.debian.net/badges/debian/testing/libvcflib-dev/version.svg)](https://packages.debian.org/testing/libvcflib-dev) [![C++0x](https://img.shields.io/badge/Language-C++0x-steelblue.svg)](https://www.cprogramming.com/c++11/what-is-c++0x.html) [![Chat on Matrix](https://matrix.to/img/matrix-badge.svg)](https://matrix.to/#/#vcflib:matrix.org)

Vcflib and related tools are the workhorses in bioinformatics for processing the VCF variant calling format. See

[A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009123).
Garrison E, Kronenberg ZN, Dawson ET, Pedersen BS, Prins P (2022), PLoS Comput Biol 18(5): e1009123. https://doi.org/10.1371/journal.pcbi.1009123

(bibtex reference at the end of this page)

## news

![Humpty Logo - CC0 license](./logos/humpty-dumpty.jpg)

**January 2023**: this is vcflib's first *Humpty Dumpty* release: [vcfcreatemulti](./doc/vcfcreatemulti.md) is the natural companion to [vcfwave](./doc/vcfwave.md).
Often variant callers are not perfect.
**vcfwave** with its companion tool **vcfcreatemulti** can take an existing VCF file that contains multiple complex overlapping and even nested alleles and, unlike Humpty Dumpty, take them apart and put them together again.
Thereby, hopefully, creating sane VCF output that is useful for analysis and getting rid of false positives.

We created these tools by including the state-of-the-art [biWFA](https://github.com/smarco/WFA2-lib) wavefront aligner.
The tools are particularly useful for the output from structural variation callers and pangenome genotypers, such as used by the Human Pangenome Reference Consortium (HPRC) because of overlapping ALT segments.

See also [RELEASE_NOTES.md](./RELEASE_NOTES.md)

**May 2022**: the [vcflib paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009123) has been published on PLoS Computational Biology!

**April 2022**: vcflib has just gone pangenome!

By introducing the wavefront algorithm we can now realign long sequences and reduce call complexity (and FPs!) introduced by pangenome variant callers using the new [vcfwave](./doc/vcfwave.md) tool.

See also [RELEASE_NOTES.md](./RELEASE_NOTES.md)

## overview

The [Variant Call Format
(VCF)](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)
is a flat-file, tab-delimited textual format that
describes reference-indexed variations between individuals.  VCF
provides a common interchange format for the description of variation
in individuals and populations of samples, and has become the
*de facto* standard reporting format for a wide array of genomic
variant detectors.

vcflib provides methods to manipulate and interpret sequence variation
described by VCF.  It is both:

 * an API for parsing and operating on records of genomic variation as
   it can be described by the VCF format
 * a collection of command-line utilities for executing complex
   manipulations on VCF files

vclib is both a library (with an API) and a collection of useful
tools. The API provides a quick and extremely permissive method to
read and write VCF files.  Extensions and applications of the library
provided in the included utilities (*.cpp) comprise the vast bulk of
the library's utility.

We have also added infrastructure to write Python bindings. See [below](#python).

---

Short index:

- [Install](#INSTALL)
- [Usage](#USAGE)
- [TOOLS](#TOOLS)
  * [Filter](#filter)
  * [Metrics](#metrics)
  * [Phenotype](#phenotype)
  * [Genotype](#genotype)
  * [Transformation](#transformation)
  * [Statistics](#statistics)
  * [Scripts](#scripts)
  * [Python bindings](#python)
- [Link library](#link-library)
- [Build from source](#build-from-source)
- [Development](#Development)
- [LICENSE](#LICENSE)
- [Credit work](#Credit)
  * [Bibtex reference](#Bibtex-reference)

---

## INSTALL

For latest updates see [RELEASE NOTES](./RELEASE_NOTES.md).

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

We develop against guix and vcflib is packaged as

```sh
guix package -i vcflib
```

See also the Guix shell below.

## USAGE

Users are encouraged to drive the utilities in the library in a
streaming fashion, using Unix pipes to fully utilize resources on
multi-core systems.  Piping provides a convenient method to interface
with other libraries (vcf-tools, BedTools, GATK, htslib,
[bio-vcf](https://github.com/vcflib/bio-vcf), bcftools,
[freebayes](https://github.com/freebayes)) which interface via VCF
files, allowing the composition of an immense variety of processing
functions. Examples can be found in the scripts,
e.g. [./scripts](./scripts/).


# TOOLS

<!--

    >>> from pytest.rtest import run_stdout, head, cat

-->


<!--
  Created with ./scripts/bin2md.rb --index
-->


## filter

| filter command | description |
| :-------------- | :---------- |
 | [vcfuniq](./doc/vcfuniq.md) | List unique genotypes. Similar to GNU uniq, but aimed at VCF records. **vcfuniq** removes records which have the same position, ref, and alt as the previous record on a sorted VCF file. Note that it does not adjust/combine genotypes in the output, but simply takes the first record. See also vcfcreatemulti for combining records. |
 | [vcfuniqalleles](./doc/vcfuniqalleles.md) | List unique alleles For each record, remove any duplicate alternate alleles that may have resulted from merging separate VCF files. |
 | [vcffilter](./doc/vcffilter.md) | VCF filter the specified vcf file using the set of filters |

## metrics

| metrics command | description |
| :-------------- | :---------- |
 | [vcfhethomratio](./doc/vcfhethomratio.md) | Generates the het/hom ratio for each individual in the file |
 | [vcfhetcount](./doc/vcfhetcount.md) | Calculate the heterozygosity rate: count the number of alternate alleles in heterozygous genotypes in all records in the vcf file |
 | [vcfdistance](./doc/vcfdistance.md) | Adds a tag to each variant record which indicates the distance to the nearest variant. (defaults to BasesToClosestVariant if no custom tag name is given. |
 | [vcfentropy](./doc/vcfentropy.md) | Annotate VCF records with the Shannon entropy of flanking sequence. Anotates the output VCF file with, for each record, EntropyLeft, EntropyRight, EntropyCenter, which are the entropies of the sequence of the given window size to the left, right, and center of the record. Also adds EntropyRef and EntropyAlt for each alt. |
 | [vcfcheck](./doc/vcfcheck.md) | Validate integrity and identity of the VCF by verifying that the VCF record's REF matches a given reference file. |

## phenotype

| phenotype command | description |
| :-------------- | :---------- |
 | [permuteGPAT++](./doc/permuteGPAT++.md) | **permuteGPAT++** is a method for adding empirical p-values to a GPAT++ score. |

## genotype

| genotype command | description |
| :-------------- | :---------- |
 | [hapLrt](./doc/hapLrt.md) | HapLRT is a likelihood ratio test for haplotype lengths. The lengths are modeled with an exponential distribution. The sign denotes if the target has longer haplotypes (1) or the background (-1). |
 | [normalize-iHS](./doc/normalize-iHS.md) | normalizes iHS or XP-EHH scores. |
 | [abba-baba](./doc/abba-baba.md) | **abba-baba** calculates the tree pattern for four indviduals. This tool assumes reference is ancestral and ignores non **abba-baba** sites. The output is a boolian value: 1 = true , 0 = false for abba and baba. the tree argument should be specified from the most basal taxa to the most derived. |

## transformation

| transformation command | description |
| :-------------- | :---------- |
 | [vcfintersect](./doc/vcfintersect.md) | VCF set analysis |
 | [vcfleftalign](./doc/vcfleftalign.md) | Left-align indels and complex variants in the input using a pairwise ref/alt alignment followed by a heuristic, iterative left realignment process that shifts indel representations to their absolute leftmost (5') extent. |
 | [vcfglxgt](./doc/vcfglxgt.md) | Set genotypes using the maximum genotype likelihood for each sample. |
 | [vcfkeepsamples](./doc/vcfkeepsamples.md) | outputs each record in the vcf file, removing samples not listed on the command line |
 | [vcfallelicprimitives](./doc/vcfallelicprimitives.md) | Note that this tool is considered legacy and will emit a warning! Use [vcfwave](./doc/vcfwave.md) instead. |
 | [vcfcommonsamples](./doc/vcfcommonsamples.md) | Generates each record in the first file, removing samples not present in the second |
 | [vcfcreatemulti](./doc/vcfcreatemulti.md) | Go through sorted VCF and if overlapping alleles are represented across multiple records, merge them into a single record. Currently only for indels. |
 | [vcfld](./doc/vcfld.md) | Compute LD |
 | [vcfsamplenames](./doc/vcfsamplenames.md) | List sample names |
 | [vcfsample2info](./doc/vcfsample2info.md) | Take annotations given in the per-sample fields and add the mean, median, min, or max to the site-level INFO. |
 | [vcflength](./doc/vcflength.md) | Add length info field |
 | [vcfstreamsort](./doc/vcfstreamsort.md) | Sorts the input (either stdin or file) using a streaming sort algorithm. Guarantees that the positional order is correct provided out-of-order variants are no more than 100 positions in the VCF file apart. |
 | [dumpContigsFromHeader](./doc/dumpContigsFromHeader.md) | Dump contigs from header |
 | [vcfgeno2haplo](./doc/vcfgeno2haplo.md) | Convert genotype-based phased alleles within --window-size into haplotype alleles. Will break haplotype construction when encountering non-phased genotypes on input. |
 | [vcfinfo2qual](./doc/vcfinfo2qual.md) | Sets QUAL from info field tag keyed by [key]. The VCF file may be omitted and read from stdin. The average of the field is used if it contains multiple values. |
 | [vcfevenregions](./doc/vcfevenregions.md) | Generates a list of regions, e.g. chr20:10..30 using the variant density information provided in the VCF file to ensure that the regions have even numbers of variants. This can be use to reduce the variance in runtime when dividing variant detection or genotyping by genomic coordinates. |
 | [vcfbreakmulti](./doc/vcfbreakmulti.md) | If multiple alleles are specified in a single record, break the record into multiple lines, preserving allele-specific INFO fields. |
 | [vcfwave](./doc/vcfwave.md) | Realign reference and alternate alleles with WFA, parsing out the primitive alleles into multiple VCF records. New records have IDs that reference the source record ID. Genotypes are handled. Deletions generate haploid/missing genotypes at overlapping sites. |
 | [vcfafpath](./doc/vcfafpath.md) | Display genotype paths |
 | [vcfremoveaberrantgenotypes](./doc/vcfremoveaberrantgenotypes.md) | strips samples which are homozygous but have observations implying heterozygosity. Remove samples for which the reported genotype (GT) and observation counts disagree (AO, RO). |
 | [vcf2tsv](./doc/vcf2tsv.md) | Converts VCF to per-allelle or per-genotype tab-delimited format, using null string to replace empty values in the table. Specifying -g will output one line per sample with genotype information. When there is more than one alt allele there will be multiple rows, one for each allele and, the info will match the 'A' index |
 | [vcfqual2info](./doc/vcfqual2info.md) | Puts QUAL into an info field tag keyed by [key]. |
 | [vcfnumalt](./doc/vcfnumalt.md) | outputs a VCF stream where NUMALT has been generated for each record using sample genotypes |
 | [vcfprimers](./doc/vcfprimers.md) | For each VCF record, extract the flanking sequences, and write them to stdout as FASTA records suitable for alignment. |
 | [vcf2fasta](./doc/vcf2fasta.md) | Generates sample_seq:N.fa for each sample, reference sequence, and chromosomal copy N in [0,1... ploidy]. Each sequence in the fasta file is named using the same pattern used for the file name, allowing them to be combined. |
 | [vcffixup](./doc/vcffixup.md) | Generates a VCF stream where AC and NS have been generated for each record using sample genotypes |
 | [vcfkeepinfo](./doc/vcfkeepinfo.md) | To decrease file size remove INFO fields not listed on the command line |
 | [vcfsamplediff](./doc/vcfsamplediff.md) | Establish putative somatic variants using reported differences between germline and somatic samples. Tags each record where the listed sample genotypes differ with <tag>. The first sample is assumed to be germline, the second somatic. Each record is tagged with <tag>={germline,somatic,loh} to specify the type of variant given the genotype difference between the two samples. |
 | [vcfflatten](./doc/vcfflatten.md) | Removes multi-allelic sites by picking the most common alternate. Requires allele frequency specification 'AF' and use of 'G' and 'A' to specify the fields which vary according to the Allele or Genotype. VCF file may be specified on the command line or piped as stdin. |
 | [vcfglbound](./doc/vcfglbound.md) | Adjust GLs so that the maximum GL is 0 by dividing all GLs for each sample by the max. |
 | [vcfecho](./doc/vcfecho.md) | Echo VCF to stdout (simple demo) |
 | [vcfgeno2alleles](./doc/vcfgeno2alleles.md) | modifies the genotypes field to provide the literal alleles rather than indexes |
 | [vcf2dag](./doc/vcf2dag.md) | Modify VCF to be able to build a directed acyclic graph (DAG) |
 | [vcfannotategenotypes](./doc/vcfannotategenotypes.md) | Examine genotype correspondence. Annotate genotypes in the first file with genotypes in the second adding the genotype as another flag to each sample filed in the first file. annotation-tag is the name of the sample flag which is added to store the annotation. also adds a 'has_variant' flag for sites where the second file has a variant. |
 | [vcfremap](./doc/vcfremap.md) | For each alternate allele, attempt to realign against the reference with lowered gap open penalty. If realignment is possible, adjust the cigar and reference/alternate alleles. Observe how different alignment parameters, including context and entropy-dependent ones, influence variant classification and interpretation. |
 | [vcfremovesamples](./doc/vcfremovesamples.md) | outputs each record in the vcf file, removing samples listed on the command line |
 | [vcfcleancomplex](./doc/vcfcleancomplex.md) | Removes reference-matching sequence from complex alleles and adjusts records to reflect positional change. |
 | [vcfindex](./doc/vcfindex.md) | Adds an index number to the INFO field (id=position) |
 | [vcfannotate](./doc/vcfannotate.md) | Intersect the records in the VCF file with targets provided in a BED file. Intersections are done on the reference sequences in the VCF file. If no VCF filename is specified on the command line (last argument) the VCF read from stdin. |
 | [smoother](./doc/smoother.md) | smoothes is a method for window smoothing many of the GPAT++ formats. |
 | [vcfkeepgeno](./doc/vcfkeepgeno.md) | Reduce file size by removing FORMAT fields not listed on the command line from sample specifications in the output |
 | [vcfcombine](./doc/vcfcombine.md) | Combine VCF files positionally, combining samples when sites and alleles are identical. Any number of VCF files may be combined. The INFO field and other columns are taken from one of the files which are combined when records in multiple files match. Alleles must have identical ordering to be combined into one record. If they do not, multiple records will be emitted. |
 | [vcfcat](./doc/vcfcat.md) | Concatenates VCF files |
 | [vcfinfosummarize](./doc/vcfinfosummarize.md) | Take annotations given in the per-sample fields and add the mean, median, min, or max to the site-level INFO. |
 | [vcfclassify](./doc/vcfclassify.md) | Creates a new VCF where each variant is tagged by allele class: snp, ts/tv, indel, mnp |
 | [vcfoverlay](./doc/vcfoverlay.md) | Overlay records in the input vcf files with order as precedence. |
 | [vcfaddinfo](./doc/vcfaddinfo.md) | Adds info fields from the second file which are not present in the first vcf file. |
 | [vcfnullgenofields](./doc/vcfnullgenofields.md) | Makes the FORMAT for each variant line the same (uses all the FORMAT fields described in the header). Fills out per-sample fields to match FORMAT. Expands GT values of '.' with number of alleles based on ploidy (eg: './doc/.' for dipolid). |
 | [vcfgenosamplenames](./doc/vcfgenosamplenames.md) | Get samplenames |

## statistics

| statistics command | description |
| :-------------- | :---------- |
 | [permuteSmooth](./doc/permuteSmooth.md) | **permuteSmooth** is a method for adding empirical p-values smoothed wcFst scores. |
 | [plotHaps](./doc/plotHaps.md) | **plotHaps** provides the formatted output that can be used with 'bin/plotHaplotypes.R'. |
 | [pVst](./doc/pVst.md) | **pVst** calculates vst, a measure of CNV stratification. |
 | [wcFst](./doc/wcFst.md) | **wcFst** is Weir & Cockerham's Fst for two populations. Negative values are VALID, they are sites which can be treated as zero Fst. For more information see Evolution, Vol. 38 N. 6 Nov 1984. Specifically **wcFst** uses equations 1,2,3,4. |
 | [popStats](./doc/popStats.md) | General population genetic statistics for each SNP |
 | [vcfsitesummarize](./doc/vcfsitesummarize.md) | Summarize by site |
 | [vcfrandomsample](./doc/vcfrandomsample.md) | Randomly sample sites from an input VCF file, which may be provided as stdin. Scale the sampling probability by the field specified in KEY. This may be used to provide uniform sampling across allele frequencies, for instance. |
 | [vcfroc](./doc/vcfroc.md) | Generates a pseudo-ROC curve using sensitivity and specificity estimated against a putative truth set. Thresholding is provided by successive QUAL cutoffs. |
 | [vcfstats](./doc/vcfstats.md) | Prints statistics about variants in the input VCF file. |
 | [bFst](./doc/bFst.md) | **bFst** is a Bayesian approach to Fst. Importantly **bFst** accounts for genotype uncertainty in the model using genotype likelihoods. For a more detailed description see: `A Bayesian approach to inferring population structure from dominant markers' by Holsinger et al. Molecular Ecology Vol 11, issue 7 2002. The likelihood function has been modified to use genotype likelihoods provided by variant callers. There are five free parameters estimated in the model: each subpopulation's allele frequency and Fis (fixation index, within each subpopulation), a free parameter for the total population's allele frequency, and Fst. |
 | [vcfcountalleles](./doc/vcfcountalleles.md) | Count alleles |
 | [vcfparsealts](./doc/vcfparsealts.md) | Alternate allele parsing method. This method uses pairwise alignment of REF and ALTs to determine component allelic primitives for each alternate allele. |
 | [genotypeSummary](./doc/genotypeSummary.md) | Generates a table of genotype counts. Summarizes genotype counts for bi-allelic SNVs and indel |
 | [vcfgenosummarize](./doc/vcfgenosummarize.md) | Adds summary statistics to each record summarizing qualities reported in called genotypes. Uses: RO (reference observation count), QR (quality sum reference observations) AO (alternate observation count), QA (quality sum alternate observations) |
 | [vcfaltcount](./doc/vcfaltcount.md) | count the number of alternate alleles in all records in the vcf file |
 | [meltEHH](./doc/meltEHH.md) |  |
 | [iHS](./doc/iHS.md) | **iHS** calculates the integrated haplotype score which measures the relative decay of extended haplotype homozygosity (EHH) for the reference and alternative alleles at a site (see: voight et al. 2006, Spiech & Hernandez 2014). |
 | [segmentFst](./doc/segmentFst.md) | **segmentFst** creates genomic segments (bed file) for regions with high wcFst |
 | [pFst](./doc/pFst.md) | **pFst** is a probabilistic approach for detecting differences in allele frequencies between two populations. |
 | [sequenceDiversity](./doc/sequenceDiversity.md) | The **sequenceDiversity** program calculates two popular metrics of haplotype diversity: pi and extended haplotype homozygoisty (eHH). Pi is calculated using the Nei and Li 1979 formulation. eHH a convenient way to think about haplotype diversity. When eHH = 0 all haplotypes in the window are unique and when eHH = 1 all haplotypes in the window are identical. |
 | [vcfgenotypecompare](./doc/vcfgenotypecompare.md) | adds statistics to the INFO field of the vcf file describing the amount of discrepancy between the genotypes (GT) in the vcf file and the genotypes reported in the <other-genotype-tag>. use this after vcfannotategenotypes to get correspondence statistics for two vcfs. |
 | [vcfrandom](./doc/vcfrandom.md) | Generate a random VCF file |
 | [segmentIhs](./doc/segmentIhs.md) | Creates genomic segments (bed file) for regions with high wcFst |
 | [vcfgenotypes](./doc/vcfgenotypes.md) | Report the genotypes for each sample, for each variant in the VCF. Convert the numerical represenation of genotypes provided by the GT field to a human-readable genotype format. |

See also [vcflib.md](./doc/vcflib.md).

## scripts

The vcflib source repository contains a number of additional scripts.
Click on the link to see the source code.
| script | description |
| :-------------- | :---------- |
| [vcfclearinfo](./scripts/vcfclearinfo) | clear INFO field |
| [vcfqualfilter](./scripts/vcfqualfilter) | quality filter |
| [vcfnulldotslashdot](./scripts/vcfnulldotslashdot) | rewrite null genotypes to ./. |
| [vcfprintaltdiscrepancy.r](./scripts/vcfprintaltdiscrepancy.r) | show ALT discrepancies in a table |
| [vcfremovenonATGC](./scripts/vcfremovenonATGC) | remove non-nucleotides in REF or ALT |
| [plotSmoothed.R](./scripts/plotSmoothed.R) | smooth plot of wcFst, pFst or abba-baba |
| [vcf_strip_extra_headers](./scripts/vcf_strip_extra_headers) | strip headers |
| [plotHapLrt.R](./scripts/plotHapLrt.R) | plot results of pFst |
| [vcfbiallelic](./scripts/vcfbiallelic) | remove anything that is not biallelic |
| [vcfsort](./scripts/vcfsort) | sort VCF using shell script |
| [vcfnosnps](./scripts/vcfnosnps) | remove SNPs |
| [vcfmultiwayscripts](./scripts/vcfmultiwayscripts) | more multiway comparisons |
| [vcfgtcompare.sh](./scripts/vcfgtcompare.sh) | annotates records in the first file with genotypes and sites from the second |
| [plotPfst.R](./scripts/plotPfst.R) | plot pFst |
| [vcfregionreduce_and_cut](./scripts/vcfregionreduce_and_cut) | reduce, gzip, and tabix |
| [plotBfst.R](./scripts/plotBfst.R) | plot results of pFst |
| [vcfnobiallelicsnps](./scripts/vcfnobiallelicsnps) | remove biallelic SNPs |
| [vcfindels](./scripts/vcfindels) | show INDELS |
| [vcfmultiway](./scripts/vcfmultiway) | multiway comparison |
| [vcfregionreduce](./scripts/vcfregionreduce) | reduce VCFs using a BED File, gzip them up and create tabix index |
| [vcfprintaltdiscrepancy.sh](./scripts/vcfprintaltdiscrepancy.sh) | runner |
| [vcfclearid](./scripts/vcfclearid) | clear ID field |
| [vcfcomplex](./scripts/vcfcomplex) | remove all SNPs but keep SVs |
| [vcffirstheader](./scripts/vcffirstheader) | show first header |
| [plotXPEHH.R](./scripts/plotXPEHH.R) | plot XPEHH |
| [vcfregionreduce_pipe](./scripts/vcfregionreduce_pipe) | reduce, gzip and tabix in a pipe |
| [vcfplotaltdiscrepancy.sh](./scripts/vcfplotaltdiscrepancy.sh) | plot ALT discrepancy runner |
| [vcfplottstv.sh](./scripts/vcfplottstv.sh) | runner |
| [vcfnoindels](./scripts/vcfnoindels) | remove INDELs |
| [bgziptabix](./scripts/bgziptabix) | runs bgzip on the input and tabix indexes the result |
| [plotHaplotypes.R](./scripts/plotHaplotypes.R) | plot results |
| [vcfplotsitediscrepancy.r](./scripts/vcfplotsitediscrepancy.r) | plot site discrepancy |
| [vcfindelproximity](./scripts/vcfindelproximity) | show SNPs around an INDEL |
| [bed2region](./scripts/bed2region) | convert VCF CHROM column in VCF file to region |
| [vcfplotaltdiscrepancy.r](./scripts/vcfplotaltdiscrepancy.r) | plot ALT discrepancies |
| [plot_roc.r](./scripts/plot_roc.r) | plot ROC |
| [vcfmultiallelic](./scripts/vcfmultiallelic) | remove anything that is not multiallelic |
| [vcfsnps](./scripts/vcfsnps) | show SNPs |
| [vcfvarstats](./scripts/vcfvarstats) | use fastahack to get stats |
| [vcfregionreduce_uncompressed](./scripts/vcfregionreduce_uncompressed) | reduce, gzip and tabix |
| [plotWCfst.R](./scripts/plotWCfst.R) | plot wcFst |
| [vcf2bed.py](./scripts/vcf2bed.py) | transform VCF to BED file |
| [vcfjoincalls](./scripts/vcfjoincalls) | overlay files using QUAL and GT from a second VCF |
| [vcf2sqlite.py](./scripts/vcf2sqlite.py) | push VCF file into SQLite3 database using dbname |


## python

vcflib has rudimentary python bindings, but the are easy to build up on. See [pyvcflib](./test/pytest/pyvcflib.md).

# Development

## build from source

VCFLIB uses the cmake build system, after a recursive checkout of the sources make the files in the 'build' directory with:

```sh
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
mkdir -p build && cd build
cmake  -DCMAKE_BUILD_TYPE=Debug -DZIG=OFF -DOPENMP=OFF ..
cmake --build .
cmake --install .
```

and to run the tests

```sh
ctest --verbose
```

Executables are built into the `./build` directory in the repository.

Note, if you have an existing repo update submodules with

```sh
git submodule update --init --recursive --progress
cd build
cmake --build . --target clean
```

Build dependencies can be viewed in the github-CI
scripts (see badges above), as well as [guix.scm](./guix.scm) used by
us to create the build environment (for instructions see the header of
guix.scm). Essentially:

- cmake
- C++ compiler
- htslib
- tabixpp
- WFA2
- pybind11 (for testing)

For include files add

- libhts-dev
- libtabixpp-dev
- libtabixpp0

And for some of the VCF executables

- zig 0.9.1 (disable with cmake -DZIG=OFF)
- python
- perl

Make sure the `zig` tool is in the `PATH`.
The `zig` dependency is for the recently upgraded [vcfcreatemulti](./doc/vcfcreatemulti.md) that can be run after vcfwave to combine overlapping records.
More on `zig` can be found in the souce code [README](./src/zig/README.md).

### Using a different htslib

Check out htslib in tabixpp (recursively) and

    cmake -DHTSLIB_LOCAL:STRING=./htslib/ ..
    cmake --build .

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

# adding tests

vcflib uses different test systems. The most important one is the
[doctest](https://docs.python.org/3/library/doctest.html) because it
doubles as documentation. For an example see
[vcf2tsv.md](./test/pytest/vcf2tsv.md) which can be run from the
command line with

```sh
cd test
python3 -m doctest -o NORMALIZE_WHITESPACE -o REPORT_UDIFF pytest/vcf2tsv.md
```

We also added support for python bindings and unit tests. See [realign.py](./test/tests/realign.py) for an example.

# Support

The developers are on the vcflib
[matrix channel](https://matrix.to/#/#vcflib:matrix.org). Please do not use the github issue tracker for support issues!

# Contributing

To contribute code to vcflib send a github pull request. We may ask
you to add a working test case as described in 'adding tests'.

## LICENSE

This software is distributed under the free software [MIT LICENSE](./LICENSE).

## CREDIT

Citations are the bread and butter of Science.  If you are using this
software in your research and want to support our future work, please
cite the following publication:

Please cite:

[A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009123).
Garrison E, Kronenberg ZN, Dawson ET, Pedersen BS, Prins P (2022), PLoS Comput Biol 18(5): e1009123. https://doi.org/10.1371/journal.pcbi.1009123

## Bibtex reference

```bibtex
@article{10.1371/journal.pcbi.1009123,
    doi = {10.1371/journal.pcbi.1009123},
    author = {Garrison, Erik AND Kronenberg, Zev N. AND Dawson, Eric T. AND Pedersen, Brent S. AND Prins, Pjotr},
    journal = {PLOS Computational Biology},
    publisher = {Public Library of Science},
    title = {A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar},
    year = {2022},
    month = {05},
    volume = {18},
    url = {https://doi.org/10.1371/journal.pcbi.1009123},
    pages = {1-15}
}
```

Below the prepublished version of our paper

```bibtex
@article {Garrison2021.05.21.445151,
	author = {Garrison, Erik and Kronenberg, Zev N. and Dawson, Eric T. and Pedersen, Brent S. and Prins, Pjotr},
	title = {Vcflib and tools for processing the VCF variant call format},
	elocation-id = {2021.05.21.445151},
	year = {2021},
	doi = {10.1101/2021.05.21.445151},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/05/23/2021.05.21.445151},
	eprint = {https://www.biorxiv.org/content/early/2021/05/23/2021.05.21.445151.full.pdf},
	journal = {bioRxiv}
}
```
