For contributions see
[contributors](https://github.com/vcflib/vcflib/graphs/contributors)
and
[commits](https://github.com/vcflib/vcflib/commits/master).

## TODO

- [ ] allow Debian to use packaged wfa2 lib
- [ ] vcfcreatemulti: fix problem with slow and wrong complex regions
- [ ] complete vcfcreatemulti merge multiple rows
      + [ ] check for indels which are really the same
      + [ ] combine vcfwave duplicated functionality
- bgzip
- tabix -p vcf my_file.vcf.gz
- pangenie, vg deconstruct, vcfbub

## ChangeLog v1.0.6 (2023xxxx)

Vcflib maintenance release

+ Fixed zig complaining about leaking memory

## ChangeLog v1.0.5 (20230116)

Vcflib's first *Humpty Dumpty* release: [vcfcreatemulti](./doc/vcfcreatemulti.md) is the natural companion to [vcfwave](./doc/vcfwave.md).

Often variant callers are not perfect.
**vcfwave** with its companion tool **vcfcreatemulti** can take an existing VCF file that contains multiple complex overlapping and even nested alleles and, unlike Humpty Dumpty, take them apart and put them together again.
Thereby, hopefully, creating sane VCF output that is useful for analysis and getting rid of false positives.

We created these tools by including the state-of-the-art [biWFA](https://github.com/smarco/WFA2-lib) wavefront aligner.
The tools are particularly useful for the output from structural variation callers and pangenome genotypers, such as used by the Human Pangenome Reference Consortium (HPRC) because of overlapping ALT segments.

Important changes:

+ vcfwave is introduced and vcfallelicprimitimes is now considered obsolete
+ INFO fields output order is now the same with every tool as on input parsing
+ vcfwave check merging of genotypes - write tests
+ vcfwave recompute AC, AFs from merged record
+ Added python bindings with pybind11
+ introduced the zig compiler with vcfcreatemulti.cpp as a first target (use cmake ZIG=OFF to disable). At this point the zig version (-n switch) gives identical results.

Introduction of O(n) wavefront algorithm WF to replace O(n^2) Smith-Waterman SW. Note that the output is different from the original SW implementation. SW is still optionally available but considered obsolete. Use the bi-directional vcfwave instead of vcfallelicprimitives.

+ Added realignment using the wavefront algorithm (now the default). See [vcfwave](./doc/vcfwave.md) (thank you Erik Garrison https://github.com/ekg and Santiago Marco-Sola  https://github.com/smarco).
+ WFA2-lib fix is merged upstream. Fixes bleeding in of macros https://github.com/vcflib/vcflib/issues/359
+ Support longer read inversions in vcfwave!
+ vcfallelicprimitives now considered legacy/obsolete
+ Fixed allowing for different field order - see https://github.com/vcflib/vcflib/issues/365
+ Improved CMake configuration
+ vcflib compiles with both gcc and clang++ and tests pass, see [guix-clang.scm](./guix-clang.scm) - mind that git submodules such as WFA2-lib still override to gcc
+ Fixed local build for tabixpp+htslib - note that htslib should be an upstream released version (currently 1.15.1). Unfortunately git submodule does not handle tags.
+ Fix for -L switch for vcfallelicprimitives
+ Added libasan and lto support
+ Removed useless googletest submodule
+ Moved git submodules into ./contrib
+ Added python testing framework
+ Added tabixpp back in as a submodule, fixes https://github.com/vcflib/vcflib/issues/305
+ Optimizations and bug fixes. (thanks @mphschmitt)
+ vcfcreatemulti merge multiple rows
+ rewrite of vcfcreatemulti using zig
+ vcfcreatemulti merge genotypes correctly, with tests
+ vcfcreatemulti adjust info and genotypes for variants that have multiple alts already (now errors)
+ vcfcreatemulti handle phase
+ vcfcreatemulti document building with zig
+ vcfcreatemulti added progress bar to vcfwave and vcfcreatemulti with update to tabixpp
+ vcfcreatemulti default vcfwave and vcfcreatemulti to nextgen mode
+ vcfcreatemulti check file is sorted for vcfcreatemulti and improve suggestions

## ChangeLog v1.0.4

Never properly released. Merged with v1.0.5.

## ChangeLog v1.0.3 (20220122)

This is a maintenance release of vcflib.

+ Merge intervaltree changes (thanks @jnunn and @timmassingham)
+ Built with gcc-11
+ Fix issue #251 hapLrt: fix segfault when accessing genotype field. (thanks @mphschmitt)
+ Fix vcfflatten: fix segfault when no 'AF' field is present (#47, thanks @mphschmitt)
+ Fixes on vcfnulldotslashdot #310 (thanks @WinterFor)
+ Fix issue #301: Replace raw pointer usage with std::unique_ptr #306 (thanks @Glebanister)
+ Fix man page installation #321 (thanks @alexreg)
+ Use `guix shell` instead of `guix environment` for development
+ Regenerated online docs
+ README: add matrix badge (removed gitter badge)

## ChangeLog v1.0.2 (20210104)

This is a maintenance release of vcflib, mostly improving the build
system, CI and generating markdown docs as well as man pages.

+ Removed tabixpp and htslib source dependencies, i.e., we are now using
  the distro provided libraries and include files through pkg-config.
  See also the [README](README.md#build-from-source)
+ Removed the tabixpp+htslib git submodules
+ Generalise and document the cmake build system
+ Added tests to the cmake build system and build instructions to README
+ Added support for ARM64 and PowerPC, see #292 (thanks @genisysram and @mr-c)
+ Added github actions for the issue tracker
+ Added githum CI
+ Updated header files in src with copyright/license info, see #16
+ Created markdown [docs](./doc/vcflib.md) and [man pages](./man/) for
  all utilities. Created a script bin2md for markdown generation and
  use pandoc for the man page generation.

## Older changes

For older changes view the git [log](https://github.com/vcflib/vcflib/commits/master).
