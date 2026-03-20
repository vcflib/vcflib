# TODO

- [ ] vcfcreatemulti: fix problem with slow and wrong complex regions (implement backtrack)
      + [ ] check for indels which are really the same
      + [ ] combine vcfwave duplicated functionality
- bgzip
- tabix -p vcf my_file.vcf.gz
- pangenie, vg deconstruct, vcfbub

# Debian

  Debian version: 1.0.12 (current vcflib HEAD is 1.0.15)

  Key findings from debian/rules:
  - -DZIG=OFF -- Zig is completely disabled
  - -DWFA_GITMODULE=OFF -- uses system wfa2-lib
  - -DHTSLIB_FOUND=ON -DTABIX_FOUND=ON -- uses system htslib/tabixpp
  - -DOPENMP=ON
  - Tests are deactivated in the patch series (#tests ## FIXME: deactivated for the moment)

  Patches (9 total, mostly about using Debian-packaged dependencies instead of bundled ones):
  1. no_fsom -- removes unused fsom dependency
  2. use_debian_packaged_smithwaterman.patch -- system libsmithwaterman
  3. use_debian_packaged_fastahack.patch -- system libfastahack
  4. use_debian_packaged_libssw.patch -- system libssw
  5. shared_lib.patch -- builds shared lib instead of static
  6. pkg-config.patch -- adds a .pc file for freebayes
  7. tabix_linking -- link against system tabixpp
  8. wfa2_linking -- workaround for wfa2 library detection
  9. cmakelists-fixups.patch -- minor CMakeLists fixes

  So yes -- zig is simply turned off (-DZIG=OFF). The patches are all about unbundling dependencies to use Debian
  system packages, plus building a shared library. No zig-related patches at all. They're also 3 minor versions behind
  and have their tests disabled.
