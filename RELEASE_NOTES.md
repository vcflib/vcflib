## ChangeLog v1.0.2 (2020????)

This is a maintenance release of Freebayes:

+ Added support for ARM64 and PowerPC, see #292 (thanks @genisysram and @mr-c)
+ Added github actions for the issue tracker
+ Added tests to the cmake build system and build instructions to README
+ Removed tabixpp and htslib source dependencies, i.e., we are now using
  the distro provided libraries and include files through pkg-config.
  The `git submodules` are left in, just in case someone wants them.
