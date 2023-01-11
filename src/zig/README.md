# Zig code in vcflib

Some new vcflib functionality is written in the Zig programming language.
The first tool to use Zig code is `vcfcreatemulti`.
Zig is very fast and performance often beats C and C++ code.

## Build

To add zig to your build environment simply download the binary and add it to the path. Next setup the Cmake build environment in the usual way:

```
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build .
ctest .
```

The current [build workflow](../../.github/workflows/ci_test.yml) on github CI builds with zig. So you can check how it is done on Debian/Ubuntu.

## Coding against C++

Zig uses a C-ABI. To exchange data with C++ requires some data wrangling as is shown in [vcf-c-api](../vcf-c-api.cpp) and [vcf.zig](vcf.zig).
Nothing too serious, just remember to use the 0 sentinel when passing C strings.

### Error handling and warnings

Zig has a tight system for handling errors. For within zig errors are handled with `catch unreachable` so the system will simply croak.

For warnings we build up a list of them. Before exiting the code (from C++) be sure to call `zig_display_warnings()`.
