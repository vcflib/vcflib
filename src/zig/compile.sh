#!/bin/sh

# zig build test
# creates library in zig-out
echo "ZIG PATH $PATH"
zig version
export ZIG_GLOBAL_CACHE_DIR=$PWD # see https://github.com/ziglang/zig/issues/19400
zig build $*
exit $?

# C example:
zig cc -o hello test_zig.c zig-out/lib/libzig.a
./hello
zig c++ -g -o hello test_zig.cpp zig-out/lib/libzig.a
