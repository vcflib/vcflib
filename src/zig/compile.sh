#!/bin/sh

# zig build test
# creates library in zig-out
zig build $*
exit $?

# C example:
zig cc -o hello test_zig.c zig-out/lib/libzig.a
./hello
zig c++ -g -o hello test_zig.cpp zig-out/lib/libzig.a
