#!/bin/sh

# zig build test
zig build # creates library in zig-out
exit $?

# C example:
zig cc -o hello test_zig.c zig-out/lib/libzig.a
./hello
zig c++ -g -o hello test_zig.cpp zig-out/lib/libzig.a
