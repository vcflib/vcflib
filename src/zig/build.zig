const std = @import("std");

pub fn build(b: *std.build.Builder) void {
    // Standard release options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall.
    const mode = b.standardReleaseOptions();

    const lib = b.addStaticLibrary("zig", "vcf.zig");
    lib.setBuildMode(mode);
    switch (mode) {
        .Debug, .ReleaseSafe => lib.bundle_compiler_rt = true,
        .ReleaseFast, .ReleaseSmall => lib.disable_stack_probing = true,
    }
    lib.force_pic = true;
    // lib.emit_h = true; future version of zig?
    lib.install();

    const main_tests = b.addTest("vcf.zig");
    main_tests.setBuildMode(mode);
    main_tests.addLibPath("../../build");
    main_tests.addObjectFile("../../build/libvcflib.a");

    const test_step = b.step("test", "Run library tests");
    test_step.dependOn(&main_tests.step);
}
