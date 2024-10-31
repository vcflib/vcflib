const zig_version = @import("builtin").zig_version;
const std = @import("std");

const test_targets = [_]std.Target.Query{
    .{}, // native
};

pub fn build(b: *std.Build) void {
    // Standard release options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall.
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    // const mode = b.standardReleaseOptions();

    const lib = b.addStaticLibrary(.{
        .optimize = optimize,
        .name = "zig",
        .target = target,
        .root_source_file = b.path("vcf.zig"),
        .pic = true,
    });
    // lib.setBuildMode(optimize);
    // lib.addObjectFile("../../build/libvcflib.a"); circular dependency
    switch (optimize) {
        .Debug, .ReleaseSafe => lib.bundle_compiler_rt = true,
        .ReleaseFast, .ReleaseSmall => {},
    }
    // lib.force_pic = true;
    // lib.emit_h = true; future version of zig?
    b.installArtifact(lib);
    // lib.install();

    const test_step = b.step("test", "Run unit tests");

    for (test_targets) |ttarget| {
        const unit_tests = b.addTest(.{
            .root_source_file = b.path("vcf.zig"),
            .target = b.resolveTargetQuery(ttarget),
        });

        const run_unit_tests = b.addRunArtifact(unit_tests);
        test_step.dependOn(&run_unit_tests.step);
    }
    
}
