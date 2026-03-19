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

    const lib = b.addLibrary(.{
        .name = "zig",
        .linkage = .static,
        .root_module = b.createModule(.{
            .optimize = optimize,
            .target = target,
            .root_source_file = b.path("vcf.zig"),
            .pic = true,
        }),
    });
    // lib.addObjectFile("../../build/libvcflib.a"); circular dependency
    b.installArtifact(lib);

    const test_step = b.step("test", "Run unit tests");

    for (test_targets) |ttarget| {
        const unit_tests = b.addTest(.{
            .root_module = b.createModule(.{
                .root_source_file = b.path("vcf.zig"),
                .target = b.resolveTargetQuery(ttarget),
            }),
        });

        const run_unit_tests = b.addRunArtifact(unit_tests);
        test_step.dependOn(&run_unit_tests.step);
    }

}
