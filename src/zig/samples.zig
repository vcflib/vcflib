/// This zig module handles VCF samples and genotypes

const std = @import("std");
const ArrayList = std.ArrayList;
const expect = @import("std").testing.expect;
const expectEqual = @import("std").testing.expectEqual;
const p = @import("std").debug.print;
const var_id = @import("vcf.zig").var_id;

const hello = "Hello World from Zig";

// extern fn call_c(*anyopaque) void;
// extern fn get_name([*] const u8) [*] const u8;

export fn hello_zig(msg: [*] const u8) [*]const u8 {
    const result = msg;
    return result;
}

export fn zig_process_vector(vsize: u64, v: [*][*] const u8) void {
    // vector<string> genotypes{ "1/0", "2/.", "3/1" };
    const s = v[1][0];
    const s1 = v[1][1];
    const s2 = v[1][2];
    p("HELLO! {any} {any} {any} {c}:{c}:{c}\n",.{vsize, &v[0], &v[1], s,s1,s2}); // expect 1 and dot
    // std.testing.expectEqual(&s, "1/0") catch unreachable;
}

export fn zig_process_opaque_ptr(variant: *anyopaque) void {
    // call_c(variant);
    const c_str = var_id(variant);
    const s = @ptrCast([*c]const u8, c_str);
    p("And yes, we are back in zig: {s}\n\n",.{s});
}

fn split_genotypes(str: []const u8) *ArrayList([] const u8) {
    const test_allocator = std.testing.allocator;

    var list = ArrayList([] const u8).init(test_allocator);
    defer list.deinit();
    var splits = std.mem.split(u8, str, " ");
    while (splits.next()) |chunk| {
            list.append(chunk) catch |err| {
                std.debug.print("out of memory {e}\n", .{err});
            };
        }
    return &list;
}

/// Take a 0-1 indexed genotype and convert it to an indexed number
pub fn renumber_genotypes(idx: usize, str: [] const u8) [] const u8 {
    _ = idx;
    return str;
}

test "hello zig" {
    try expectEqual(hello_zig(hello),hello);
}

test "split genotypes" {
    const input_genotypes = "1|0 .|1 0|1 1|1";
    try std.testing.expectEqual(split_genotypes(input_genotypes).items.len, 4);
}
