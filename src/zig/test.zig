/// Test module

const var_id = @import("vcf.zig").var_id;

const hello = "Hello World from Zig";

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

export fn zig_variant_window() *void {
    // p("And yes, we are back in zig: {s}\n\n",.{s});
    p("And yes, we are back in zig: {s}\n\n",.{"YES"});
}

test "hello zig" {
    try expectEqual(hello_zig(hello),hello);
}
