// Zig variant handler

const std = @import("std");
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;
const p = @import("std").debug.print;

pub fn VarList() ArrayList(u8) {
    return struct {
        list: ArrayList(u8),

        const Self = @This();

        pub fn init(allocator: Allocator) Self {
            return Self{ .stack = ArrayList(u8).init(allocator) };
        }
    }
}

// void *zig_variant_window()

export fn zig_variant_window() *void {
    // p("And yes, we are back in zig: {s}\n\n",.{s});
    p("And yes, we are back in zig: {s}\n\n",.{"YES"});
}
