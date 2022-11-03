// Zig variant handler

const std = @import("std");
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;
const p = std.debug.print;

pub fn VarList() ArrayList(u8) {
    return struct {
        list: ArrayList(u8),

        const Self = @This();

        pub fn init(allocator: Allocator) Self {
            return Self{ .stack = ArrayList(u8).init(allocator) };
        }
    }
}

