/// This zig module handles VCF samples and genotypes

const std = @import("std");
const mem = std.mem;
const fmt = std.fmt;
const ArrayList = std.ArrayList;
const expect = std.testing.expect;
const expectEqual = std.testing.expectEqual;
const expectEqualSlices = std.testing.expectEqualSlices;
const p = std.debug.print;
const test_allocator = std.testing.allocator;


const GENOTYPE_MISSING = -256;

/// Genotypes come as a list of integers separated by | (phased) or /
/// (unphased).  You can't really mix phased/unphased so we can track
/// that as a boolean.
const Genotypes = struct {
    genos: ArrayList(i64),
    phased: bool,

    const Self = @This();

    /// Convert a genotype (sample) string to a list of numbers
    fn to_num(str: []const u8) !ArrayList(i64) {
        var list = ArrayList(i64).init(test_allocator);
        var splits = std.mem.split(u8, str, "|");
        while (splits.next()) |chunk| {
            const i: i64 = 
                if (chunk[0] == '.') 
                GENOTYPE_MISSING
                else
                try fmt.parseInt(i64,chunk,10);
            list.append(i) catch unreachable;
        }
        // p("{s}",.{list});
        return list;
    }

    /// Take a 0-n indexed genotype and add offset idx. When the
    /// genotype is 0 (ref) or missing it is not changed.
    fn renumber(idx: usize, list: ArrayList(i64)) ArrayList(i64) {
        _ = idx;
        for (list.items) | g,i | {
            list.items[i] = 
                switch (g) {
                    0 => 0,
                    GENOTYPE_MISSING => GENOTYPE_MISSING,
                    else => g+@intCast(i64,idx)
            };
        }
        return list;
    }
    
    fn to_s(self: *const Self) !ArrayList([] const u8) {
        _ = self;
        var s = ArrayList([] const u8).init(test_allocator);
        for (self.genos.items) |g| {
                // parseInt to go to str
                // charDigit to int
                const result = try fmt.allocPrint(test_allocator, "{}", .{g});
                p("Result is {s}!\n", .{result});
                try s.append(result);
            }
        return s;
    }
};

fn split_genotypes(str: []const u8) *ArrayList([] const u8) {
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


/// Walks all genotypes in the list of variants and reduces them to
/// the genotypes of the combined variant. For examples 0|1 genotypes
/// get translated to their list numbers 0|6. This works for
/// heterozygous only (at this point).
pub fn reduce_renumber_genotypes(comptime T: type, vs: ArrayList(T)) !ArrayList([] const u8) {
    var ngenos = ArrayList([] const u8).init(test_allocator);
    _ = vs;
    // for (vs.items) |v,i| {
        // Fetch the genotypes from each variant
    //     for (v.genotypes().items) | geno | {
    //        // const geno2 = Genotypes.renumber_genotypes(i,geno);
    //        p("({s})",.{geno2});
    //        ngenos.append(geno2) catch unreachable;
    //    }
    //}
    return ngenos;
}

test "genotypes" {
    var list = ArrayList(i64).init(test_allocator);
    defer list.deinit();
    try list.append(0);
    try list.append(1);
    var gs = Genotypes{.genos = list, .phased = true};
    var genos = try gs.to_s();
    defer {
        for (genos.items) |item| {
                test_allocator.free(item);
            }
        genos.deinit();
    }
    p("YES {s}",.{genos.items});

    const gs2 = try Genotypes.to_num("1|0");
    defer gs2.deinit();
    try expect(gs2.items[0]==1);
    try expectEqual(gs2.items[1],0);
    const gs3 = try Genotypes.to_num("1|.");
    defer gs3.deinit();
    try expectEqual(gs3.items.len,2);
    try expectEqual(gs3.items[1],GENOTYPE_MISSING);
    const add3 = Genotypes.renumber(1,gs3);
    try expectEqual(add3.items[0],2);
    try expectEqual(add3.items[1],GENOTYPE_MISSING);
    const gs4 = try Genotypes.to_num(".|2");
    defer gs4.deinit();
    try expectEqualSlices(i64, gs4.items, &.{ GENOTYPE_MISSING, 2 });
    const add4 = Genotypes.renumber(1,gs4);
    try expectEqualSlices(i64, add4.items, &.{ GENOTYPE_MISSING, 3 });
}

test "split genotypes" {
    const input_genotypes = "1|0 .|1 0|1 1|1";
    try std.testing.expectEqual(split_genotypes(input_genotypes).items.len, 4);
}
