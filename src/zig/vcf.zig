//! This zig module provides functionality for cleaning genotypes -
//! mostly for combining multiple records into one. This code is partly
//! in C++ in vcfwave and the older vcfallelicprimitives as well as C++
//! version of vcfcreatemulti. None of these implementations did a great
//! job and I felt a zig version would be faster and cleaner.

const std = @import("std");
// const variant = @import("variant");
const mem = @import("std").mem;
const fmt = std.fmt;
const samples = @import("./samples.zig");
const expectEqual = @import("std").testing.expectEqual;
const expect = @import("std").testing.expect;
const ArrayList = std.ArrayList;
const StringList = ArrayList([] const u8);
const p = @import("std").debug.print;

const VCFError = error{
    UnexpectedOrder,
    MultiAltNotSupported
};

const hello = "Hello World from Zig";

// C++ constructor
extern fn var_parse(line: [*c] const u8, parse_samples: bool) *anyopaque;

// C++ accessors for Variant object
pub extern fn var_id(* anyopaque) [*c] const u8;
extern fn var_pos(* anyopaque) u64;
extern fn var_ref(* anyopaque) [*c] const u8;
extern fn var_alt_num(variant: *anyopaque) usize;
extern fn var_samples_num(variant: *anyopaque) usize;
extern fn var_info_num(variant: *anyopaque, name: [*c] const u8) usize;
extern fn var_clear_alt(variant: *anyopaque) void;
extern fn var_alt(variant: * anyopaque, buf: [*c]* anyopaque) [*c][*c] const u8;
extern fn var_info(variant: * anyopaque, name: [*c] const u8, buf: [*c]* anyopaque) [*c][*c] const u8;
extern fn var_geno(variant: * anyopaque, buf: [*c]* anyopaque) [*c][*c] const u8;
extern fn var_clear_info(variant: *anyopaque, name: [*c] const u8) void;
extern fn var_set_id(?* anyopaque, [*c] const u8) void;
extern fn var_set_ref(?* anyopaque, [*c] const u8) void;
extern fn var_set_alt(?* anyopaque, [*c] const u8, usize) void;
extern fn var_set_info(?* anyopaque, name: [*c] const u8, value: [*c] const u8, int: usize) void;
extern fn var_clear_sample(variant: *anyopaque, usize) void;
extern fn var_set_sample(?* anyopaque, [*c] const u8, usize) void;
// extern fn var_set_geno(?* anyopaque, value: [*c] const u8, int: usize) void;
extern fn call_c([*] const u8) void;

export fn hello_zig2(msg: [*] const u8) [*]const u8 {
    const result = msg;
    return result;
}

// const allocator = std.testing.allocator;
var gpa = std.heap.GeneralPurposeAllocator(.{}){};
const allocator = gpa.allocator();

var warnings = std.StringHashMap(bool).init(allocator); // Note: not thread safe

pub fn warning(str: [] const u8) !void {
    try warnings.put(str, true);
}

/// Called by C++
export fn zig_display_warnings() void {
    // warnings.append("NO WARNINGS+") catch unreachable;
    var iterator = warnings.iterator();

    while (iterator.next()) |msg| {
        p("WARNING: {s}\n",.{msg.key_ptr.*});
    }
}

/// Variant struct maps over the equivalent C++ version - copying data
/// forth to zig and back to C++.
const Variant = struct {
    v: *anyopaque,

    const Self = @This();

    inline fn to_slice0(c_str: [*c] const u8) [:0] const u8 {
        return std.mem.span(@as([*:0]const u8, @ptrCast(c_str)));
    }

    inline fn to_slice(c_str: [*c] const u8) [] const u8 {
        return std.mem.span(c_str);
    }

    inline fn to_cstr(str: [:0] const u8) [*c]const u8 {
        return @as([*c]const u8,@ptrCast(str));
    }

    inline fn to_cstr0(str: [] const u8) [*c]const u8 { // not sure this works because we need final zero
        //var s0 = str.toOwnedSliceSentinel(0);
        return @as([*c]const u8,@ptrCast(str));
    }

    pub fn id(self: *const Self) [:0]const u8 {
        // const buf: [*c]const u8 = var_id(self.v);
        // const str = std.mem.span(@ptrCast([*:0]const u8, buf));
        // return str;
        return to_slice0(var_id(self.v));
    }

    /// Get the C++ pos
    pub fn pos(self: *const Self) u64 {
        return var_pos(self.v);
    }

    /// Get the C++ ref
    pub fn ref(self: *const Self) [] const u8 {
        // const buf: [*c]const u8 = var_ref(self.v);
        // const str = std.mem.span(@ptrCast([*:0]const u8, buf));
        // return str;
        return to_slice(var_ref(self.v));
    }

    /// Get the C++ alts
    pub fn alt(self: *const Self) ArrayList([] const u8) {
        var list: ArrayList([] const u8) = .{};
        const altsize = var_alt_num(self.v);
        const buffer = allocator.alloc(*anyopaque, altsize) catch unreachable;
        defer allocator.free(buffer);
        const res = var_alt(self.v,@as([*c]* anyopaque,@ptrCast(buffer)));
        var i: usize = 0;
        while (i < altsize) : (i += 1) {
            const s = res[i];
            const s2 = to_slice(s);
            list.append(allocator, s2) catch unreachable;
        }
        return list;
    }

    /// Get the C++ infos
    pub fn info(self: *const Self, name: [] const u8) ArrayList([] const u8) {
        const c_name = to_cstr0(name);
        var list: ArrayList([] const u8) = .{};
        const size = var_info_num(self.v,c_name);
        const buffer = allocator.alloc(*anyopaque, size) catch unreachable;
        defer allocator.free(buffer);
        const res = var_info(self.v,c_name,@as([*c]* anyopaque,@ptrCast(buffer)));
        var i: usize = 0;
        while (i < size) : (i += 1) {
            // list.append(buffer[i]) catch unreachable;
            const s = res[i];
            // const s1 = to_slice(s);
            // p("<{d}:{d}><{any}--{s}--{any}>\n",.{i,altsize,s,s1,buffer[i]});
            const s2 = to_slice(s);
            // p("{s}\n",.{s2});
            list.append(allocator, s2) catch unreachable;
        }

        return list;
    }

    /// Get the C++ genotypes as a list
    pub fn genotypes(self: *const Self) ArrayList([] const u8) {
        // p("Inside genotypes:\n",.{});
        var list: ArrayList([] const u8) = .{};
        const size = var_samples_num(self.v);
        const buffer = allocator.alloc(*anyopaque, size) catch unreachable;
        defer allocator.free(buffer);
        const res = var_geno(self.v,@as([*c]* anyopaque,@ptrCast(buffer)));
        var i: usize = 0;
        while (i < size) : (i += 1) {
            const s = res[i];
            const s2 = to_slice(s);
            // p("<{s}>",.{s2});
            list.append(allocator, s2) catch unreachable;
        }
        return list;
    }

    /// Set C++ ref
    pub fn set_ref(self: *const Self, nref: [:0] const u8) void {
        var_set_ref(self.v,@as([*c]const u8,@ptrCast(nref)));
    }

    /// Set C++ alts
    pub fn set_alt(self: *const Self, nalt: ArrayList([*:0] const u8)) void {
        // Create ptrlist
        var_clear_alt(self.v);
        var i: usize = 0;
        // p("<{s}>", .{nalt.items});
        while (i < nalt.items.len) : (i += 1) {
            // p("<{s}>\n",.{nalt.items[i]});
            // var x = to_cstr(nalt.items[i]);
            // var x: [:0] const u8 =
            //     nalt.items[i];

            // var_set_alt(self.v,@ptrCast([*c] const u8,x),i);
            // var_set_alt(self.v,x,i);
            var_set_alt(self.v,nalt.items[i],i);
        }
    }

    /// Set C++ infos
    pub fn set_info(self: *const Self, name: [] const u8, data: ArrayList([] const u8)) void {
        const c_name = to_cstr0(name);
        var_clear_info(self.v,c_name);
        var i: usize = 0;
        while (i < data.items.len) : (i += 1) {
            var_set_info(self.v,c_name,to_cstr0(data.items[i]),i);
        }
    }

    pub fn set_samples(self: *const Self, nsamples: ArrayList([] const u8)) void {
        var i: usize = 0;
        while (i < nsamples.items.len) : (i += 1) {
            var_clear_sample(self.v,i);
            const s = nsamples.items[i];
            var buffer = allocator.alloc(u8, s.len + 1) catch unreachable;
            defer allocator.free(buffer);

            for (s, 0..)  | c,j | {
                buffer[j] = c;
            }
            buffer[s.len] = 0;
            var_set_sample(self.v,to_cstr0(buffer),i);
        }
    }
};

// by @Cimport:
// pub extern fn zig_create_multi_allelic(retvar: ?*anyopaque, varlist: [*c]?*anyopaque, size: c_long) ?*anyopaque;

// Obsolete test version of multi_allelic
export fn zig_create_multi_allelic2(variant: ?*anyopaque, varlist: [*c]?* anyopaque, size: usize) ?*anyopaque {
    const v1 = var_parse("TEST\t1\t2\t3\t4\tt5\t6",false);
    _ = v1;
    const c_var = var_parse("a\t281\t>1>9\tAGCCGGGGCAGAAAGTTCTTCCTTGAATGTGGTCATCTGCATTTCAGCTCAGGAATCCTGCAAAAGACAG\tCTGTCTTTTGCAGGATTCCTGTGCTGAAATGCAGATGACCGCATTCAAGGAAGAACTATCTGCCCCGGCT\t60.0\t.\tAC=1;AF=1;AN=1;AT=>1>2>3>4>5>6>7>8>9,>1<8>10<6>11<4>12<2>9;NS=1;LV=0\tGT\t1",false);
    var v2 = Variant{.v = c_var};
    p("---->{s}\n",.{v2.id().ptr});
    expect(mem.eql(u8, v2.id(), ">1>9")) catch unreachable;
    expectEqual(v2.id().len,">1>9".len) catch |err| {
        std.debug.print("{} <-> {s}\n", .{err,v2.id()});
    };
    const c_str = var_id(variant.?);
    // const s = @ptrCast([*c]const u8, c_str);
    const s = @as([*c]const u8,@ptrCast(c_str));
    p("And yes, we are back in zig: {s} -- {}\n\n",.{s,size});

    const p3 = @as(* anyopaque, @ptrCast(varlist[3]));
    const s3 = var_id(p3);
    var v = Variant{.v = varlist[3].?};
    p("id={s} !{s}! pos={d} ref={s}\n",.{s3,v.id(),v.pos(),v.ref()});


    const as_slice: [:0]const u8 = std.mem.span(s3); // makes 0 terminated slice (sentinel value is zero byte)
    std.testing.expectEqualStrings(as_slice, ">3655>3662_4") catch |err| {
        std.debug.print("{} {s}\n", .{err,as_slice});
    };

    // expectEqual(variant,@intToPtr(*anyopaque,varlist[0])) catch unreachable;
    // var vars = @ptrCast([*] u8, varlist);
    // Now walk the list
    var i:u64 = 0;
    for (varlist[0..size]) |ptr| {
             i = i + 1;
             const p2 = @as(* anyopaque, @ptrCast(ptr));
             const s2 = var_id(p2);
             p("num = {}",.{i});
             p("id = {s}, pos = {d}\n",.{s2,var_pos(p2)});
         }
    return variant;
}

/// This function is the main entry point and called from C++ to
/// reduce a set of variants to a single VCF record and adjusting
/// genotypes accordingly. Essentially a list of variants is passed
/// that overlap. This code simplifies ref and alts for each variant
/// and adjusts the metrics for AF, AC, sample genotype index etc.
///
export fn zig_create_multi_allelic(variant: ?*anyopaque, varlist: [*c]?* anyopaque, size: usize) *anyopaque {
    // Create vs as a list of variants
    // var hanging_pointer: ArrayList([] const u8) = .{};
    // hanging_pointer.append(allocator, "C") catch unreachable;

    var mvar = Variant{.v = variant.?}; // FIXME: we need to clean this small struct up from C++
    var vs: ArrayList(Variant) = .{};
    defer vs.deinit(allocator);

    var i: usize = 0;
    while (i < size) : (i += 1) { // use index to access *anyopaque
        const v = Variant{.v = varlist[i].?};
        vs.append(allocator, v) catch unreachable;
    }

    // Get the reference and update mvar (multi VCF record containing multiple variants)
    var nref = expand_ref(Variant,vs) catch unreachable;
    defer nref.deinit(allocator);
    const c_nref = nref.toOwnedSliceSentinel(allocator, 0) catch unreachable;
    mvar.set_ref(c_nref);

    // Get the alts and update mvar
    const first = vs.items[0];
    var nalt = expand_alt(Variant,first.pos(),c_nref,vs) catch unreachable;
    defer nalt.deinit(allocator);
    mvar.set_alt(nalt);

    // Get infos and update mvar
    const list = [_][] const u8{ "AN","AT","AC","AF","INV","TYPE" };
    for (list) |name| {
            var at = expand_info(Variant,name,vs) catch unreachable;
            defer at.deinit(allocator);
            mvar.set_info(name,at);
        }

    // Get genotypes and update mvar
    var genotypes_result = samples.reduce_renumber_genotypes(Variant,vs) catch unreachable;
    defer genotypes_result.s_samples.deinit(allocator);

    mvar.set_samples(genotypes_result.s_samples);
    var ninfo: ArrayList([] const u8) = .{};
    defer ninfo.deinit(allocator);

    if (genotypes_result.g_err != samples.VcfSampleError.None) {
        ninfo.append(allocator, "ALTPROBLEM") catch unreachable;
        mvar.set_info("MULTI",ninfo);
    }

    return mvar.v;
}

/// The C++ code should call this to cleanup

export fn zig_cleanup() void {
    warnings.deinit();

    // ---- Not cleaning up the GPA unless we are debugging
    // std.debug.assert(!gpa.deinit());
}


fn refs_maxpos(comptime T: type, list: ArrayList(T)) usize {
    var mpos = list.items[0].pos();
    for (list.items) |v| {
            const npos = v.pos() + v.ref().len;
            if (npos > mpos)
                mpos = npos;
        }
    return mpos;
}

//    // set maxpos to the most outward allele position + its reference size
//    auto maxpos = first.position + first.ref.size();
//    for (auto v: vars) {
//         if (maxpos < v.position + v.ref.size()) {
//             maxpos = v.position + v.ref.size();
//         }
//     }

const MockVariant = struct {
    id_: [] const u8 = "TEST",
    pos_: u64,
    ref_: [] const u8,
    alt_: ArrayList([] u8) = undefined,

    const Self = @This();

    pub fn id(self: *const Self) []const u8 {
        return self.id_;
    }

    pub fn pos(self: *const Self) usize {
        return self.pos_;
    }

    pub fn ref(self: *const Self) [] const u8 {
        return self.ref_;
    }

    pub fn alt(self: *const Self) ArrayList([] u8) {
        return self.alt_;
    }

};

/// Expands reference to overlap all variants
fn expand_ref(comptime T: type, list: ArrayList(T)) !ArrayList(u8) {
    var res: ArrayList(u8) = .{};
    // defer res.deinit(allocator);
    // try res.append(allocator, 'T');
    const first = list.items[0];
    // concat2(&res,first.ref);
    try res.appendSlice(allocator, first.ref());
    // p("!{s}!",.{res});
    // defer allocator.free(result);

    const left0 = first.pos();
    for (list.items) |v| {
        const right0 = left0 + res.items.len;
        const left1 = v.pos();
        const right1 = left1 + v.ref().len;
        //            ref    sdiff
        // ref0     |AAAAA|------->|
        // ref1      |AAAAAAAAAAAAA|
        //           |--->| append |
        //            pdiff

        if (right1 > right0) {

            const sdiff = right1 - right0; // diff between ref0 and ref1 right positions
            const pdiff = right0 - left1; // diff between ref0 right and ref1 left
            // newref = ref + append
            try res.appendSlice(allocator, v.ref()[pdiff..pdiff+sdiff]);
        }
    }
    return res;
}

fn expand_alt(comptime T: type, pos: usize, ref_seq: [] const u8, list: ArrayList(T)) !ArrayList([*:0] const u8) {
    // add alternates and splice them into the reference. It does not modify the ref.
    var nalt: ArrayList([*:0] const u8) = .{};

    for (list.items) |v| {
        const p5diff = v.pos() - pos; // always >= 0 - will raise error otherwise
        const before = ref_seq[0..p5diff]; // leading ref

        // ref0 has been expanded in a previous step to cover the full variant.
        // the original code only deals with p3diff > 0.
        //
        // SNP
        //            ref
        // ref0     |AAAAAAAA|
        //        p5diff    p3diff = +3
        // SNP           C--->
        //
        // Insertion:
        //            ref
        // ref0     |AAAAA|------->|
        //        p5diff    p3diff = -8 (start = 5--8 = 13
        // ref1      |AAAAAAAAAAAAA|
        //
        // Deletion:
        //            ref
        // ref0     |AAAAA|
        //        p5diff    p3diff = +2 (start = 5-2 = 3
        // ref1      |AA|--

        const right0 = pos + ref_seq.len;
        const right1 = v.pos() + v.ref().len;
        const p3diff:i64 = @as(i64, @intCast(right0)) - @as(i64, @intCast(right1));

        var after: [] const u8 = undefined;
        if (p3diff > 0 and p3diff < ref_seq.len) {
            const last  = ref_seq.len - @as(usize, @intCast(p3diff));
            after = ref_seq[last..];
        }
        else after = "";
        if (v.alt().items.len > 1) {
            warning("This code only supports one ALT allele per record: bailing out\nTry normalising the data with `bcftools norm -m-`") catch unreachable;
            // p("Error: this code only supports one ALT allele per record (WIP/FIXME)\n",.{});
            return error.MultiAltNotSupported;
        }

        for (v.alt().items) | a | {
            var n: ArrayList(u8) = .{};
            defer n.deinit(allocator);
            if (p3diff != 0 or p5diff != 0) {
                // p("{any}-{s},{s}\n",.{p3diff,before,after});
                try n.appendSlice(allocator, before);
                try n.appendSlice(allocator, a);
                try n.appendSlice(allocator, after);
                // try nalt.append(allocator, n.items);
                // n copied to nalt and emptied (no longer in care of n)
                try nalt.append(allocator, n.toOwnedSliceSentinel(allocator, 0) catch unreachable);
                // p("new alt={s}\n",.{new.items});
            } else {
                try n.appendSlice(allocator, a);
                // try n.toOwnedslice(alt);
                try nalt.append(allocator, n.toOwnedSliceSentinel(allocator, 0) catch unreachable);
                // try nalt.append(allocator, n.items);
            }
        }
    }
    return nalt; // caller needs to clean up
}

fn expand_info(comptime T: type, name: [] const u8, list: ArrayList(T)) !ArrayList([] const u8) {
    var ninfo: ArrayList([] const u8) = .{};
    for (list.items) |v| {
        for (v.info(name).items) | info_item | {
            // try ninfo.append(allocator, info);
            // p("{s}",.{info_item});
            ninfo.append(allocator, info_item) catch unreachable;
        }
    }
    return ninfo;
}

test "hello zig" {
    try expectEqual(hello_zig2(hello),hello);
}


test "variant ref expansion" {

    // var c_var = var_parse("a\t281\t>1>9\tAGCCGGGGCAGAAAGTTCTTCCTTGAATGTGGTCATCTGCATTTCAGCTCAGGAATCCTGCAAAAGACAG\tCTGTCTTTTGCAGGATTCCTGTGCTGAAATGCAGATGACCGCATTCAAGGAAGAACTATCTGCCCCGGCT\t60.0\t.\tAC=1;AF=1;AN=1;AT=>1>2>3>4>5>6>7>8>9,>1<8>10<6>11<4>12<2>9;NS=1;LV=0\tGT\t1",false);
    // var v2 = Variant{.v = c_var};
    // p("---->{s}\n",.{v2.id()});
    // expect(mem.eql(u8, v2.id(), ">1>9")) catch |err| {
    //     std.debug.print("{e} <-> {s}\n", .{err,v2.id()});
    // };

}

test "mock variant" {
    const talloc = std.testing.allocator;
    var list: ArrayList(MockVariant) = .{};
    defer list.deinit(talloc);

    const v1 = MockVariant{ .pos_ = 10, .ref_ = "AAAA" };
    try expect(std.mem.eql(u8, v1.id(), "TEST"));
    try list.append(talloc, v1);
    const v2 = MockVariant{ .pos_ = 10, .ref_ = "AAAAA" };
    try list.append(talloc, v2);
    const v3 = MockVariant{ .pos_ = 10, .ref_ = "AAAAACC" };
    try list.append(talloc, v3);
    const maxpos = refs_maxpos(MockVariant,list);
    p("<{any}>",.{maxpos});
    try expect(maxpos == 17);

    const nref = try expand_ref(MockVariant,list);
    // defer std.testing.allocator.free(nref);
    // p("<{s}>",.{nref});
    // p("!{s}!",.{nref});
    try expect(nref.items.len == 7);
    try expect(std.mem.eql(u8, nref.items, "AAAAACC"));
    nref.deinit(allocator);
}

test "variant alt expansion" {
    // var hanging_pointer: ArrayList([] const u8) = .{};
    // hanging_pointer.append(allocator, "C") catch unreachable;

    const talloc = std.testing.allocator;
    var list: ArrayList(MockVariant) = .{};
    defer {
        list.deinit(talloc);
    }

    var alt1: std.ArrayList([] u8) = .{};
    defer alt1.deinit(talloc);
    var a1 = [_]u8{'c', 'c'};
    try alt1.append(talloc, a1[0..]);
    const v1 = MockVariant{ .pos_ = 10, .ref_ = "AAAA", .alt_ = alt1 };
    try expect(std.mem.eql(u8, v1.id(), "TEST"));
    try list.append(talloc, v1);

    var alt2: std.ArrayList([] u8) = .{};
    defer alt2.deinit(talloc);
    var a2 = [_]u8{'c'};
    try alt2.append(talloc, a2[0..]);
    const v2 = MockVariant{ .pos_ = 10, .ref_ = "AAAAA", .alt_ = alt2 };
    try list.append(talloc, v2);

    var alt3: std.ArrayList([] u8) = .{};
    defer alt3.deinit(talloc);
    var a3 = [_]u8{'c', 'c', 'c', 'c'};
    try alt3.append(talloc, a3[0..]);
    const v3 = MockVariant{ .pos_ = 10, .ref_ = "CC", .alt_ = alt3 };
    try list.append(talloc, v3);
    // const nalt = try expand_alt(MockVariant,10,"AAAAACC",list);
    // defer {
        // for (nalt.items) |item| {
        //         std.testing.allocator.free(item);
        //     }
        // std.testing.allocator.free(nalt);
    //    nalt.deinit(allocator);
    // }
    // expect(nalt.items.len == 3) catch |e| {
    //     p("{e}: {d}",.{e,nalt.items.len});
    //    return;
    //};
}

test {
    _ = @import("samples.zig");
}
