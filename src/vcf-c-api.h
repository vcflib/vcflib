/*
  C API for vcflib
 */

// do not use extern "C" because this file ought to be ready for a C compiler (!C++)

// VCF constructors
#pragma once

void *var_parse(const char *line, bool parse_samples);

// VCF variant accessors
const char *var_id(void *var);
const long var_pos(void *var);
const char *var_ref(void *var);
const unsigned long var_alt_num(void *var);
const char **var_alt(void *var, const char ** ret);
const unsigned long var_info_num(void *variant, const char *name);
const char **var_info(void *var, const char *name, const char **ret);
const char **var_geno(void *var, const char **ret);

void var_set_id(void *var, const char *);
void var_set_ref(void *var, const char *);

void var_clear_alt(void *var);
void var_set_alt(void *var, const char *alt, long idx);

void var_clear_info(void *var, const char *name);
void var_set_info(void *var, const char *name, const char *values, long idx);

void var_clear_sample(void *var, long idx);
void var_set_sample(void *var, const char *name, long idx);

const unsigned long var_samples_num(void *var);

// Zig functionality

void *zig_create_multi_allelic(void *retvar, void *varlist[], long size);

// Some test functions
void testme();
void *zig_variant_window();
void zig_variant_window_cleanup(void *varwin);
void win_push(void *varwin, void *var);
long win_size(void *varwin);


char *hello_zig2(char *s);

void *zig_variant(void *var);

void zig_display_warnings();

void zig_cleanup();
