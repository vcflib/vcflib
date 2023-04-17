/*
 * Simple progress bar implementation in C.
 *
 *
 */

#pragma once

#include <sys/stat.h>

#ifdef __cplusplus
extern "C" {
#endif

    extern off_t get_file_size(const char *filename);
    extern void print_progress(double percentage, uint64_t start);
    extern uint64_t get_timestamp(void);

#ifdef __cplusplus
}
#endif
