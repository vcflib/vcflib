/*
 * Simple progress bar implementation in C.
 *
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/time.h>
#include <progress.h>

/* Specify how wide the progress bar should be. */
#define PROGRESS_BAR_WIDTH 50

/*
 * Alternative progress bar where each block grows
 * vertically instead of horizontally.
 */
/* #define VERTICAL */

/* Various unicode character definitions. */
#define BAR_START "\u2595"
#define BAR_STOP  "\u258F"
#define PROGRESS_BLOCK     "\u2588"

#ifdef VERTICAL
static const char * subprogress_blocks[] = { " ",
                                             "\u2581",
                                             "\u2582",
                                             "\u2583",
                                             "\u2584",
                                             "\u2585",
                                             "\u2586",
                                             "\u2587"
};
#else
static const char * subprogress_blocks[] = { " ",
                                             "\u258F",
                                             "\u258E",
                                             "\u258D",
                                             "\u258C",
                                             "\u258B",
                                             "\u258A",
                                             "\u2589"
};
#endif

#define NUM_SUBBLOCKS (sizeof(subprogress_blocks) / sizeof(subprogress_blocks[0]))

/* Get file size */
off_t get_file_size(const char *filename)
{
    struct stat stat_buf;
    off_t rc = stat(filename, &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

/* Helper function to get the current time in usecs past the epoch. */
uint64_t get_timestamp(void) {
    struct timeval tv;
    uint64_t stamp = 0;
    gettimeofday(&tv, NULL);
    stamp = tv.tv_sec * 1000000 + tv.tv_usec;
    return stamp;
}

/* Helper function to print a usecs value as a duration. */
static void print_timedelta(uint64_t delta) {

    uint64_t delta_secs = delta / 1000000;
    uint64_t hours    = delta_secs / 3600;
    uint64_t minutes  = (delta_secs - hours * 3600) / 60;
    uint64_t seconds  = (delta_secs - hours * 3600 - minutes * 60);
    uint64_t mseconds = (delta / 100000) % 10;

    if (hours) {
        fprintf(stderr,"%lluh %llum %llus    ", hours, minutes, seconds);
    }
    else if (minutes) {
        fprintf(stderr,"%llum %02llus        ", minutes, seconds);
    }
    else {
        fprintf(stderr,"%llu.%llus           ", seconds, mseconds);
    }
}

/*
 * Main interface function for updating the progress bar. This
 * function doesn't print a newline, so you can call it iteratively
 * and have the progress bar grow across the screen. The caller can
 * print a newline when the're ready to go to a new line.
 *
 * percentage: a double between 0.0 and 100.0 indicating the progress.

 * start: usecs timestamp for when the task started, for calculating
 *        remaining time.
 */
void print_progress(double percentage, uint64_t start) {
    size_t i;
    size_t total_blocks = PROGRESS_BAR_WIDTH * NUM_SUBBLOCKS;
    size_t done = round(percentage / 100.0 * total_blocks);
    size_t num_blocks = done / NUM_SUBBLOCKS;
    size_t num_subblocks = done % NUM_SUBBLOCKS;

    uint64_t now = get_timestamp();
    uint64_t elapsed = now - start;
    uint64_t estimated_total = elapsed / (percentage / 100.0);
    uint64_t remaining = estimated_total - elapsed;

    fprintf(stderr,"   Progress: %6.2f%% \t%s", percentage, BAR_START);

    for (i = 0; i < num_blocks; i++) {
        fprintf(stderr,"%s", PROGRESS_BLOCK);
    }

    if (num_subblocks) {
        fprintf(stderr,"%s", subprogress_blocks[num_subblocks]);
        i++;
    }

    for (; i < PROGRESS_BAR_WIDTH; i++) {
        fprintf(stderr," ");
    }


    fprintf(stderr,"%s\t", BAR_STOP);

    if (percentage < 100.0) {
        fprintf(stderr,"ETA: ");
        print_timedelta(remaining);
    }
    else {
        fprintf(stderr,"                          ");
    }
    fprintf(stderr,"\r");
    fflush(stdout);
}

#ifdef WITH_MAIN

int main(int argc, char **argv) {
    int i;
    double amount = 0.0;
    uint64_t start = get_timestamp();

    for (i = 0; i < 10000; i++) {
        amount += 0.01;
        print_progress(amount, start);
        usleep(2000000 / (i + 1));
    }

    fprintf(stderr,"\n");
}

#endif
