/*
 * Simple progress bar implementation in C. 
 *
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

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

/* Helper function to get the current time in usecs past the epoch. */
static uint64_t get_timestamp(void) {
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
        printf("%lluh %llum %llus    ", hours, minutes, seconds);
    }
    else if (minutes) {
        printf("%llum %02llus        ", minutes, seconds);
    }
    else {
        printf("%llu.%llus           ", seconds, mseconds);
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

    printf("   Progress: %6.2f%% \t%s", percentage, BAR_START);

    for (i = 0; i < num_blocks; i++) {
        printf("%s", PROGRESS_BLOCK);
    }

    if (num_subblocks) {
        printf("%s", subprogress_blocks[num_subblocks]);
        i++;
    }

    for (; i < PROGRESS_BAR_WIDTH; i++) {
        printf(" ");
    }

    
    printf("%s\t", BAR_STOP);

    if (percentage < 100.0) {
        printf("ETA: ");
        print_timedelta(remaining);
    }
    else {
        printf("                          ");
    }
    printf("\r");
    fflush(stdout);
}

int main(int argc, char **argv) {
    int i;
    double amount = 0.0;
    uint64_t start = get_timestamp();

    for (i = 0; i < 10000; i++) {
        amount += 0.01;
        print_progress(amount, start);
        usleep(2000000 / (i + 1));
    }

    printf("\n");
}
