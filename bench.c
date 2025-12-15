#include "bench.h"
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>

double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}
