#ifndef KERNELS_H
#define KERNELS_H

#include "matrix.h"
#include <stddef.h>

typedef struct {
    int nt;
    int tile;
} KCfg;

int mv_mt(const Mat *A, const Vec *x, Vec *y, KCfg cfg);

int mm_mt(const Mat *A, const Mat *B, Mat *C, KCfg cfg);

int dt_mt(const Vec *x, const Vec *y, double *out, int nt);

int ax_mt(double a, const Vec *x, Vec *y, int nt);

#endif
