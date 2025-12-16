#include "kernels.h"
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

extern volatile sig_atomic_t g_stop;

static void row_range(size_t nrows, int tid, int nt, size_t *i0, size_t *i1) {
    size_t base = nrows / (size_t)nt;
    size_t rem  = nrows % (size_t)nt;
    size_t start = (size_t)tid * base + (size_t)((tid < (int)rem) ? tid : rem);
    size_t extra = (tid < (int)rem) ? 1u : 0u;
    size_t end = start + base + extra;
    *i0 = start; *i1 = end;
}


typedef struct {
    const Mat *A;
    const Vec *x;
    Vec *y;
    int tid, nt;
} MVJob;

static void* mv_worker(void *p) {
    MVJob *j = (MVJob*)p;
    size_t i0, i1;
    row_range(j->A->rows, j->tid, j->nt, &i0, &i1);

    const size_t n = j->A->cols;
    for (size_t i = i0; i < i1; i++) {
        if (g_stop) break;
        double sum = 0.0;
        const double *row = &j->A->data[i * n];
        for (size_t k = 0; k < n; k++) sum += row[k] * j->x->data[k];
        j->y->data[i] = sum;
    }
    return NULL;
}

int mv_mt(const Mat *A, const Vec *x, Vec *y, KCfg cfg) {
    if (!A || !x || !y || !A->data || !x->data || !y->data) return -1;
    if (A->cols != x->len || A->rows != y->len) return -1;
    int nt = cfg.nt <= 0 ? 1 : cfg.nt;

    pthread_t *ths = (pthread_t*)calloc((size_t)nt, sizeof(pthread_t));
    MVJob  *jobs = (MVJob*)calloc((size_t)nt, sizeof(MVJob));
    if (!ths || !jobs) { free(ths); free(jobs); return -1; }

    for (int t = 0; t < nt; t++) {
        jobs[t] = (MVJob){ .A=A, .x=x, .y=y, .tid=t, .nt=nt };
        if (pthread_create(&ths[t], NULL, mv_worker, &jobs[t]) != 0) {
            nt = t; break;
        }
    }
    for (int t = 0; t < nt; t++) pthread_join(ths[t], NULL);

    free(ths); free(jobs);
    return 0;
}


typedef struct {
    const Mat *A;
    const Mat *B;
    Mat *C;
    int tid, nt;
    int tile;
} MMJob;

static inline size_t min_sz(size_t a, size_t b) { return a < b ? a : b; }

static void* mm_worker(void *p) {
    MMJob *j = (MMJob*)p;
    const Mat *A = j->A;
    const Mat *B = j->B;
    Mat *C = j->C;

    size_t i0, i1;
    row_range(A->rows, j->tid, j->nt, &i0, &i1);

    size_t M = A->rows, K = A->cols, N = B->cols;
    (void)M;

    if (j->tile <= 0) {
        for (size_t i = i0; i < i1; i++) {
            if (g_stop) break;
            for (size_t k = 0; k < K; k++) {
                double a = A->data[i*K + k];
                const double *brow = &B->data[k*N];
                double *crow = &C->data[i*N];
                for (size_t jj = 0; jj < N; jj++) crow[jj] += a * brow[jj];
            }
        }
        return NULL;
    }

    size_t T = (size_t)j->tile;
    for (size_t i = i0; i < i1; i++) {
        if (g_stop) break;
        for (size_t jj = 0; jj < N; jj += T) {
            size_t jend = min_sz(jj + T, N);
            for (size_t kk = 0; kk < K; kk += T) {
                size_t kend = min_sz(kk + T, K);
                for (size_t k = kk; k < kend; k++) {
                    double a = A->data[i*K + k];
                    const double *brow = &B->data[k*N];
                    double *crow = &C->data[i*N];
                    for (size_t col = jj; col < jend; col++) {
                        crow[col] += a * brow[col];
                    }
                }
            }
        }
    }
    return NULL;
}

int mm_mt(const Mat *A, const Mat *B, Mat *C, KCfg cfg) {
    if (!A || !B || !C || !A->data || !B->data || !C->data) return -1;
    if (A->cols != B->rows) return -1;
    if (C->rows != A->rows || C->cols != B->cols) return -1;

    memset(C->data, 0, sizeof(double) * C->rows * C->cols);

    int nt = cfg.nt <= 0 ? 1 : cfg.nt;
    pthread_t *ths = (pthread_t*)calloc((size_t)nt, sizeof(pthread_t));
    MMJob  *jobs = (MMJob*)calloc((size_t)nt, sizeof(MMJob));
    if (!ths || !jobs) { free(ths); free(jobs); return -1; }

    for (int t = 0; t < nt; t++) {
        jobs[t] = (MMJob){ .A=A, .B=B, .C=C, .tid=t, .nt=nt, .tile=cfg.tile };
        if (pthread_create(&ths[t], NULL, mm_worker, &jobs[t]) != 0) {
            nt = t; break;
        }
    }
    for (int t = 0; t < nt; t++) pthread_join(ths[t], NULL);

    free(ths); free(jobs);
    return 0;
}


typedef struct {
    const Vec *x;
    const Vec *y;
    double partial;
    int tid, nt;
} DotJob;

static void* dt_worker(void *p) {
    DotJob *j = (DotJob*)p;
    size_t n = j->x->len;
    size_t i0, i1;
    row_range(n, j->tid, j->nt, &i0, &i1);
    double s = 0.0;
    for (size_t i = i0; i < i1; i++) {
        if (g_stop) break;
        s += j->x->data[i] * j->y->data[i];
    }
    j->partial = s;
    return NULL;
}

int dt_mt(const Vec *x, const Vec *y, double *out, int nt) {
    if (!x || !y || !out || !x->data || !y->data) return -1;
    if (x->len != y->len) return -1;
    nt = nt <= 0 ? 1 : nt;

    pthread_t *ths = (pthread_t*)calloc((size_t)nt, sizeof(pthread_t));
    DotJob   *jobs = (DotJob*)calloc((size_t)nt, sizeof(DotJob));
    if (!ths || !jobs) { free(ths); free(jobs); return -1; }

    for (int t = 0; t < nt; t++) {
        jobs[t] = (DotJob){ .x=x, .y=y, .partial=0.0, .tid=t, .nt=nt };
        if (pthread_create(&ths[t], NULL, dt_worker, &jobs[t]) != 0) {
            nt = t; break;
        }
    }
    double sum = 0.0;
    for (int t = 0; t < nt; t++) {
        pthread_join(ths[t], NULL);
        sum += jobs[t].partial;
    }
    *out = sum;

    free(ths); free(jobs);
    return 0;
}

typedef struct {
    double a;
    const Vec *x;
    Vec *y;
    int tid, nt;
} AXJob;

static void* ax_worker(void *p) {
    AXJob *j = (AXJob*)p;
    size_t n = j->x->len;
    size_t i0, i1;
    row_range(n, j->tid, j->nt, &i0, &i1);
    for (size_t i = i0; i < i1; i++) {
        if (g_stop) break;
        j->y->data[i] = j->a * j->x->data[i] + j->y->data[i];
    }
    return NULL;
}

int ax_mt(double a, const Vec *x, Vec *y, int nt) {
    if (!x || !y || !x->data || !y->data) return -1;
    if (x->len != y->len) return -1;
    nt = nt <= 0 ? 1 : nt;

    pthread_t *ths = (pthread_t*)calloc((size_t)nt, sizeof(pthread_t));
    AXJob  *jobs = (AXJob*)calloc((size_t)nt, sizeof(AXJob));
    if (!ths || !jobs) { free(ths); free(jobs); return -1; }

    for (int t = 0; t < nt; t++) {
        jobs[t] = (AXJob){ .a=a, .x=x, .y=y, .tid=t, .nt=nt };
        if (pthread_create(&ths[t], NULL, ax_worker, &jobs[t]) != 0) {
            nt = t; break;
        }
    }
    for (int t = 0; t < nt; t++) pthread_join(ths[t], NULL);

    free(ths); free(jobs);
    return 0;
}
