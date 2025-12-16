#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    size_t rows;
    size_t cols;
    double *data;
} Mat;

typedef struct {
    size_t len;
    double *data;
} Vec;

typedef enum { FMT_TEXT, FMT_BIN } FileFmt;

Mat  m_alloc(size_t r, size_t c);
void m_free(Mat *m);
Vec  v_alloc(size_t n);
void v_free(Vec *v);

static inline double m_get(const Mat *A, size_t i, size_t j) {
    return A->data[i * A->cols + j];
}
static inline void m_set(Mat *A, size_t i, size_t j, double x) {
    A->data[i * A->cols + j] = x;
}

int m_load(const char *path, FileFmt fmt, Mat *out);
int m_save(const char *path, FileFmt fmt, const Mat *m);

int v_load(const char *path, FileFmt fmt, Vec *out);
int v_save(const char *path, FileFmt fmt, const Vec *v);

#endif
