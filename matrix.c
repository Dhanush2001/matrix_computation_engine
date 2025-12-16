#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

static int read_u64(FILE *f, uint64_t *x) {
    return fread(x, sizeof(uint64_t), 1, f) == 1 ? 0 : -1;
}
static int write_u64(FILE *f, uint64_t x) {
    return fwrite(&x, sizeof(uint64_t), 1, f) == 1 ? 0 : -1;
}

Mat m_alloc(size_t r, size_t c) {
    Mat m = {0};
    m.rows = r; m.cols = c;
    if (r == 0 || c == 0) return m;
    m.data = (double*)calloc(r * c, sizeof(double));
    if (!m.data) { m.rows = m.cols = 0; }
    return m;
}

void m_free(Mat *m) {
    if (!m) return;
    free(m->data);
    m->data = NULL;
    m->rows = m->cols = 0;
}

Vec v_alloc(size_t n) {
    Vec v = {0};
    v.len = n;
    if (n == 0) return v;
    v.data = (double*)calloc(n, sizeof(double));
    if (!v.data) v.len = 0;
    return v;
}

void v_free(Vec *v) {
    if (!v) return;
    free(v->data);
    v->data = NULL;
    v->len = 0;
}

int m_load(const char *path, FileFmt fmt, Mat *out) {
    if (!path || !out) return -1;
    FILE *f = fopen(path, fmt == FMT_BIN ? "rb" : "r");
    if (!f) { perror("fopen"); return -1; }

    Mat m = {0};

    if (fmt == FMT_TEXT) {
        size_t r, c;
        if (fscanf(f, "%zu %zu", &r, &c) != 2) { fclose(f); return -1; }
        m = m_alloc(r, c);
        if (!m.data) { fclose(f); return -1; }
        for (size_t i = 0; i < r * c; i++) {
            if (fscanf(f, "%lf", &m.data[i]) != 1) {
                m_free(&m); fclose(f); return -1;
            }
        }
    } else {
        uint64_t r64, c64;
        if (read_u64(f, &r64) || read_u64(f, &c64)) { fclose(f); return -1; }
        if (r64 == 0 || c64 == 0) { fclose(f); return -1; }
        m = m_alloc((size_t)r64, (size_t)c64);
        if (!m.data) { fclose(f); return -1; }
        size_t n = m.rows * m.cols;
        if (fread(m.data, sizeof(double), n, f) != n) {
            m_free(&m); fclose(f); return -1;
        }
    }

    fclose(f);
    *out = m;
    return 0;
}

int m_save(const char *path, FileFmt fmt, const Mat *m) {
    if (!path || !m || !m->data) return -1;
    FILE *f = fopen(path, fmt == FMT_BIN ? "wb" : "w");
    if (!f) { perror("fopen"); return -1; }

    if (fmt == FMT_TEXT) {
        fprintf(f, "%zu %zu\n", m->rows, m->cols);
        for (size_t i = 0; i < m->rows; i++) {
            for (size_t j = 0; j < m->cols; j++) {
                fprintf(f, "%.17g%s", m_get(m, i, j), (j + 1 == m->cols) ? "" : " ");
            }
            fprintf(f, "\n");
        }
    } else {
        if (write_u64(f, (uint64_t)m->rows) || write_u64(f, (uint64_t)m->cols)) {
            fclose(f); return -1;
        }
        size_t n = m->rows * m->cols;
        if (fwrite(m->data, sizeof(double), n, f) != n) { fclose(f); return -1; }
    }

    fclose(f);
    return 0;
}

int v_load(const char *path, FileFmt fmt, Vec *out) {
    if (!path || !out) return -1;
    FILE *f = fopen(path, fmt == FMT_BIN ? "rb" : "r");
    if (!f) { perror("fopen"); return -1; }

    Vec v = {0};

    if (fmt == FMT_TEXT) {
        size_t n;
        if (fscanf(f, "%zu", &n) != 1) { fclose(f); return -1; }
        v = v_alloc(n);
        if (!v.data) { fclose(f); return -1; }
        for (size_t i = 0; i < n; i++) {
            if (fscanf(f, "%lf", &v.data[i]) != 1) {
                v_free(&v); fclose(f); return -1;
            }
        }
    } else {
        uint64_t n64;
        if (read_u64(f, &n64)) { fclose(f); return -1; }
        v = v_alloc((size_t)n64);
        if (!v.data) { fclose(f); return -1; }
        if (fread(v.data, sizeof(double), v.len, f) != v.len) {
            v_free(&v); fclose(f); return -1;
        }
    }

    fclose(f);
    *out = v;
    return 0;
}

int v_save(const char *path, FileFmt fmt, const Vec *v) {
    if (!path || !v || !v->data) return -1;
    FILE *f = fopen(path, fmt == FMT_BIN ? "wb" : "w");
    if (!f) { perror("fopen"); return -1; }

    if (fmt == FMT_TEXT) {
        fprintf(f, "%zu\n", v->len);
        for (size_t i = 0; i < v->len; i++) {
            fprintf(f, "%.17g\n", v->data[i]);
        }
    } else {
        if (write_u64(f, (uint64_t)v->len)) { fclose(f); return -1; }
        if (fwrite(v->data, sizeof(double), v->len, f) != v->len) { fclose(f); return -1; }
    }

    fclose(f);
    return 0;
}
