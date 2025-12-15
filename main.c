#define _POSIX_C_SOURCE 200809L
#include "matrix.h"
#include "kernels.h"
#include "bench.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <signal.h>
#include <limits.h>

volatile sig_atomic_t g_stop = 0;

static void on_sigint(int signo) {
    (void)signo;
    g_stop = 1;
}

typedef enum { OP_NONE, OP_MM, OP_MV, OP_DOT, OP_AXPY, OP_ALL } Op;

static FileFmt parse_fmt(const char *s) {
    if (!s) return FMT_TEXT;
    if (strcmp(s, "text") == 0) return FMT_TEXT;
    if (strcmp(s, "bin")  == 0) return FMT_BIN;
    return FMT_TEXT;
}

static Op parse_op(const char *s) {
    if (!s) return OP_NONE;
    if (strcmp(s, "mm") == 0)   return OP_MM;
    if (strcmp(s, "mv") == 0)   return OP_MV;
    if (strcmp(s, "dot") == 0)  return OP_DOT;
    if (strcmp(s, "axpy") == 0) return OP_AXPY;
    if (strcmp(s, "all") == 0)  return OP_ALL;
    return OP_NONE;
}

static const char* op_name(Op op) {
    switch (op) {
        case OP_MM:   return "mm";
        case OP_MV:   return "mv";
        case OP_DOT:  return "dot";
        case OP_AXPY: return "axpy";
        case OP_ALL:  return "all";
        default:      return "unknown";
    }
}

static void usage(const char *argv0) {
    fprintf(stderr,
        "Usage:\n"
    "  %s --op {mm|mv|dot|axpy|all} --format {text|bin} --threads N --out OUT [--repeat R]\n"
        "\n"
        "Ops:\n"
        "  mm:   --A Afile --B Bfile [--tile T]\n"
        "  mv:   --A Afile --x xfile\n"
        "  dot:  --x xfile --y yfile\n"
        "  axpy: --alpha a --x xfile --y yfile\n"
        "\n"
        "all:\n"
        "  Runs mm -> mv -> dot -> axpy in that order.\n"
        "  It will SKIP ops whose required inputs are missing.\n"
        "\n"
    "Output: Results printed to terminal only.\n"
        "\n",
        argv0);
}

static double gflops_for(Op op, size_t m, size_t n, size_t k, size_t len, double sec) {
    if (sec <= 0) return 0.0;
    double flops = 0.0;
    switch (op) {
        case OP_MM:   flops = 2.0 * (double)m * (double)n * (double)k; break;
        case OP_MV:   flops = 2.0 * (double)m * (double)n; break;
        case OP_DOT:  flops = 2.0 * (double)len; break;
        case OP_AXPY: flops = 2.0 * (double)len; break;
        default: break;
    }
    return (flops / 1e9) / sec;
}

static void print_csv_header_stdout(void) {
    printf("op,m,n,k,threads,seconds,gflops,speedup,efficiency,format\n");
}

static void print_csv_row_stdout(const char *op,
                                 size_t m, size_t n, size_t k,
                                 int threads,
                                 double seconds,
                                 double gflops,
                                 double speedup,
                                 double efficiency,
                                 const char *fmt_name) {
    printf("%s,%zu,%zu,%zu,%d,%.9f,%.6f,%.4f,%.2f,%s\n",
           op, m, n, k, threads, seconds, gflops, speedup, efficiency, fmt_name ? fmt_name : "");
}




static void print_vector_preview(const Vec *v, size_t maxn) {
    size_t n = v->len < maxn ? v->len : maxn;
    printf("[");
    for (size_t i = 0; i < n; i++) {
        printf("%.6g%s", v->data[i], (i + 1 == n) ? "" : ", ");
    }
    if (v->len > n) printf(", ...");
    printf("]\n");
}

static void print_matrix_preview(const Mat *m, size_t maxr, size_t maxc) {
    size_t r = m->rows < maxr ? m->rows : maxr;
    size_t c = m->cols < maxc ? m->cols : maxc;
    for (size_t i = 0; i < r; i++) {
        printf("[");
        for (size_t j = 0; j < c; j++) {
            printf("%.6g%s", m->data[i*m->cols + j], (j + 1 == c) ? "" : ", ");
        }
        if (m->cols > c) printf(", ...");
        printf("]\n");
    }
    if (m->rows > r) printf("...\n");
}

static int prep_logging(const char *out_base, const char *op, char *csv_path, size_t csv_path_sz) {
    (void)out_base;  // Unused: CSV files disabled
    (void)csv_path_sz;
    csv_path[0] = '\0';  // Clear path
    printf("\n[%s] Results:\n", op);
    print_csv_header_stdout();
    return 0;
}

static void run_bench(const char *op, size_t m, size_t n, size_t k, size_t len,
                      int nt1, int ntN, double sec1, double secN,
                      const char *csv_path, const char *fmt_name) {
    (void)csv_path;  // Unused: CSV output disabled
    double g1 = gflops_for((strcmp(op,"mm")==0 ? OP_MM : strcmp(op,"mv")==0 ? OP_MV : strcmp(op,"dot")==0 ? OP_DOT : OP_AXPY), m, n, k, len, sec1);
    double gN = gflops_for((strcmp(op,"mm")==0 ? OP_MM : strcmp(op,"mv")==0 ? OP_MV : strcmp(op,"dot")==0 ? OP_DOT : OP_AXPY), m, n, k, len, secN);
    double sp = (secN > 0.0) ? (sec1 / secN) : 0.0;
    double eff = (ntN > 0) ? (100.0 * sp / (double)ntN) : 0.0;
    print_csv_row_stdout(op, m, n, k, nt1, sec1, g1, 1.0, 100.0, fmt_name);
    print_csv_row_stdout(op, m, n, k, ntN, secN, gN, sp, eff, fmt_name);
}

static int do_mm(const char *out_base, FileFmt fmt,
                 const char *Apath, const char *Bpath,
                 int nt, int rep, int tile) {
    if (!Apath || !Bpath) {
        printf("\n[mm] Skipped: need --A and --B\n");
        return 0;
    }

    char csv_path[PATH_MAX];
    if (prep_logging(out_base, "mm", csv_path, sizeof(csv_path)) != 0) return -1;

    const char *fmt_name = (fmt == FMT_BIN) ? "bin" : "text";
    KCfg cfg = { .nt = nt, .tile = tile };
    KCfg cfg1 = cfg; cfg1.nt = 1;

    Mat A={0}, B={0};
    if (m_load(Apath, fmt, &A) || m_load(Bpath, fmt, &B)) {
        fprintf(stderr, "[mm] Failed to load A/B\n");
        return -1;
    }
    if (A.cols != B.rows) {
        fprintf(stderr, "[mm] Dimension mismatch: A=%zux%zu, B=%zux%zu\n", A.rows, A.cols, B.rows, B.cols);
        m_free(&A); m_free(&B);
        return -1;
    }

    Mat C1 = m_alloc(A.rows, B.cols);
    Mat CN = m_alloc(A.rows, B.cols);
    if (!C1.data || !CN.data) {
        fprintf(stderr, "[mm] Allocation failure\n");
        m_free(&A); m_free(&B);
        m_free(&C1); m_free(&CN);
        return -1;
    }

    double total1 = 0.0;
    for (int r = 0; r < rep && !g_stop; r++) {
        double t0 = now_s();
        mm_mt(&A, &B, &C1, cfg1);
        total1 += now_s() - t0;
    }
    if (g_stop) {
        fprintf(stderr, "[mm] Interrupted\n");
        m_free(&A); m_free(&B); m_free(&C1); m_free(&CN);
        return 2;
    }

    double totalN = 0.0;
    for (int r = 0; r < rep && !g_stop; r++) {
        double t0 = now_s();
        mm_mt(&A, &B, &CN, cfg);
        totalN += now_s() - t0;
    }
    if (g_stop) {
        fprintf(stderr, "[mm] Interrupted\n");
        m_free(&A); m_free(&B); m_free(&C1); m_free(&CN);
        return 2;
    }

    double sec1 = total1 / (double)rep;
    double secN = totalN / (double)rep;
    run_bench("mm", A.rows, B.cols, A.cols, 0, 1, nt, sec1, secN, csv_path, fmt_name);

    printf("C preview (top-left):\n");
    print_matrix_preview(&CN, 4, 4);

    m_free(&A); m_free(&B); m_free(&C1); m_free(&CN);
    return g_stop ? 2 : 0;
}

static int do_mv(const char *out_base, FileFmt fmt,
                 const char *Apath, const char *xpath,
                 int nt, int rep) {
    if (!Apath || !xpath) {
        printf("\n[mv] Skipped: need --A and --x\n");
        return 0;
    }

    char csv_path[PATH_MAX];
    if (prep_logging(out_base, "mv", csv_path, sizeof(csv_path)) != 0) return -1;

    const char *fmt_name = (fmt == FMT_BIN) ? "bin" : "text";
    KCfg cfg = { .nt = nt, .tile = 0 };
    KCfg cfg1 = cfg; cfg1.nt = 1;

    Mat A={0};
    Vec x={0}, y1={0}, yN={0};

    if (m_load(Apath, fmt, &A) || v_load(xpath, fmt, &x)) {
        fprintf(stderr, "[mv] Failed to load A/x\n");
        return -1;
    }
    if (A.cols != x.len) {
        fprintf(stderr, "[mv] Dimension mismatch: A=%zux%zu, x=%zu\n", A.rows, A.cols, x.len);
        m_free(&A); v_free(&x);
        return -1;
    }

    y1 = v_alloc(A.rows);
    yN = v_alloc(A.rows);

    double total1=0.0;
    for (int r=0; r<rep && !g_stop; r++) {
        double t0 = now_s();
        mv_mt(&A, &x, &y1, cfg1);
        total1 += now_s() - t0;
    }
    if (g_stop) {
        fprintf(stderr, "[mv] Interrupted\n");
        m_free(&A); v_free(&x); v_free(&y1); v_free(&yN);
        return 2;
    }

    double totalN=0.0;
    for (int r=0; r<rep && !g_stop; r++) {
        double t0 = now_s();
        mv_mt(&A, &x, &yN, cfg);
        totalN += now_s() - t0;
    }
    if (g_stop) {
        fprintf(stderr, "[mv] Interrupted\n");
        m_free(&A); v_free(&x); v_free(&y1); v_free(&yN);
        return 2;
    }

    double sec1 = total1 / (double)rep;
    double secN = totalN / (double)rep;
    run_bench("mv", A.rows, A.cols, 0, 0, 1, nt, sec1, secN, csv_path, fmt_name);

    printf("y preview:\n");
    print_vector_preview(&yN, 10);


    m_free(&A); v_free(&x); v_free(&y1); v_free(&yN);
    return g_stop ? 2 : 0;
}

static int do_dot(const char *out_base, FileFmt fmt,
                  const char *xpath, const char *ypath,
                  int nt, int rep) {
    if (!xpath || !ypath) {
        printf("\n[dot] Skipped: need --x and --y\n");
        return 0;
    }

    char csv_path[PATH_MAX];
    if (prep_logging(out_base, "dot", csv_path, sizeof(csv_path)) != 0) return -1;

    const char *fmt_name = (fmt == FMT_BIN) ? "bin" : "text";
    Vec x={0}, y={0};

    if (v_load(xpath, fmt, &x) || v_load(ypath, fmt, &y)) {
        fprintf(stderr, "[dot] Failed to load x/y\n");
        return -1;
    }
    if (x.len != y.len) {
        fprintf(stderr, "[dot] Dimension mismatch: x=%zu, y=%zu\n", x.len, y.len);
        v_free(&x); v_free(&y);
        return -1;
    }

    double out1=0.0, outN=0.0;
    double total1=0.0;
    for (int r=0; r<rep && !g_stop; r++) {
        double t0 = now_s();
        dt_mt(&x, &y, &out1, 1);
        total1 += now_s() - t0;
    }
    if (g_stop) {
        fprintf(stderr, "[dot] Interrupted\n");
        v_free(&x); v_free(&y);
        return 2;
    }

    double totalN=0.0;
    for (int r=0; r<rep && !g_stop; r++) {
        double t0 = now_s();
        dt_mt(&x, &y, &outN, nt);
        totalN += now_s() - t0;
    }
    if (g_stop) {
        fprintf(stderr, "[dot] Interrupted\n");
        v_free(&x); v_free(&y);
        return 2;
    }

    double sec1 = total1 / (double)rep;
    double secN = totalN / (double)rep;
    run_bench("dot", 0, 0, 0, x.len, 1, nt, sec1, secN, csv_path, fmt_name);

    printf("dot = %.17g (1t), %.17g (%dt)\n", out1, outN, nt);

    v_free(&x); v_free(&y);
    return g_stop ? 2 : 0;
}

static int do_axpy(const char *out_base, FileFmt fmt,
                   double a, const char *xpath, const char *ypath,
                   int nt, int rep) {
    if (!xpath || !ypath) {
        printf("\n[axpy] Skipped: need --x and --y\n");
        return 0;
    }

    char csv_path[PATH_MAX];
    if (prep_logging(out_base, "axpy", csv_path, sizeof(csv_path)) != 0) return -1;

    const char *fmt_name = (fmt == FMT_BIN) ? "bin" : "text";
    Vec x={0}, y1={0}, yN={0};

    if (v_load(xpath, fmt, &x) || v_load(ypath, fmt, &y1)) {
        fprintf(stderr, "[axpy] Failed to load x/y\n");
        return -1;
    }
    if (x.len != y1.len) {
        fprintf(stderr, "[axpy] Dimension mismatch: x=%zu, y=%zu\n", x.len, y1.len);
        v_free(&x); v_free(&y1);
        return -1;
    }

    yN = v_alloc(y1.len);
    if (!yN.data) {
        fprintf(stderr, "[axpy] Allocation failure\n");
        v_free(&x); v_free(&y1);
        return -1;
    }
    memcpy(yN.data, y1.data, sizeof(double) * y1.len);

    double total1=0.0;
    for (int r=0; r<rep && !g_stop; r++) {
        Vec ytmp = v_alloc(y1.len);
        memcpy(ytmp.data, y1.data, sizeof(double) * y1.len);
        double t0 = now_s();
        ax_mt(a, &x, &ytmp, 1);
        total1 += now_s() - t0;
        v_free(&ytmp);
    }
    if (g_stop) {
        fprintf(stderr, "[axpy] Interrupted\n");
        v_free(&x); v_free(&y1); v_free(&yN);
        return 2;
    }

    double totalN=0.0;
    for (int r=0; r<rep && !g_stop; r++) {
        memcpy(yN.data, y1.data, sizeof(double) * y1.len);
        double t0 = now_s();
        ax_mt(a, &x, &yN, nt);
        totalN += now_s() - t0;
    }
    if (g_stop) {
        fprintf(stderr, "[axpy] Interrupted\n");
        v_free(&x); v_free(&y1); v_free(&yN);
        return 2;
    }

    double sec1 = total1 / (double)rep;
    double secN = totalN / (double)rep;
    run_bench("axpy", 0, 0, 0, x.len, 1, nt, sec1, secN, csv_path, fmt_name);

    printf("alpha=%.6g, y preview:\n", a);
    print_vector_preview(&yN, 10);

    v_free(&x); v_free(&y1); v_free(&yN);
    return g_stop ? 2 : 0;
}

int main(int argc, char **argv) {
    signal(SIGINT, on_sigint);

    Op op = OP_NONE;
    FileFmt fmt = FMT_TEXT;
    const char *Apath=NULL, *Bpath=NULL, *xpath=NULL, *ypath=NULL;
    const char *out_base=NULL;
    int nt = 1;
    int rep = 1;
    int tile = 64;
    double alpha = 1.0;

    static struct option longopts[] = {
        {"op", required_argument, 0, 'o'},
        {"format", required_argument, 0, 'f'},
        {"A", required_argument, 0, 'A'},
        {"B", required_argument, 0, 'B'},
        {"x", required_argument, 0, 'x'},
        {"y", required_argument, 0, 'y'},
        {"alpha", required_argument, 0, 'a'},
        {"threads", required_argument, 0, 't'},
        {"repeat", required_argument, 0, 'r'},
        {"tile", required_argument, 0, 'T'},
        {"result", required_argument, 0, 'R'},
        {"out", required_argument, 0, 'O'},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "o:f:A:B:x:y:a:t:r:T:R:O:h", longopts, NULL)) != -1) {
        switch (c) {
            case 'o': op = parse_op(optarg); break;
            case 'f': fmt = parse_fmt(optarg); break;
            case 'A': Apath = optarg; break;
            case 'B': Bpath = optarg; break;
            case 'x': xpath = optarg; break;
            case 'y': ypath = optarg; break;
            case 'a': alpha = strtod(optarg, NULL); break;
            case 't': nt = atoi(optarg); break;
            case 'r': rep = atoi(optarg); break;
            case 'T': tile = atoi(optarg); break;
            case 'O': out_base = optarg; break;
            case 'h': usage(argv[0]); return 0;
            default: usage(argv[0]); return 1;
        }
    }

    if (op == OP_NONE || !out_base || nt <= 0 || rep <= 0) {
        usage(argv[0]); return 1;
    }

    if (op == OP_ALL) {
        printf("[Mode] --op all\n");
         printf("[Inputs] A=%s B=%s x=%s y=%s alpha=%.6g threads=%d repeat=%d tile=%d format=%s\n",
               Apath?Apath:"(null)", Bpath?Bpath:"(null)", xpath?xpath:"(null)", ypath?ypath:"(null)",
             alpha, nt, rep, tile, (fmt==FMT_BIN?"bin":"text"));

        int rc;

        rc = do_mm(out_base, fmt, Apath, Bpath, nt, rep, tile);
        if (rc == 2) return 2;
        if (rc < 0) return 1;

        rc = do_mv(out_base, fmt, Apath, xpath, nt, rep);
        if (rc == 2) return 2;
        if (rc < 0) return 1;

        rc = do_dot(out_base, fmt, xpath, ypath, nt, rep);
        if (rc == 2) return 2;
        if (rc < 0) return 1;

        rc = do_axpy(out_base, fmt, alpha, xpath, ypath, nt, rep);
        if (rc == 2) return 2;
        if (rc < 0) return 1;

        return 0;
    }

    if (op == OP_MM)   return do_mm(out_base, fmt, Apath, Bpath, nt, rep, tile) < 0 ? 1 : (g_stop ? 2 : 0);
    if (op == OP_MV)   return do_mv(out_base, fmt, Apath, xpath, nt, rep) < 0 ? 1 : (g_stop ? 2 : 0);
    if (op == OP_DOT)  return do_dot(out_base, fmt, xpath, ypath, nt, rep) < 0 ? 1 : (g_stop ? 2 : 0);
    if (op == OP_AXPY) return do_axpy(out_base, fmt, alpha, xpath, ypath, nt, rep) < 0 ? 1 : (g_stop ? 2 : 0);

    fprintf(stderr, "Unknown op: %s\n", op_name(op));
    return 1;
}
