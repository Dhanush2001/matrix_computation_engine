# Multithreaded Linear Algebra Benchmark

A high-performance C program for benchmarking common linear algebra operations with multithreading support using pthreads.

## Overview

This program implements and benchmarks four fundamental linear algebra operations:

- **Matrix Multiplication (mm)**: `C = A × B`
- **Matrix-Vector Product (mv)**: `y = A × x`
- **Dot Product (dot)**: `result = x · y`
- **AXPY (axpy)**: `y = α×x + y` (scaled vector addition)

Each operation is benchmarked with both single-threaded and multi-threaded implementations to measure parallel performance, speedup, and efficiency.

## Features

- **Multithreaded computation** using POSIX threads (pthreads)
- **Tiled matrix multiplication** for improved cache performance
- **High-precision timing** using `clock_gettime(CLOCK_MONOTONIC)`
- **Performance metrics**: execution time, GFLOPS, speedup, parallel efficiency
- **Flexible I/O**: supports both text and binary file formats
- **Graceful interruption** via SIGINT (Ctrl+C)
- **Batch mode**: run all operations sequentially with `--op all`

## Requirements

- **Compiler**: GCC with C11 support
- **OS**: POSIX-compliant system (Linux, macOS, Unix)
- **Libraries**: pthread

## Building

```bash
make
```

This compiles the program with `-O3 -march=native` optimizations for maximum performance.

To clean build artifacts:
```bash
make clean
```

## Usage

### Basic Syntax

```bash
./main --op <operation> --format <text|bin> --threads <N> --out <output_path> [options]
```

### Operations

#### Matrix Multiplication
```bash
./main --op mm --format text \
       --A A.txt --B B.txt \
       --threads 4 --repeat 10 --tile 128 \
       --out /dev/null
```

#### Matrix-Vector Product
```bash
./main --op mv --format text \
       --A A.txt --x x.txt \
       --threads 4 --repeat 10 \
       --out /dev/null
```

#### Dot Product
```bash
./main --op dot --format text \
       --x x.txt --y y.txt \
       --threads 4 --repeat 10 \
       --out /dev/null
```

#### AXPY (Scaled Vector Addition)
```bash
./main --op axpy --format text \
       --x x.txt --y y.txt --alpha 2.0 \
       --threads 4 --repeat 10 \
       --out /dev/null
```

#### Run All Operations
```bash
./main --op all --format text \
       --A A.txt --B B.txt --x x.txt --y y.txt \
       --alpha 2.0 --threads 4 --repeat 10 --tile 128 \
       --out /dev/null
```

### Command-Line Options

| Option | Description | Required |
|--------|-------------|----------|
| `--op` | Operation: `mm`, `mv`, `dot`, `axpy`, or `all` | Yes |
| `--format` | File format: `text` or `bin` | Yes |
| `--threads` | Number of threads to use | Yes |
| `--out` | Output path (use `/dev/null` for no output) | Yes |
| `--A` | Path to matrix A (for `mm`, `mv`) | Conditional |
| `--B` | Path to matrix B (for `mm`) | Conditional |
| `--x` | Path to vector x (for `mv`, `dot`, `axpy`) | Conditional |
| `--y` | Path to vector y (for `dot`, `axpy`) | Conditional |
| `--alpha` | Scalar value for AXPY (default: 1.0) | Optional |
| `--repeat` | Number of repetitions for timing (default: 1) | Optional |
| `--tile` | Tile size for matrix multiplication (default: 64) | Optional |
| `--help` | Display help message | Optional |

## File Formats

### Text Format

**Matrix** (rows cols followed by data):
```
3 2
1.0 2.0
3.0 4.0
5.0 6.0
```

**Vector** (length followed by data):
```
3
1.0
2.0
3.0
```

### Binary Format

- **Matrix**: `uint64_t rows`, `uint64_t cols`, followed by `rows×cols` doubles
- **Vector**: `uint64_t length`, followed by `length` doubles

All values are in little-endian format.

## Output

The program outputs performance metrics to stdout showing:
- Operation name
- Matrix/vector dimensions
- Number of threads used
- Average execution time
- GFLOPS (billion floating-point operations per second)
- Speedup relative to single-threaded execution
- Parallel efficiency (speedup / threads × 100%)

Each operation also displays a preview of the result (matrix/vector values or scalar).

## Project Structure

```
.
├── main.c          # Main program and benchmark orchestration
├── kernels.c       # Multithreaded kernel implementations
├── kernels.h       # Kernel function prototypes
├── matrix.c        # Matrix/vector I/O and memory management
├── matrix.h        # Data structures and interfaces
├── bench.c         # High-precision timing utilities
├── bench.h         # Timing function prototypes
├── Makefile        # Build configuration
├── A.txt, B.txt    # Sample matrix data files
├── x.txt, y.txt    # Sample vector data files
└── *_l.txt         # Large test data files
```

## Implementation Details

### Parallelization Strategy

- **Work Distribution**: Row-based partitioning with load balancing
- **Thread Safety**: Each thread operates on independent data regions
- **Synchronization**: Uses pthread join for barrier synchronization

### Matrix Multiplication Optimization

- **Tiling**: Configurable tile size for improved cache locality
- **Memory Access**: Row-major storage with cache-friendly access patterns
- **Algorithm**: `C[i][j] += A[i][k] × B[k][j]` with tiled blocking

### Performance Considerations

1. **Compiler Optimizations**: Built with `-O3 -march=native`
2. **Cache Performance**: Tiling reduces cache misses
3. **Thread Count**: Use a portable method to detect logical CPUs (e.g., `nproc`, `sysctl -n hw.logicalcpu`, or `getconf _NPROCESSORS_ONLN`); example snippets below.
4. **Repetitions**: Use `--repeat` for more accurate timing measurements

## Example Sessions

### Running with Small Test Files (A.txt, B.txt, x.txt, y.txt)

```bash
# Detect logical CPU count (portable)
THREADS=$(
       { command -v nproc >/dev/null 2>&1 && nproc; } ||
       { command -v sysctl >/dev/null 2>&1 && sysctl -n hw.logicalcpu 2>/dev/null; } ||
       { command -v getconf >/dev/null 2>&1 && getconf _NPROCESSORS_ONLN; } ||
       echo 1
)

# Run comprehensive benchmark
./main --op all --format text \
       --A A.txt --B B.txt --x x.txt --y y.txt \
       --alpha 2.0 --threads $THREADS --repeat 5 --tile 128 \
       --out /dev/null
```

### Running with Large Test Files (A_l.txt, B_l.txt, x_l.txt, y_l.txt)

#### Single-Threaded Execution

Run each operation individually with 1 thread:

```bash
# Matrix multiplication (single-threaded)
./main --op mm --format text \
       --A A_l.txt --B B_l.txt \
       --threads 1 --repeat 3 --tile 128 \
       --out /dev/null

# Matrix-vector product (single-threaded)
./main --op mv --format text \
       --A A_l.txt --x x_l.txt \
       --threads 1 --repeat 3 \
       --out /dev/null

# Dot product (single-threaded)
./main --op dot --format text \
       --x x_l.txt --y y_l.txt \
       --threads 1 --repeat 3 \
       --out /dev/null

# AXPY (single-threaded)
./main --op axpy --format text \
       --x x_l.txt --y y_l.txt --alpha 2.0 \
       --threads 1 --repeat 3 \
       --out /dev/null
```

#### Multi-Threaded Execution

Run each operation with maximum available threads:

```bash
# Detect logical CPU count (portable)
THREADS=$(
       { command -v nproc >/dev/null 2>&1 && nproc; } ||
       { command -v sysctl >/dev/null 2>&1 && sysctl -n hw.logicalcpu 2>/dev/null; } ||
       { command -v getconf >/dev/null 2>&1 && getconf _NPROCESSORS_ONLN; } ||
       echo 1
)

# Matrix multiplication (multi-threaded)
./main --op mm --format text \
       --A A_l.txt --B B_l.txt \
       --threads $THREADS --repeat 3 --tile 128 \
       --out /dev/null

# Matrix-vector product (multi-threaded)
./main --op mv --format text \
       --A A_l.txt --x x_l.txt \
       --threads $THREADS --repeat 3 \
       --out /dev/null

# Dot product (multi-threaded)
./main --op dot --format text \
       --x x_l.txt --y y_l.txt \
       --threads $THREADS --repeat 3 \
       --out /dev/null

# AXPY (multi-threaded)
./main --op axpy --format text \
       --x x_l.txt --y y_l.txt --alpha 2.0 \
       --threads $THREADS --repeat 3 \
       --out /dev/null
```

#### Run All Operations (Large Files)

Execute all four operations sequentially:

```bash
# Single-threaded (baseline)
./main --op all --format text \
       --A A_l.txt --B B_l.txt --x x_l.txt --y y_l.txt \
       --alpha 2.0 --threads 1 --repeat 3 --tile 128 \
       --out /dev/null

# Multi-threaded (parallel)
# Using THREADS from the detection snippet above
./main --op all --format text \
       --A A_l.txt --B B_l.txt --x x_l.txt --y y_l.txt \
       --alpha 2.0 --threads $THREADS --repeat 3 --tile 128 \
       --out /dev/null
```

#### Performance Comparison Script

Compare single vs multi-threaded performance:

```bash
#!/bin/bash
THREADS=$(
       { command -v nproc >/dev/null 2>&1 && nproc; } ||
       { command -v sysctl >/dev/null 2>&1 && sysctl -n hw.logicalcpu 2>/dev/null; } ||
       { command -v getconf >/dev/null 2>&1 && getconf _NPROCESSORS_ONLN; } ||
       echo 1
)

echo "=== Single-Threaded Baseline ==="
./main --op all --format text \
       --A A_l.txt --B B_l.txt --x x_l.txt --y y_l.txt \
       --alpha 2.0 --threads 1 --repeat 3 --tile 128 \
       --out /dev/null

echo ""
echo "=== Multi-Threaded ($THREADS threads) ==="
./main --op all --format text \
       --A A_l.txt --B B_l.txt --x x_l.txt --y y_l.txt \
       --alpha 2.0 --threads $THREADS --repeat 3 --tile 128 \
       --out /dev/null
```

The program displays performance metrics for each operation with both 1 thread and N threads, showing execution time, GFLOPS, speedup, and parallel efficiency.

## Interruption Handling

Press `Ctrl+C` to gracefully interrupt any running benchmark. The program will:
- Stop current computation
- Clean up allocated memory
- Exit with code 2


