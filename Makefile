CC=gcc
CFLAGS=-O3 -march=native -Wall -Wextra -std=c11 -pthread
LDFLAGS=-pthread

OBJS=main.o matrix.o kernels.o bench.o

all: main

main: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o main
