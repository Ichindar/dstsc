CC=g++
CFLAGS= -O3 
VPATH= ../SAIS-SK
default: test

sk-sain: sk-sain.c
	$(CC) $(CFLAGS) -c $< -o sk-sain.o

gt-alloc: gt-alloc.c
	$(CC) $(CFLAGS) -c $< -o gt-alloc.o

mmap: mmap.c
	$(CC) $(CFLAGS) -c $< -o mmap.o

fileopen: fileopen.c
	$(CC) $(CFLAGS) -c $< -o fileopen.o

lzOG: lzOG.cpp
	$(CC) $(CFLAGS) -c $< -o lzOG.o

main: main.cpp
	$(CC) $(CFLAGS) -c $< -o main.o

test: main sk-sain gt-alloc mmap fileopen lzOG
	$(CC) $(CFLAGS) main.o sk-sain.o gt-alloc.o mmap.o fileopen.o lzOG.o -o lz_factor_file

clean:
	rm *.o
