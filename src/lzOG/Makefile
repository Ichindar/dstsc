CC=g++
CFLAGS= -O3 

default: test

sk-sain: ../SAIS-SK/sk-sain.c
	$(CC) $(CFLAGS) -c ../SAIS-SK/sk-sain.c -o sk-sain.o

gt-alloc: ../SAIS-SK/gt-alloc.c
	$(CC) $(CFLAGS) -c ../SAIS-SK/gt-alloc.c -o gt-alloc.o

mmap: ../SAIS-SK/mmap.c
	$(CC) $(CFLAGS) -c ../SAIS-SK/mmap.c -o mmap.o

fileopen: ../SAIS-SK/fileopen.c
	$(CC) $(CFLAGS) -c ../SAIS-SK/fileopen.c -o fileopen.o

lzOG: lzOG.cpp
	$(CC) $(CFLAGS) -c lzOG.cpp -o lzOG.o

main: main.cpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

test: main sk-sain gt-alloc mmap fileopen lzOG
	$(CC) $(CFLAGS) main.o sk-sain.o gt-alloc.o mmap.o fileopen.o lzOG.o -o lz_factor

clean:
	rm *.o
