CC= g++
CFLAGS= -O3 
VPATH= ../../SAIS-SK/src
OBJECTS= main.o sk-sain.o gt-alloc.o mmap.o fileopen.o lzOG.o 

default: test

sk-sain.o: sk-sain.c
	$(CC) $(CFLAGS) -c $< -o sk-sain.o

gt-alloc.o: gt-alloc.c
	$(CC) $(CFLAGS) -c $< -o gt-alloc.o

mmap.o: mmap.c
	$(CC) $(CFLAGS) -c $< -o mmap.o

fileopen.o: fileopen.c
	$(CC) $(CFLAGS) -c $< -o fileopen.o

lzOG.o: lzOG.cpp
	$(CC) $(CFLAGS) -c $< -o lzOG.o

main.o: main.cpp
	$(CC) $(CFLAGS) -c $< -o main.o

test: $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECTS) -o lz_factor_file

clean:
	rm $(OBJECTS)
