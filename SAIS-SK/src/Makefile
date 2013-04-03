CC = g++
CFLAGS=-Werror -Wall -O3
OBJ=sk-sain.o\
    mmap.o\
    fileopen.o\
    sk-sain-main.o\
    sufcheck.o\
    gt-alloc.o

all:sk-sain.x

sk-sain.x:${OBJ}
	${CC} -o $@ ${OBJ}	

clean:
	rm -f *.o *.x
