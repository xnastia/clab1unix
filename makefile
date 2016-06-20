#
CC=g++

#
CFLAGS=-c -Wall
#
LDFLAGS=


all: prog1

prog1: prog1.o
	$(CC) $(LDFLAGS) prog1.o -o prog1

prog1.o: prog1.c
	$(CC) $(CFLAGS) prog1.c -o prog1.o

clean:
	rm -rf *.o prog1
