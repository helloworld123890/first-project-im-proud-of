CC = gcc
CFLAGS = -Wall -std=c99 -pedantic
OBJS = mol.o
SRC = mol.c


all: libmol.so mol.so

clean:  
	rm -f *.o *.so

libmol.so: mol.o
	$(CC) mol.o -shared -o mol.so

mol.o:  mol.c mol.h
	$(CC) $(CFLAGS) -c mol.c -fPIC -o mol.o
