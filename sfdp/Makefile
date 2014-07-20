OS = $(shell uname -s)

ifeq ($(OS),Darwin)
    CC = clang
else
    CC = gcc
endif

CFLAGS = -O2
OBJECTS = QuadTree.o LinkedList.o general.o
TARGET = libsfdp.a

libsfdp.a: $(OBJECTS)
	ar cr $(TARGET) $(OBJECTS)

QuadTree.o: QuadTree.h QuadTree.c
	$(CC) $(CFLAGS) -c QuadTree.c

LinkedList.o: LinkedList.h LinkedList.c
	$(CC) $(CFLAGS) -c LinkedList.c

general.o: general.h general.c
	$(CC) $(CFLAGS) -c general.c

clean:
	rm -f $(TARGET) $(OBJECTS) *~ 
