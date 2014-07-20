OS = $(shell uname -s)

ifeq ($(OS),Darwin)
    CC = clang
else
    CC = gcc
endif

CFLAGS = -O2
OBJECTS = dijkstra.o graph.o layout.o linalg.o marsopts.o
LIBS = -lcgraph -llapack -lcdt libsfdp.a
TARGET = mars

mars: $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) marsmain.c $(OBJECTS) $(LIBS) $(LDFLAGS)

dijkstra.o: dijkstra.h dijkstra.c
	$(CC) $(CFLAGS) -c dijkstra.c

graph.o: graph.h graph.c
	$(CC) $(CFLAGS) -c graph.c

layout.o: layout.h layout.c
	$(CC) $(CFLAGS) -c layout.c

linalg.o: linalg.h linalg.c
	$(CC) $(CFLAGS) -c linalg.c

marsopts.o: marsopts.h marsopts.c
	$(CC) $(CFLAGS) -c marsopts.c

clean:
	rm -f $(TARGET) $(OBJECTS) *~ 
