OS = $(shell uname -s)

ifeq ($(OS),Darwin)
    CC = clang
    LIBS = -lcgraph -llapack -lcdt libsfdp.a -framework GLUT -framework OpenGL -framework Cocoa
else
    CC = gcc
    LIBS = -lcgraph -llapack -lcdt libsfdp.a -lGL -lglut -lGLU
endif

CFLAGS = -g -O2
OBJECTS = dijkstra.o graph.o layout.o linalg.o marsopts.o glcontext.o marsviewer.o
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

glcontext.o: glcontext.h glcontext.c
	$(CC) $(CFLAGS) -c glcontext.c

marsviewer.o: marsviewer.h marsviewer.c
	$(CC) $(CFLAGS) -c marsviewer.c

clean:
	rm -f $(TARGET) $(OBJECTS) *~ 
