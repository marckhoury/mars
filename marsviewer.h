#ifndef VIEWER_H
#define VIEWER_H

#include <stdlib.h>
#include <graphviz/cgraph.h>

#include "linalg.h"

void init_viewer(Agraph_t* graph, int length);
void append_layout(mat m);
void viewer(int argc, char* argv[]);

#endif
