/**
 *
 * Author: Marc Khoury <khoury@eecs.berkeley.edu>
 *
 * Copyright (C) 2014 Marc Khoury
 *
 * This version is released under the Eclipse Public License
 * with the Graphviz distribution.
 */

#ifndef VIEWER_H
#define VIEWER_H

#include <stdlib.h>
#include <graphviz/cgraph.h>

#include "linalg.h"

void init_viewer(Agraph_t* graph, int length);
void append_layout(mat m);
void viewer(int argc, char* argv[]);

#endif
