/**
 *
 * Author: Marc Khoury <khoury@eecs.berkeley.edu>
 *
 * Copyright (C) 2014 Marc Khoury
 *
 * This version is released under the Eclipse Public License
 * with the Graphviz distribution.
 *
 * A version is also available on GitHub: https://github.com/marckhoury/mars
 * If you make improvements or bug fixes to this code it would be much
 * appreciated if you could also contribute those changes back to the
 * GitHub repository.
 */

#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <stdio.h>
#include <stdlib.h>
#include <graphviz/cgraph.h>
#include "graph.h"

double* dijkstra(Agraph_t* g, Agnode_t* n);

#endif
