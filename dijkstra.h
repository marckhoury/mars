#ifndef DIJKSTRA_H
#define DIJKSTRA_H

/*
 * Written by Marc Khoury
 */

#include <stdio.h>
#include <stdlib.h>
#include <graphviz/cgraph.h>
#include "graph.h"

double* dijkstra(Agraph_t* g, Agnode_t* n);

#endif
