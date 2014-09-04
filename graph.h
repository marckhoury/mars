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

#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <graphviz/cgraph.h>

typedef struct
{
    Agrec_t h;
    int id;
    int done;
    double dist;
} nodedata_t; 

Agsym_t* weight;
Agsym_t* pos;
Agsym_t* color; 
Agsym_t* comment;

#define getdist(n) (((nodedata_t*)((n)->base.data))->dist)
#define setdist(n,d) (((nodedata_t*)((n)->base.data))->dist = (d))
#define isDone(n) (((nodedata_t*)((n)->base.data))->done)
#define setDone(n,d) (((nodedata_t*)((n)->base.data))->done = (d))
#define getid(n)  (((nodedata_t*)(n->base.data))->id)
#define setid(n,d) (((nodedata_t*)(n->base.data))->id = (d))

void init_graph(Agraph_t* g);
Agnode_t* get_node(int id);
void clean_up(Agraph_t* g);

#endif
