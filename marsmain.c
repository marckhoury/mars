/**
 *
 * Author: Marc Khoury <khoury@eecs.berkeley.edu>
 *
 * Copyright (C) 2014 Marc Khoury
 *
 * This version is released under the Eclipse Public License
 * with the Graphviz distribution.
 */

#include <stdio.h>
#include <stdlib.h>
#include <graphviz/cgraph.h>
#include "linalg.h"
#include "marsopts.h"
#include "graph.h"
#include "dijkstra.h"
#include "layout.h"
#include "marsviewer.h"

char* pos_to_str(double* pos, int dim)
{
	char buffer[256];
	char* s;
	int i,length;
	if(dim == 2) {
		length = sprintf(buffer, "%lf, %lf", pos[0], pos[1]);
	} else {
		length = sprintf(buffer, "%lf, %lf, %lf", pos[0], pos[1], pos[2]);
	}
	s = (char*) malloc(sizeof(char)*(length+1));
	for(i = 0; i < length; i++) {
		s[i] = buffer[i];
	}
	s[length] = '\0';
	return s;
}

int main(int argc, char* argv[])
{
	struct marsopts opts = init(argc, argv);
	Agraph_t* g = agread(opts.fin, (Agdisc_t*)NULL);
	Agnode_t* n;
	mat z;
	
	init_graph(g);
	z = mars(g, opts);
	mat_scalar_mult(z, opts.scale);

    if(opts.viewer) {
        viewer(argc, argv);
    } else {
	    for(n = agfstnode(g); n; n = agnxtnode(g,n)) {
		    int id = getid(n);
		    char* s = pos_to_str(&z->m[mindex(id, 0, z)], z->c);
		    agset(n,"pos",s);
		    free(s);
	    }
	
	    agwrite(g, opts.fout);
	}
	mat_free(z);
	clean_up(g);
	agclose(g);
	return 0;
}
