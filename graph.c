#include "graph.h"

static Agnode_t** nodes;

void init_graph(Agraph_t* g)
{
    int i = 0;
    Agnode_t* n;
    
    weight = agattr(g, AGEDGE, "weight", "1.0");
    pos = agattr(g, AGNODE, "pos", "");
    color = agattr(g, AGNODE, "color", "black");
    comment = agattr(g, AGRAPH, "cmd", "");
    
    nodes = (Agnode_t**) malloc(sizeof(Agnode_t*)*agnnodes(g));
    for(n = agfstnode(g); n; n = agnxtnode(g,n)) {
        agbindrec(n, (char*)"nodedata", sizeof(nodedata_t),1);
        setid(n,i);
        nodes[i] = n;
        i++;
    }
}

Agnode_t* get_node(int id)
{
    return nodes[id];
}

void clean_up(Agraph_t* g)
{
    agclean(g, AGNODE, "nodedata");
    free(nodes);
}
