#include "dijkstra.h"

static double getlength(Agedge_t * e)
{
    double len;
    char*  lens;
    char*  p;

    if (weight && (*(lens = agxget(e, weight)))) {
        len = strtod (lens, &p);
        if ((len < 0) || (p == lens)) {
            len = 1;
        }
    } else {
        len = 1;
    }
    return len;
}

static int cmpf(Dt_t * d, void *key1, void *key2, Dtdisc_t * disc)
{
    double t;
    t = getdist((Agnode_t *) key1) - getdist((Agnode_t *) key2);
    
    if (t < 0) {
        return -1;
    }
    if (t > 0) {
        return 1;
    }
    if (key1 < key2) {
        return -1;
    }
    if (key1 > key2) {
        return 1;
    }
    return 0;
}

static Dtdisc_t MyDisc = {
    0,              /* int key */
    0,              /* int size */
    -1,             /* int link */
    0,              /* Dtmake_f makef */
    0,              /* Dtfree_f freef */
    cmpf,           /* Dtcompar_f comparf */
    0,              /* Dthash_f hashf */
    0,              /* Dtmemory_f memoryf */
    0               /* Dtevent_f eventf */
};

static Agnode_t *extract_min(Dict_t * Q)
{
    Agnode_t *rv;
    rv = (Agnode_t*) dtfirst(Q);
    dtdelete(Q, rv);
    return rv;
}

static void update(Dict_t * Q, Agnode_t * dest, Agnode_t * src, double len)
{
    double newlen = getdist(src) + len;
    double oldlen = getdist(dest);

    if (oldlen == 0) {      /* first time to see dest */
        setdist(dest, newlen);
        dtinsert(Q, dest);
    } else if (newlen < oldlen) {
        dtdelete(Q, dest);
        setdist(dest, newlen);
        dtinsert(Q, dest);
    }
}

double* dijkstra(Agraph_t* g, Agnode_t* n)
{
    Dict_t *Q;
    Agnode_t *u;
    Agedge_t *e;
    double* dist;
    
    Q = dtopen(&MyDisc, Dtoset);
    dist = (double*) malloc(sizeof(double)*agnnodes(g));
    
    setdist(n, 0);
    dtinsert(Q, n);
    while ((u = extract_min(Q))) {
        setDone (u,1);
        for (e = agfstedge(g, u); e; e = agnxtedge(g, e, u)) {
            if (!isDone(e->node)) update(Q, e->node, u, getlength(e));
        }
    }
    dtclose(Q);
    
    for(u = agfstnode(g);  u; u = agnxtnode(g,u)) {
        dist[getid(u)] = getdist(u);
        setdist(u,0);
        setDone(u,0);
    }
    return dist;
}
