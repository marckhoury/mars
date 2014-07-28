#ifndef CMDLINE_H
#define CMDLINE_H
/*
 * Written by Marc Khoury
 */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

struct marsopts {
    char* cmd;
    FILE* fin;
    FILE* fout;
    int k;
    int power;
    int dim;
    int scale;
    int max_iter;
    int color;
    int viewer;
    int given;
};

static char* use_string = "Usage: mars [-k k] [-p power] [-d dim] [-s scale] [-i iter] [-o outfile] [-cg?]\n\
  -k k       - sample k columns from the full laplacian matrix\n\
  -p power   - exponent of the weight matrix\n\
  -d dim     - dimension of the layout\n\
  -s scale   - scale the layout by a constant\n\
  -i iter    - set the maximum iterations to converge\n\
  -o outfile - write output graph to file\n\
  -c         - color anchor nodes\n\
  -v         - visualize layout in interactive viewer\n\
  -g         - use given initial layout, else random\n\
  -?         - print help message\n";

void usage(int status);
struct marsopts init(int argc, char* argv[]);
#endif
