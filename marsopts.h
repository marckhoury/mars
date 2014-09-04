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

#ifndef CMDLINE_H
#define CMDLINE_H

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

static char* use_string = "Usage: mars [-k k] [-p power] [-d dim] [-s scale] [-i iter] [-o outfile] [-cvg?]\n\
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
