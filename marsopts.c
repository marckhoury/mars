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

#include "marsopts.h"

void usage(int status)
{
    fprintf(stderr, "%s", use_string);
    exit(status);
}

struct marsopts init(int argc, char* argv[])
{
    int c;
    //{cmd,fin,fout,k,power,dim,scale,max_iter,verbose, viewer}
    struct marsopts opts = {argv[0], stdin, stdout, 100, 1, 2, 72, 200, 0, 0};
    while((c = getopt(argc,argv,":o:k:p:d:s:i:cvg")) != -1) {
        switch(c) {
            case 'o':
                opts.fout = fopen(optarg,"w");
                if(!opts.fout) {
                    fprintf(stderr, "Failed to open %s for writing\n", optarg);
                    exit(1);
                }
            break;
            case 'k':
                opts.k = atoi(optarg);
            break;
            case 'p':
                opts.power = atoi(optarg);
            break;
            case 'd':
                opts.dim = atoi(optarg);
            break;
            case 's':
                opts.scale = atoi(optarg);
            break;
            case 'i':
                opts.max_iter = atoi(optarg);
            break;
            case 'c':
                opts.color = 1;
            break;
            case 'v':
                opts.viewer = 1;
            break;
            case 'g':
                opts.given = 1;
            break;
            case '?':
                if(optopt == '?') {
                    usage(0);
                } else {
                    fprintf(stderr, "%s: option -%c unrecognized\n", opts.cmd, optopt);
                    usage(1);
                }
            break;
            case ':':
                fprintf(stderr, "%s: missing argument for option -%c\n", opts.cmd, optopt);
                usage(1);
            break;
        }
    }
    if(optind < argc) {
        opts.fin = fopen(argv[optind],"r");
        if(!opts.fin) {
            fprintf(stderr, "Failed to open %s for reading\n", argv[optind]);
            exit(1);
        }   
    }
    return opts;
}
