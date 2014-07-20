/* $Id$Revision: */
/* vim:set shiftwidth=4 ts=8: */

/*************************************************************************
 * Copyright (c) 2011 AT&T Intellectual Property 
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors: See CVS logs. Details at http://www.graphviz.org/
 *************************************************************************/

#ifndef GENERAL_H
#define GENERAL_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define real double

#define set_flag(a, flag) ((a)=((a)|(flag)))
#define test_flag(a, flag) ((a)&(flag))
#define clear_flag(a, flag) ((a) &=(~(flag)))

#define MALLOC malloc
#define REALLOC realloc

#define N_NEW(n,t)   (t*)malloc((n)*sizeof(t))
#define NEW(t)       (t*)malloc(sizeof(t))
#define MAX(a,b) ((a)>(b)?(a):b)
#define MIN(a,b) ((a)<(b)?(a):b)
#define ABS(a) (((a)>0)?(a):(-(a)))

#ifdef TRUE
#undef TRUE
#endif
#define TRUE 1

#ifdef FALSE
#undef FALSE
#endif
#define FALSE 0

#define MAXINT 1<<30
#define PI 3.14159

#define POINTS(inch) 72*(inch)

typedef unsigned int boolean;
extern unsigned char Verbose;

#define FREE free
#define MEMCPY memcpy

#ifndef DEBUG
#ifndef NDEBUG
#define NDEBUG /* switch off assert*/
#endif
#endif

#ifdef DEBUG
extern double _statistics[10];
#endif

#define MACHINEACC 1.0e-16
#define SQRT_MACHINEACC 1.0e-8

#define MINDIST 1.e-15

enum {UNMATCHED = -1};


real distance(real *x, int dim, int i, int j);
real distance_cropped(real *x, int dim, int i, int j);

real point_distance(real *p1, real *p2, int dim);
#endif
