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

#include "general.h"
#include <errno.h>

#ifdef DEBUG
double _statistics[10];
#endif

real distance_cropped(real *x, int dim, int i, int j){
  int k;
  real dist = 0.;
  for (k = 0; k < dim; k++) dist += (x[i*dim+k] - x[j*dim + k])*(x[i*dim+k] - x[j*dim + k]);
  dist = sqrt(dist);
  return MAX(dist, MINDIST);
}

real distance(real *x, int dim, int i, int j){
  int k;
  real dist = 0.;
  for (k = 0; k < dim; k++) dist += (x[i*dim+k] - x[j*dim + k])*(x[i*dim+k] - x[j*dim + k]);
  dist = sqrt(dist);
  return dist;
}

real point_distance(real *p1, real *p2, int dim){
  int i;
  real dist;
  dist = 0;
  for (i = 0; i < dim; i++) dist += (p1[i] - p2[i])*(p1[i] - p2[i]);
  return sqrt(dist);
}
