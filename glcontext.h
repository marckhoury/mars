/**
 *
 * Author: Marc Khoury <khoury@eecs.berkeley.edu>
 *
 * Copyright (C) 2014 Marc Khoury
 *
 * This version is released under the Eclipse Public License
 * with the Graphviz distribution.
 */

#ifndef GL_CONTEXT_H
#define GL_CONTEXT_H

#include <stdlib.h>

struct Context
{
    int fovy;
    int width, height;
    int mouse_pos[2];
    int mouse_button[3];
    int animation_mode;
    int frame_rate;
    double drag[3];
    double window[4];
    double znear, zfar;
    float mat[16]; 
    float matinv[16];
    float light_pos[4];
    float ambient[4];
    float diffuse[4];
    float specular[4];
    float global_ambient[4];

};

struct Context* context_new();
void context_free(struct Context* ctx);

#endif
