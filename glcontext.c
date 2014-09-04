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

#include "glcontext.h"

struct Context* context_new()
{
    struct Context* ctx = (struct Context*) malloc(sizeof(struct Context));
    ctx->fovy = 60;
    ctx->width = 600, ctx->height = 900;
    ctx->mouse_button[0] = 0, ctx->mouse_button[1] = 0, ctx->mouse_button[2] = 0;
    ctx->animation_mode = 0;
    ctx->frame_rate = 60;
    ctx->znear = 0.001, ctx->zfar = 100;
    ctx->light_pos[0] = 0, ctx->light_pos[1] = 0, ctx->light_pos[2] = -1, ctx->light_pos[3] = 1;
    ctx->ambient[0] = 0, ctx->ambient[1] = 0, ctx->ambient[2] = 0, ctx->ambient[3] = 1;
    ctx->diffuse[0] = 1, ctx->diffuse[1] = 1, ctx->diffuse[3] = 1, ctx->diffuse[3] = 1;
    ctx->specular[0] = 1, ctx->specular[1] = 1, ctx->specular[2] = 1, ctx->specular[3] = 1;
    ctx->global_ambient[0] = 0.4, ctx->global_ambient[1] = 0.4, ctx->global_ambient[2] = 0.4, ctx->global_ambient[3] = 1;
    int i;
    for(i = 0; i < 16; i++) {
        int entry = i == 0 || i == 5 || i == 10 || i == 15;
        ctx->mat[i] = entry;
        ctx->matinv[i] = entry;
    }
    return ctx;
}

void context_free(struct Context* ctx)
{
    free(ctx);
}
