/**
 *
 * Author: Marc Khoury <khoury@eecs.berkeley.edu>
 *
 * Copyright (C) 2014 Marc Khoury
 *
 * This version is released under the Eclipse Public License
 * with the Graphviz distribution.
 */

#include "marsviewer.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "glcontext.h"
#include "graph.h"

#define MAX(X,Y) (((X) < (Y)) ? (Y) : (X))
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))

Agraph_t* g;
mat* layouts;
int layouts_size;
int layouts_max;
int curr_layout;
struct Context* ctx;

void init_viewer(Agraph_t* graph, int length)
{
    g = graph;
    layouts_max = length;
    curr_layout = 0;
    layouts_size = 0;
    layouts = (mat*) malloc(layouts_max*sizeof(mat));
}

void scale_layout(int i)
{
    double maxp[3], minp[3];
    mat layout = layouts[i];
    int dim = layout->c;
    int j;

    if(dim == 2) {
        double x = layout->m[mindex(0, 0, layout)];
        double y = layout->m[mindex(0, 1, layout)];
        double z = 0;
        maxp[0] = minp[0] = x;
        maxp[1] = minp[1] = y;
        maxp[2] = minp[2] = z;
    } else if(dim >= 3) {
        double x = layout->m[mindex(0, 0, layout)];
        double y = layout->m[mindex(0, 1, layout)];
        double z = layout->m[mindex(0, 2, layout)];
        maxp[0] = minp[0] = x;
        maxp[1] = minp[1] = y;
        maxp[2] = minp[2] = z;
    }

    for(j = 1; j < layout->r; j++) {
        if(dim == 2) {
            double x = layout->m[mindex(j, 0, layout)];
            double y = layout->m[mindex(j, 1, layout)];
            maxp[0] = MAX(maxp[0], x);
            maxp[1] = MAX(maxp[1], y);
            minp[0] = MIN(minp[0], x);
            minp[1] = MIN(minp[1], y);
        } else if(dim >= 3) { 
            double x = layout->m[mindex(j, 0, layout)];
            double y = layout->m[mindex(j, 1, layout)];
            double z = layout->m[mindex(j, 2, layout)];
            maxp[0] = MAX(maxp[0], x);
            maxp[1] = MAX(maxp[1], y);
            maxp[2] = MAX(maxp[2], z);
            minp[0] = MIN(minp[0], x);
            minp[1] = MIN(minp[1], y);
            minp[2] = MIN(minp[2], z);
        } 
    }
    double center[3] = {(maxp[0] + minp[0])/2, (maxp[1] + minp[1])/2, (maxp[2] + minp[2])/2};
    double scale;
    if(dim == 2) {
        scale = 2.0/MAX(maxp[0] - minp[0], maxp[1] - minp[1]);
    } else if (dim >= 3) {
        scale = 2.0/MAX(maxp[0] - minp[0], MAX(maxp[1] - minp[1], maxp[2] - minp[2]));
    }
    for(j = 0; j < layout->r; j++) {
        if(dim == 2) {
            layout->m[mindex(j, 0, layout)] -= center[0];
            layout->m[mindex(j, 1, layout)] -= center[1];
            layout->m[mindex(j, 0, layout)] *= scale;
            layout->m[mindex(j, 1, layout)] *= scale;
        } else if(dim >= 3) { 
            layout->m[mindex(j, 0, layout)] -= center[0];
            layout->m[mindex(j, 1, layout)] -= center[1];
            layout->m[mindex(j, 2, layout)] -= center[2];
            layout->m[mindex(j, 0, layout)] *= scale;
            layout->m[mindex(j, 1, layout)] *= scale;
            layout->m[mindex(j, 2, layout)] *= scale;
        } 
    }
}

void append_layout(mat m)
{
    if(layouts_size < layouts_max) {
        int i;
        layouts[layouts_size] = mat_new(m->r, m->c);
        for(i = 0; i < m->r*m->c; i++) {
            layouts[layouts_size]->m[i] = m->m[i];
        }
        scale_layout(layouts_size);
        layouts_size++;
    }
}

void free_layouts()
{
    int i;
    for(i = 0; i < layouts_size; i++) {
        mat_free(layouts[i]);
    }
    free(layouts);
}

void inverse(float* m, float* res)
{
    float det;

    float d12 = m[2]*m[7] - m[3]*m[6]; 
    float d13 = m[2]*m[11] - m[3]*m[10];
    float d23 = m[6]*m[11] - m[7]*m[10];
    float d24 = m[6]*m[15] - m[7]*m[2];
    float d34 = m[10]*m[15] - m[11]*m[14];
    float d41 = m[14]*m[3] - m[15]*m[2];

    res[0] =  m[5]*d34 - m[9]*d24 + m[13]*d13;
    res[1] = -m[1]*d34 - m[9]*d41 - m[13]*d13;
    res[2] =  m[1]*d24 + m[5]*d41 + m[13]*d12;
    res[3] = -m[1]*d23 + m[5]*d13 - m[9]*d12;

    det = m[0]*res[0] + m[4]*res[1] + m[8]*res[2] + m[12]*res[3];

    if(det != 0.0) {
        float invdet = 1.0/det;

        res[0] *= invdet;
        res[1] *= invdet;
        res[2] *= invdet;
        res[3] *= invdet;
        res[4] = -(m[4]*d34 - m[8]*d24 + m[12]*d23) * invdet; 
        res[5] =  (m[0]*d34 + m[8]*d41 + m[12]*d13) * invdet;
        res[6] = -(m[0]*d24 + m[4]*d41 + m[12]*d12) * invdet;
        res[7] =  (m[0]*d23 - m[4]*d13 + m[8]*d12) * invdet;
        
        d12 = m[0]*m[5] - m[1]*m[4];
        d13 = m[0]*m[9] - m[1]*m[8];
        d23 = m[4]*m[8] - m[5]*m[8];
        d24 = m[4]*m[13] - m[5]*m[14];
        d34 = m[8]*m[13] - m[9]*m[12];
        d41 = m[12]*m[1] - m[13]*m[0];

        res[9] = -(m[3]*d34 + m[11]*d41 + m[15]*d13) * invdet;
        res[10] = (m[3]*d24 + m[7]*d41 + m[15]*d12) * invdet;
        res[11] = -(m[3]*d23 - m[7]*d13 + m[11]*d12) * invdet;
        res[12] = -(m[6]*d34 - m[10]*d24 + m[14]*d23) * invdet;
        res[13] = (m[2]*d34 + m[10]*d41 + m[14]*d13) * invdet;
        res[14] = -(m[2]*d24 + m[6]*d41 + m[14]*d12) * invdet;
        res[15] = (m[2]*d23 - m[6]*d13 + m[10]*d12) * invdet;
    }
}

void display()
{
    int i;
    Agnode_t* v;
    Agedge_t* e;
    mat layout = layouts[curr_layout];    
    int dim = layout->c;
        
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(0, 0, -4);
    glMultMatrixf(ctx->mat);

    glPointSize(5);
    glBegin(GL_POINTS);
    glColor3f(1, 0, 0);
    for(i = 0; i < layout->r; i++) {
        if(dim == 2) {
            double x = layout->m[mindex(i, 0, layout)];
            double y = layout->m[mindex(i, 1, layout)];
            glVertex2d(x, y);
        } else if(dim >= 3) { 
            double x = layout->m[mindex(i, 0, layout)];
            double y = layout->m[mindex(i, 1, layout)];
            double z = layout->m[mindex(i, 2, layout)];
            glVertex3d(x, y, z);
        }    
    }
    glEnd();
    
    glBegin(GL_LINES);
    glColor3f(0, 0, 0);
    for(v = agfstnode(g);  v; v = agnxtnode(g, v)) {
        for(e = agfstedge(g, v); e; e = agnxtedge(g, e, v)) {
            Agnode_t* u = aghead(e);
            Agnode_t* w = agtail(e);
            int uid = getid(u), wid = getid(w);
            if(dim == 2) {
                double x1 = layout->m[mindex(uid, 0, layout)];
                double y1 = layout->m[mindex(uid, 1, layout)];
                double x2 = layout->m[mindex(wid, 0, layout)];
                double y2 = layout->m[mindex(wid, 1, layout)];
                glVertex2d(x1, y1);
                glVertex2d(x2, y2);
            } else if(dim >= 3) {
                double x1 = layout->m[mindex(uid, 0, layout)];
                double y1 = layout->m[mindex(uid, 1, layout)];
                double z1 = layout->m[mindex(uid, 2, layout)];
                double x2 = layout->m[mindex(wid, 0, layout)];
                double y2 = layout->m[mindex(wid, 1, layout)];
                double z2 = layout->m[mindex(wid, 2, layout)];
                glVertex3d(x1, y1, z1);
                glVertex3d(x2, y2, z2);
            }
        }    
    }
    glEnd();

    glPopMatrix();
    glutSwapBuffers();
}

void timer(int k)
{
    if(ctx->animation_mode) {
        if(curr_layout < layouts_size - 1) {
            curr_layout++;
        } else {
            ctx->animation_mode = 0;
        }
        glutPostRedisplay();
    }
    glutTimerFunc(1000/ctx->frame_rate, timer, k);
}

void reshape(int w, int h)
{
    ctx->window[0] = -((double) w)/h;
    ctx->window[1] = -ctx->window[0];
    ctx->window[2] = 1.0;
    ctx->window[3] = -1.0;

    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(ctx->fovy, ((double) w)/h, ctx->znear, ctx->zfar);
    glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y)
{
    switch(key) {
        case 'f':
            curr_layout = 0;
        break;
        case 'l':
            curr_layout = layouts_size - 1;
        break;
        case 'd':
            if(curr_layout < layouts_size - 1) {
                curr_layout++;
            }
        break;
        case 's':
            if(curr_layout > 0) {
                curr_layout--;
            }
        break;
        case 'a':
            ctx->animation_mode = !ctx->animation_mode;
        break;
        case 27: //ESC key
            exit(0);
        break;
        default:
            fprintf(stderr, "Unassigned character: %c\n", key);
        break;
    }
}

void world_coords(int x, int y, double* p)
{
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    p[0] = ((double)(x-viewport[0]))/viewport[2];
    p[1] = ((double)(y-viewport[1]))/viewport[3];

    p[0] = ctx->window[0] + p[0]*(ctx->window[1]-ctx->window[0]);
    p[1] = ctx->window[2] + p[1]*(ctx->window[3]-ctx->window[2]);
    p[2] = ctx->znear;
}

void mouse(int button, int state, int x, int y)
{
    int cursor = GLUT_CURSOR_RIGHT_ARROW;
    if(state == GLUT_DOWN) {
        if(button == GLUT_LEFT_BUTTON) {
            cursor = GLUT_CURSOR_CYCLE;
            ctx->mouse_button[0] = 1;
        } else if(button == GLUT_MIDDLE_BUTTON) {
            cursor = GLUT_CURSOR_CROSSHAIR;
            ctx->mouse_button[1] = 1;
        } else if(button == GLUT_RIGHT_BUTTON) {
            cursor = GLUT_CURSOR_UP_DOWN;
            ctx->mouse_button[2] = 1;
        }
    } else {
        ctx->mouse_button[0] = ctx->mouse_button[1] = ctx->mouse_button[2] = 0;
    }
    glutSetCursor(cursor);
    ctx->mouse_pos[0] = x;
    ctx->mouse_pos[1] = y;
    world_coords(x, y, ctx->drag);
}

void motion(int x, int y)
{
    int dx = x - ctx->mouse_pos[0];
    int dy = y - ctx->mouse_pos[1];
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    if(dx == 0 && dy == 0) {
        return;
    } else if(ctx->mouse_button[0]) {
        double angle = sqrt(dy*dy + dx*dx)/((double)(viewport[2]+1))*180.0;
        double rx = ctx->matinv[0]*dy + ctx->matinv[4]*dx;
        double ry = ctx->matinv[1]*dy + ctx->matinv[5]*dx;
        double rz = ctx->matinv[2]*dy + ctx->matinv[6]*dx;
        glRotatef(angle, rx, ry, rz);
    } else if(ctx->mouse_button[1]) {
        double p[3];
        world_coords(x, y, p);
        glLoadIdentity();
        glTranslatef(p[0] - ctx->drag[0], p[1] - ctx->drag[1], p[2] - ctx->drag[2]);
        glMultMatrixf(ctx->mat);
        ctx->drag[0] = p[0], ctx->drag[1] = p[1], ctx->drag[2] = p[2];
    } else if(ctx->mouse_button[2]) {
        glLoadIdentity();
        glTranslatef(0, 0, dy*0.01);
        glMultMatrixf(ctx->mat);
    }

    ctx->mouse_pos[0] = x;
    ctx->mouse_pos[1] = y;
    
    glGetFloatv(GL_MODELVIEW_MATRIX, ctx->mat);
    inverse(ctx->mat, ctx->matinv);
    glutPostRedisplay();
}

void idle()
{
    glutPostRedisplay();
}

void init_opengl()
{
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(ctx->fovy, 1, ctx->znear, ctx->zfar);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void menu()
{
    printf("mars viewer menu\n");
    printf("a - Animation Mode\n");
    printf("d - Forward one frame\n");
    printf("s - Backward one frame\n");
    printf("f - First frame\n");
    printf("l - Last frame\n");
    printf("Esc - Quit\n");
}

void viewer(int argc, char* argv[])
{
    ctx = context_new();
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(ctx->width, ctx->height);
    glutCreateWindow("mars viewer");
    init_opengl();
    menu();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutIdleFunc(idle);
    glutTimerFunc(0, timer, 0);

    glutMainLoop();

    context_free(ctx);
    free_layouts();
}
