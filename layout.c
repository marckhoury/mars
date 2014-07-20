#include "layout.h"

/*Lapack declarations */
void dgesvd_(const char* jobu, const char* jobvt, const int* M, const int* N,
             double* A, const int* lda, double* S, double* U, const int* ldu,
             double* VT, const int* ldvt, double* work,const int* lwork, const
             int* info);

void dgeqrf_(const int* M, const int* N, double* A, const int* lda,
             double* TAU, double* work, const int* lwork, int* info);

void dorgqr_(const int* M, const int* N, const int* K, double* A, const int* lda,
             double* TAU, double* work, const int* lwork, int* info);

void dtrtrs_(const char* UPLO, const char* TRANS, const char* DIAG, const int* N, 
             const int* NRHS, double* A, const int* lda, double* B, const int* ldb, 
             int* info);
/*----------------------*/

static void iarrset(int* arr, int n, int val)
{
    int i;
    for(i = 0; i < n; i++) {
        arr[i] = val;   
    }
}

static void darrset(double* arr, int n, double val)
{
    int i;
    for(i = 0; i < n; i++) {
        arr[i] = val;   
    }
}

static mat some_pairs_dist(mat z, int* anchors, int k)
{
    int i, j, x, dim = z->c;
    mat d = mat_new(k,z->r);
    for(i = 0; i < k; i++) {
        int anchor = anchors[i];
        for(j = 0; j < z->r; j++) {
            double dist = 0;
            for(x = 0; x < dim; x++) {
                dist += ((z->m[mindex(anchor,x,z)]-z->m[mindex(j,x,z)])
                        * (z->m[mindex(anchor,x,z)]-z->m[mindex(j,x,z)]));
            }
            d->m[mindex(i,j,d)] = sqrt(dist);
        }
    }
    return d;   
}

static double stress(mat z, mat dij, int* anchors, int power)
{
    int i,j;
    double sum;
    mat d = some_pairs_dist(z,anchors,dij->r);
    mat_sub(d,dij);
    for(i = 0; i < d->r; i++) {
        for(j = 0; j < d->c; j++) {
            d->m[mindex(i,j,d)] = pow(d->m[mindex(i,j,d)],2);
            if(dij->m[mindex(i,j,d)] != 0) {
                d->m[mindex(i,j,d)] *= 1.0 / pow(dij->m[mindex(i,j,d)], power);
            }   
        }
    }
    sum = mat_accu(d);
    mat_free(d);
    return sum;
}

static double* prob(mat x, int power)
{
    int i,j;
    double* p = (double*) malloc(sizeof(double)*x->r);
    double frob_norm = 0, sum;
    for(i = 0; i < x->r; i++) {
        sum = 0;
        for(j = 0; j < x->c; j++) {
            if(x->m[mindex(i,j,x)] != 0) {
                sum += (1.0 / pow(x->m[mindex(i,j,x)], 2*power));
            }
        }
        p[i] = sum;
        frob_norm += sum;
    }
    frob_norm *= (((double)x->c)/x->r);
    for(i = 0; i < x->r; i++) {
        p[i] /= frob_norm;  
    }
    return p;
}

static mat sample_cols(Agraph_t* g, mat x, int power)
{
    mat m = mat_new(x->c, x->r);
    int i, j, k = x->r;
    double* p = prob(x,power);
    for(i = 0; i < x->r; i++) {
        for(j = 0; j < x->c; j++) {
            if(x->m[mindex(i,j,x)] == 0) {
                m->m[mindex(j,i,m)] = 0;
            } else {
                m->m[mindex(j,i,m)] = (-1.0/pow(x->m[mindex(i,j,x)],power))/sqrt(k*p[i]);
            }
        }
    }
    free(p);
    return m;
}

static void singular_vectors(Agraph_t* g, mat x, int power, mat u, double* s)
{
    int i, j, loc, n = x->c, k = x->r;
    int lda = n, ldu = n, ldvt = k, info, lwork = -1;
    mat c = sample_cols(g,x,power);
    double* c_lapack = (double*) malloc(sizeof(double)*n*k);
    double* s_lapack = (double*) malloc(sizeof(double)*k);
    double* u_lapack = (double*) malloc(sizeof(double)*n*k);
    double* vt_lapack = (double*) malloc(sizeof(double)*k*k);
    double* work;
    double wkopt;
    
    loc = 0;
    for(j = 0; j < c->c; j++) {//Write c into array in column major order
        for(i = 0; i < c->r; i++) {
            c_lapack[loc++] = c->m[mindex(i,j,c)];
        }
    }
    //Query for optimal size of work array
    dgesvd_("S","S", &n, &k, c_lapack, &lda, s_lapack, u_lapack, &ldu, vt_lapack, &ldvt, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work = (double*) malloc(sizeof(double)*lwork);
    //Compute svd
    dgesvd_("S","S", &n, &k, c_lapack, &lda, s_lapack, u_lapack, &ldu, vt_lapack, &ldvt, work, &lwork, &info);
    
    for(i = 0; i < n; i++) {
        for(j = 0; j < k; j++) {
            u->m[mindex(i,j,u)] = u_lapack[i+j*ldu];        
        }
    }
    for(i = 0; i < k; i++) {
        s[i] = s_lapack[i]; 
    }
    mat_free(c);
    free(c_lapack);
    free(s_lapack);
    free(u_lapack);
    free(vt_lapack);
    free(work);
}

static void upper_tri_solve(mat u, double* b)
{
    double* a = (double*) malloc(sizeof(double)*u->r*u->c);
    char uplo = 'U';
    char trans = 'N';
    char diag = 'N';
    int nrhs = 1, lda = u->r, ldb = u->c;
    int info;
    int loc, i, j;
    loc = 0;
    for(j = 0; j < u->c; j++) {
        for(i = 0; i < u->r; i++) {
            a[loc++] = u->m[mindex(i,j,u)];
        }
    }
    dtrtrs_(&uplo, &trans, &diag, &u->c, &nrhs, a, &lda, b, &ldb, &info);

    free(a);
}

static void qr_factorize(mat m, mat q, mat r)
{
    int i,j;
    double* a = (double*) malloc(sizeof(double)*m->r*m->c);
    double* tau = (double*) malloc(sizeof(double)*m->c);
    double* work;
    double workopt;
    int lda = m->r;
    int lwork = -1;
    int info;
    int loc;
    
    loc = 0;
    for(j = 0; j < m->c; j++) {
        for(i = 0; i < m->r; i++) {
            a[loc++] = m->m[mindex(i,j,m)];     
        }   
    }
    /*QR decomposition into a compact format*/
    dgeqrf_(&m->r, &m->c, a, &lda, tau, &workopt, &lwork, &info); 
    lwork = (int) workopt;
    work = (double*) malloc(sizeof(double)*lwork);
    dgeqrf_(&m->r, &m->c, a, &lda, tau, work, &lwork, &info);
    
    for(i = 0; i < m->c; i++) {
        for(j = i; j < m->c; j++) {
            r->m[mindex(i,j,r)] = a[i+j*lda];
            a[i+j*lda] = (i == j) ? 1 : 0;  
        }   
    }
    
    lwork = -1;
    dorgqr_(&m->r, &m->c, &m->c, a, &lda, tau, &workopt, &lwork, &info);
    lwork = (int) workopt;
    work = (double*) realloc(work, sizeof(double)*lwork);
    dorgqr_(&m->r, &m->c, &m->c, a, &lda, tau, work, &lwork, &info);

    for(i = 0; i < q->r; i++) {
        for(j = 0; j < q->c; j++) {
            q->m[mindex(i,j,q)] = a[i+j*lda];
        }   
    }   
    
    free(work);
    free(a);
    free(tau);
}

static void select_anchors(Agraph_t* g, mat x, int* anchors, int k)
{
    int i, j, max_dist, max_index;
    Agnode_t* n;
    double *d, *dmin;
    
    n = agfstnode(g);
    d = dijkstra(g,n);
    dmin = (double*) malloc(sizeof(double)*agnnodes(g));
    darrset(dmin, agnnodes(g), DBL_MAX);
    
    max_dist = d[0];
    max_index = 0;
    for(i = 0; i < agnnodes(g); i++) {
        if(max_dist < d[i]) {
            max_dist = d[i];
            max_index = i;
        }
    }
    for(i = 0; i < k; i++) {
        n = get_node(max_index);
        anchors[i] = max_index;
        free(d);
        d = dijkstra(g,n);
        for(j = 0; j < agnnodes(g); j++) {
            x->m[mindex(i,j,x)] = d[j];
            dmin[j] = MIN(dmin[j],d[j]);        
        }
        max_dist = dmin[0];
        max_index = 0;
        for(j = 0; j < agnnodes(g); j++) {
            if(max_dist < dmin[j]) {
                max_dist = dmin[j];
                max_index = j;
            }
        }
    }
    free(d);
    free(dmin);
}

static int* graph_cluster(Agraph_t* g, mat dij, int* anchors)
{
    int k = dij->r;
    int* clusters = (int*) malloc(sizeof(int)*dij->c);
    int i, j;
    for(i = 0; i < dij->c; i++) {
        int min_index = 0, min_dist = dij->m[mindex(0,i,dij)];
        for(j = 1; j < k; j++) {
            if(min_dist > dij->m[mindex(j,i,dij)]) {
                min_dist = dij->m[mindex(j,i,dij)];
                min_index = j;          
            }
        }
        clusters[i] = min_index;
    }
    return clusters;
}

static mat barnes_hut_cluster(mat z, mat dij, int* clusters, int power)
{
    mat lz = mat_new(z->r, z->c);
    int k = dij->r;
    int n = z->r;
    int i, j, qt_i, node;
    real** x = (real**) malloc(sizeof(real*)*k);
    real** wgts = (real**) malloc(sizeof(real*)*k);
    real* x_all = (real*) malloc(sizeof(real)*z->r*z->c);
    real *center = NULL, *supernode_wgts = NULL, *distances = NULL;
    real bh = 0.65, counts = 0;
    int flag, nsuper, nsupermax, dim = z->c, max_qtree_level = 10;
    QuadTree* qt = (QuadTree*) malloc(sizeof(QuadTree)*2*k);
    int* sizes = (int*) malloc(sizeof(int)*k);
    int* curr_loc = (int*) malloc(sizeof(int)*k);
    iarrset(sizes,k,0);
    iarrset(curr_loc,k,0);
    
    for(i = 0; i < n; i++) {
        sizes[clusters[i]]++;   
    }
    for(i = 0; i < k; i++) {
        x[i] = (real*) malloc(sizeof(real)*sizes[i]*z->c);
        wgts[i] = (real*) malloc(sizeof(real)*sizes[i]);
    }
    for(i = 0; i < z->r; i++) {
        for(j = 0; j < z->c; j++) {
            x_all[i*z->c+j] = z->m[mindex(i,j,z)];      
        }   
    }
    for(i = 0; i < n; i++) {
        int cluster = clusters[i];
        for(j = 0; j < z->c; j++) {
            x[cluster][curr_loc[cluster]*z->c+j] = z->m[mindex(i,j,z)];
        }
        curr_loc[cluster]++;
    }
    iarrset(curr_loc,k,0);
    for(i = 0; i < n; i++) {
        int cluster = clusters[i];
        wgts[cluster][curr_loc[cluster]] = dij->m[mindex(cluster,i,dij)] < EPSILON 
                                           ? 0 : 1.0/pow(dij->m[mindex(cluster,i,dij)],power-1);
        curr_loc[cluster]++;
    }
    for(i = 0; i < k; i++) {
        qt[i] = QuadTree_new_from_point_list(dim, sizes[i], max_qtree_level, x[i], NULL);
    }
    for(i = 0; i < k; i++) {
        qt[i+k] = QuadTree_new_from_point_list(dim, sizes[i], max_qtree_level, x[i], wgts[i]);
    }
    for(i = 0; i < z->r; i++) {
        for(qt_i = 0; qt_i < k; qt_i++) {
            double d = dij->m[mindex(qt_i,i,dij)];
            if(d < EPSILON) { //Node i is the center of some cluster
                QuadTree_get_supernodes(qt[qt_i+k], bh, &(x_all[dim*i]), -1, &nsuper, &nsupermax, 
                                        &center, &supernode_wgts, &distances, &counts, &flag);
                d = 1;
            } else {
                QuadTree_get_supernodes(qt[qt_i], bh, &(x_all[dim*i]), -1, &nsuper, &nsupermax, 
                                        &center, &supernode_wgts, &distances, &counts, &flag);
                d = 1.0 / pow(d,power-1);
            }
            for(j = 0; j < dim; j++) {
                for(node = 0; node < nsuper; node++) {
                    if(distances[node] > EPSILON && d > EPSILON)
                        lz->m[mindex(i,j,lz)] += d*supernode_wgts[node] 
                                                 * (z->m[mindex(i,j,z)]-center[node*dim+j]) 
                                                 / distances[node];
                }
            }
        }   
    }
    for(i = 0; i < k; i++) {
        free(x[i]);
        free(wgts[i]);  
    }
    for(i = 0; i < 2*k; i++) {
        QuadTree_delete(qt[i]); 
    }
    free(qt);
    free(sizes);
    free(curr_loc);
    free(x_all);
    free(center);
    free(distances);
    free(supernode_wgts);
    return lz;
}

static mat barnes_hut(mat z)
{
    mat lz = mat_new(z->r, z->c);
    real* x = (real*) malloc(sizeof(real)*z->r*z->c);
    real *center = NULL, *supernode_wgts = NULL, *distances = NULL;
    real bh = 0.65, counts = 0; /* Barnes-Hut constant, if width(snode)/dist[i,snode] < bh, treat snode as a supernode.*/   
    int flag, nsuper, nsupermax, dim = z->c, max_qtree_level = 10; //max_qtree_level, like bh, should also be a command line argument
    QuadTree qt = NULL;
    int i, j, node;
    
    for(i = 0; i < z->r; i++) {
        for(j = 0; j < z->c; j++) {
            x[i*z->c + j] = z->m[mindex(i,j,z)];        
        }
    }
    qt = QuadTree_new_from_point_list(dim, z->r, max_qtree_level, x, NULL);
    for(i = 0; i < z->r; i++) {
        QuadTree_get_supernodes(qt, bh, &(x[dim*i]), i, &nsuper, &nsupermax, 
                                &center, &supernode_wgts, &distances, &counts, &flag);
        
        for(j = 0; j < dim; j++) {
            for(node = 0; node < nsuper; node++) {
                if(distances[node] > EPSILON) {
                    lz->m[mindex(i,j,lz)] += supernode_wgts[node]
                                             * (z->m[mindex(i,j,z)]-center[node*dim+j]) 
                                             / distances[node];
                } else {
                    lz->m[mindex(i,j,lz)] = 0;
                }
            }
        }
    }
    free(x);
    free(center);
    free(supernode_wgts);
    free(distances);
    QuadTree_delete(qt);
    return lz;
}

static double* inv_mul_ax(double* a, double* x, int n)
{
    double* right = (double*) malloc(sizeof(double)*n);
    double* result = (double*) malloc(sizeof(double)*n);
    double right_side = 0, inner_inverse = 0, scalar;
    int i;
    for(i = 0; i < n; i++) { //1^TA^{-1}x
        right_side += x[i]/a[i];
        inner_inverse += 1.0/a[i];
    }
    inner_inverse += n;
    inner_inverse = 1.0/inner_inverse; //(I + 1^TA^{-1})^{-1}
    scalar = inner_inverse*right_side;//(I + 1^TA^{-1})^{-1}*1^TA^{-1}x
    for(i = 0; i < n; i++) { //A^{-1}*1*(I + 1^TA^{-1})^{-1}*1^TA^{-1}x
        right[i] = scalar/a[i];
    }
    for(i = 0; i < n; i++) { //A^{-1}x - A^{-1}*1*(I + 1^TA^{-1})^{-1}*1^TA^{-1}x
        result[i] = x[i]/a[i] - right[i];
    }
    free(right);
    return result;
}

static double* inv_mul_full(double* a, double* x, int n, mat u, mat v, mat q_t, mat r)
{
    int i;
    double* a_x = inv_mul_ax(a,x,n); //A^{-1}x
    double* tmp1 = mat_vec_mult(v,a_x); //VA^{-1}x
    double* tmp2 = mat_vec_mult(q_t,tmp1); //Q^TVA^{-1]x
    upper_tri_solve(r,tmp2); //R^{-1}Q^TVA^{-1}x
    free(tmp1);
    tmp1 = mat_vec_mult(u,tmp2); //UR^{-1}Q^TVA^{-1}x
    free(tmp2);
    tmp2 = inv_mul_ax(a,tmp1,n); //A^{-1}UR^{-1}Q^TVA^{-1}x
    for(i = 0; i < n; i++) { //A^{-1}x-A^{-1}UR^{-1}Q^TVA^{-1}x
        a_x[i] = a_x[i] - tmp2[i];  
    }
    free(tmp1);
    free(tmp2);
    return a_x;
}

static void adjust_anchors(Agraph_t* g, int* anchors, int k, mat z)
{
    int i, j;
    double* centroid = (double*) malloc(sizeof(double)*z->c);
    for(i = 0; i < k; i++) {
        Agnode_t* n = get_node(anchors[i]);
        Agedge_t* e;
        int degree = 0;
        for(j = 0; j < z->c; j++) {
            centroid[j] = 0;        
        }
        for(e = agfstedge(g,n); e; e = agnxtedge(g,e,n)) {
            int v = (n == aghead(e)) ? getid(agtail(e)) : getid(aghead(e));
            for(j = 0; j < z->c; j++) {
                centroid[j] += z->m[mindex(v,j,z)];         
            }
            degree++;       
        }
        for(j = 0; j < z->c; j++) {
            if(degree != 0)
                centroid[j] /= degree;
            z->m[mindex(anchors[i],j,z)] = centroid[j];
        }
    }
    free(centroid);
}

static void update_anchors(mat z, mat dij, int* anchors, int power)
{
    int i, j, k;
    double* xi = (double*) malloc(sizeof(double)*z->c);
    double* xi_next = (double*) malloc(sizeof(double)*z->c);
    double* xj = (double*) malloc(sizeof(double)*z->c);
    for(i = 0; i < dij->r; i++) {
        int n = anchors[i];
        double wgt_sum = 0;
        darrset(xi_next, z->c, 0);
        for(j = 0; j < z->c; j++) {
            xi[j] = z->m[mindex(n,j,z)];        
        }
        for(j = 0; j < z->r; j++) {
            if(n != j) {
                double wgt = 1.0/pow(dij->m[mindex(i,j,dij)],power);
                double dist = 0;
                for(k = 0; k < z->c; k++) {
                    xj[k] = z->m[mindex(j,k,z)];                
                }
                for(k = 0; k < z->c; k++) {
                    dist += (xi[k]-xj[k])*(xi[k]-xj[k]);
                }
                dist = sqrt(dist);
                if(dist > EPSILON) {
                    for(k = 0; k < z->c; k++) {
                        xi_next[k] += wgt*(xj[k]+dij->m[mindex(i,j,dij)]*(xi[k]-xj[k])/dist);
                    }
                    wgt_sum += wgt;     
                }
            }
        }
        for(j = 0; j < z->c; j++) {
            z->m[mindex(n,j,z)] = xi_next[j]/wgt_sum;       
        }
    }
    free(xi);
    free(xi_next);
    free(xj);
}

static void translate_by_centroid(mat z)
{
    int i, j;
    double* centroid = (double*) malloc(sizeof(double)*z->c);
    darrset(centroid, z->c, 0);
    
    for(i = 0; i < z->r; i++) {
        for(j = 0; j < z->c; j++) {
            centroid[j] += z->m[mindex(i,j,z)];     
        }
    }
    for(i = 0; i < z->c; i++) {
        centroid[i] /= z->r;
    }
    for(i = 0; i < z->r; i++) {
        for(j = 0; j < z->c; j++) {
            z->m[mindex(i,j,z)] -= centroid[j];     
        }
    }
    free(centroid);
}

static double* mat_mult_for_d(mat u, double* s, mat u_t, double* ones)
{
    double* tmp = (double*) malloc(sizeof(double)*u->c);
    double* d = (double*) malloc(sizeof(double)*u->r);
    int i,j;
    
    for(i = 0; i < u_t->r; i++) { //tmp = u_t*ones
        double sum = 0;
        for(j = 0; j < u_t->c; j++) {
            sum += u_t->m[mindex(i,j,u_t)]*ones[j];     
        }
        tmp[i] = sum;
    }
    for(i = 0; i < u->c; i++) { //tmp = diagmat(s)*tmp
        tmp[i] *= s[i]; 
    }
    
    for(i = 0; i < u->r; i++) { //d = u*tmp
        double sum = 0;
        for(j = 0; j < u->c; j++) {
            sum += u->m[mindex(i,j,u)]*tmp[j];      
        }
        d[i] = sum;
    }
    
    free(tmp);
    return d;
}

static mat get_positions(Agraph_t* g, int dim)
{
    int i;
    mat z = mat_new(agnnodes(g),dim);
    for(i = 0; i < z->r; i++) {
        Agnode_t* n = get_node(i);
        char* pos = agget(n,"pos");
        if(dim == 2) {
            sscanf(pos,"%lf,%lf",&z->m[mindex(i,0,z)],&z->m[mindex(i,1,z)]);
        } else {
            sscanf(pos,"%lf,%lf,%lf",&z->m[mindex(i,0,z)],&z->m[mindex(i,1,z)],&z->m[mindex(i,2,z)]);
        }
    }
    return z;
}

mat mars(Agraph_t* g, struct marsopts opts)
{
    int i, j, n = agnnodes(g), k = MIN(n, MAX(opts.k, 2)), iter = 0;
    mat dij, u, u_trans, q, r, q_t, tmp, tmp2, z;
    double* s = (double*) malloc(sizeof(double)*k);
    double* ones = (double*) malloc(sizeof(double)*n);
    double* d;
    int* anchors = (int*) malloc(sizeof(int)*k);
    int* clusters = NULL;
    double change = 1, old_stress = -1;
    dij = mat_new(k, n);
    u = mat_new(n,k);
    tmp = mat_new(n,k);
    darrset(ones,n,-1);
    
    select_anchors(g, dij, anchors, k);
    if(opts.color) {
        for(i = 0; i < k; i++) {
            Agnode_t* anchor = get_node(anchors[i]);
            agset(anchor, "color", "red");
        }
    }
    if(opts.power != 1) {
        clusters = graph_cluster(g,dij,anchors);
    }

    singular_vectors(g, dij, opts.power, u, s);
    vec_scalar_mult(s, k, -1);
    u_trans = mat_trans(u);
    d = mat_mult_for_d(u, s, u_trans, ones);
    for(i = 0; i < u->c; i++) {
        double* col = mat_col(u,i);
        double* b = inv_mul_ax(d,col,u->r);
        for(j = 0; j < u->r; j++) {
            tmp->m[mindex(j,i,tmp)] = b[j];     
        }
        free(b);
        free(col);
    }
    tmp2 = mat_mult(u_trans,tmp);
    for(i = 0; i < k; i++) {
        tmp2->m[mindex(i,i,tmp2)] += (1.0/s[i]);
    }
    q = mat_new(tmp2->r, tmp2->c);
    r = mat_new(tmp2->c, tmp2->c);
    qr_factorize(tmp2,q,r);
    q_t = mat_trans(q);

    if(opts.given) {
        z = get_positions(g, opts.dim);
    } else {
        z = mat_rand(n, opts.dim);
    }
    translate_by_centroid(z);
    
    old_stress = stress(z, dij, anchors, opts.power);
    while(change > EPSILON && iter < opts.max_iter) {
        mat right_side;
        double new_stress;
        
        if(opts.power == 1) {
            right_side = barnes_hut(z);
        } else {
            right_side = barnes_hut_cluster(z, dij, clusters, opts.power);
        }
        for(i = 0; i < opts.dim; i++) {
            double sum = 0;         
            double* x;
            double* b = mat_col(right_side,i);
            for(j = 0; j < right_side->r; j++) {
                sum += b[j];
            }
            x = inv_mul_full(d, b, right_side->r, u, u_trans, q_t, r);
            for(j = 0; j < z->r; j++) {
                z->m[mindex(j,i,z)] = x[j] - sum/right_side->r;
            }
            free(x);
            free(b);
        }
        
        adjust_anchors(g, anchors, k, z);
        update_anchors(z, dij, anchors, opts.power);
        translate_by_centroid(z);
    
        new_stress = stress(z, dij, anchors, opts.power);
        change = fabs(new_stress-old_stress)/old_stress;
        old_stress = new_stress;
        
        mat_free(right_side);
        iter++;
    }
    
    mat_free(dij);
    mat_free(u);
    mat_free(u_trans);
    mat_free(q);
    mat_free(r);
    mat_free(q_t);
    mat_free(tmp);
    mat_free(tmp2);
    free(s);
    free(ones);
    free(d);
    free(anchors);
    free(clusters);
    
    return z;
}
