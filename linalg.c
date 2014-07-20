#include "linalg.h"

mat mat_new(int r, int c)
{
	mat m = (mat) malloc(sizeof(struct matrix));
	m->r = r;
	m->c = c;
	m->m = (double*) malloc(r*c*sizeof(double));
	memset(m->m, 0, sizeof(double)*r*c);
	return m;
}

mat mat_rand(int r, int c)
{
	int i;
	mat m = (mat) malloc(sizeof(struct matrix));
	m->r = r;
	m->c = c;
	m->m = (double*) malloc(r*c*sizeof(double));
	for(i = 0; i < r*c; i++) {
		m->m[i] = rand() / ((double)(RAND_MAX));
	}
	return m;
}

void mat_free(mat m)
{
	free(m->m);
	free(m);
}

void mat_set(mat m, double val)
{
	int i;
	for(i = 0; i < m->r*m->c; i++) {
		m->m[i] = val;
	}
}

mat mat_trans(mat m)
{
	int i,j;
	mat m_trans = mat_new(m->c,m->r);
	for(i = 0; i < m->r; i++) {
		for(j = 0; j < m->c; j++) {
			m_trans->m[mindex(j,i,m_trans)] = m->m[mindex(i,j,m)];
		}	
	}
	return m_trans;
}

void mat_print(mat m, FILE* out)
{
	int i,j;
	for(i = 0; i < m->r; i++) {
		for(j = 0; j < m->c; j++) {
			fprintf(out,"%f ",m->m[mindex(i,j,m)]);
		}
		fprintf(out,"\n");
	}
}

void vec_print(double* v, int n, FILE* out)
{
	int i;
	for(i = 0; i < n; i++) {
		fprintf(out, "%f ", v[i]);
	}
	fprintf(out,"\n");
}

void mat_sub(mat a, mat b)
{
	int i,j;
	for(i = 0; i < a->r; i++) {
		for(j = 0; j < a->c; j++) {
			a->m[mindex(i,j,a)] = a->m[mindex(i,j,a)]-b->m[mindex(i,j,b)];		
		}	
	}
}

void mat_scalar_mult(mat m, double s)
{
	int i;
	for(i = 0; i < m->r*m->c; i++) {
		m->m[i] = s*m->m[i];
	}
}

mat mat_mult(mat a, mat b)
{
	int i,j,k;
	mat res = mat_new(a->r, b->c);
	
	if(a->c != b->r) {
		fprintf(stderr, "Invalid matrix multiplication. Cols of param 1 do not equal rows of param 2.");
		exit(1);
	}
	
	for(i = 0; i < res->r; i++) {
		for(j = 0; j < res->c; j++) {
			double sum = 0;
			for(k = 0; k < a->c; k++) {
				sum += 	a->m[mindex(i,k,a)]*b->m[mindex(k,j,b)];	
			}
			res->m[mindex(i,j,res)] = sum;
		}
	}
	return res;
}

double* mat_vec_mult(mat a, double* b)
{
	double* res = (double*) malloc(sizeof(double)*a->r);
	int i,j;
	for(i = 0; i < a->r; i++) {
		double sum = 0;
		for(j = 0; j < a->c; j++) {
			sum += a->m[mindex(i,j,a)]*b[j];
		}
		res[i] = sum;
	}
	return res;
}

double* mat_col(mat a, int j)
{
	double* col = (double*) malloc(sizeof(double)*a->r);
	int i;
	for(i = 0; i < a->r; i++) {
		col[i] = a->m[mindex(i,j,a)];
	}
	return col;
}

double mat_accu(mat m)
{
	double sum = 0;
	int i;
	for(i = 0; i < m->r*m->c; i++) {
		sum += m->m[i];
	}
	return sum;
}

void vec_scalar_mult(double* v, int n, double s)
{
	int i;
	for(i = 0; i < n; i++) {
		v[i] *= s;
	}
}

int mindex(int i, int j, mat m)
{
	return i*m->c + j;
}
