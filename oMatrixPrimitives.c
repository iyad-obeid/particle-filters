#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include "oMatrixPrimitives.h"

unsigned long int random_seed();

void o_matrix_normpdf(gsl_matrix *prior, gsl_matrix *sim , double sd){
  int i,j;
  for (i=0;i<sim->size1;i++)
    for (j=0;j<sim->size2;j++)
      gsl_matrix_set(prior,i,j,gsl_ran_gaussian_pdf(gsl_matrix_get(sim,i,j),sd));
}

void o_matrix_rep_rows(gsl_matrix *dest, gsl_vector *src){
    int i;
    for (i=0;i<dest->size1;i++)
        gsl_matrix_set_row(dest,i,src);
}

void o_matrix_rep_cols(gsl_matrix *dest, gsl_vector *src){
    int i;
    for (i=0;i<dest->size2;i++)
        gsl_matrix_set_col(dest,i,src);
}

void o_matrix_randn(gsl_rng *r, double m, double s, gsl_matrix *x){
    int i,j;
    for (i=0;i<x->size1;i++)
        for(j=0;j<x->size2;j++)
            gsl_matrix_set(x,i,j,gsl_ran_gaussian(r,s)+m);
}

void o_vector_randn(gsl_rng *r, double m, double s, gsl_vector *x){
    int i;
    for (i=0;i<x->size;i++)
        gsl_vector_set(x,i,gsl_ran_gaussian(r,s)+m);
}

double o_vector_sd(gsl_vector *v){
    return(gsl_stats_sd(v->data,1,v->size));
}

double o_vector_mean(gsl_vector *v){
    return(gsl_stats_mean(v->data,1,v->size));
}

double o_matrix_mean_row(gsl_matrix *m,int i){
    return(gsl_stats_mean(m->data+i*m->size2,1,m->size2));
}

double o_matrix_mean_col(gsl_matrix *m,int i){
    return(gsl_stats_mean(m->data+i,m->size2,m->size1));
}

double o_matrix_std_col(gsl_matrix *m,int i){
    return(gsl_stats_sd(m->data+i,m->size2,m->size1));
}

void o_matrix_std1(gsl_vector *dest, gsl_matrix *src){
    int i;
    for (i=0;i<src->size2;i++)
        gsl_vector_set(dest,i,o_matrix_std_col(src,i));
}

void o_matrix_mean1(gsl_vector *dest, gsl_matrix *src){
    int i;
    for(i=0;i<src->size2;i++)
        gsl_vector_set(dest,i,o_matrix_mean_col(src,i));
}

void o_matrix_mean2(gsl_vector *dest, gsl_matrix *src){
    int i;
    for(i=0;i<src->size1;i++)
        gsl_vector_set(dest,i,o_matrix_mean_row(src,i));
}

void o_matrix_fread(FILE *fid, gsl_matrix *m){
    gsl_matrix *n = gsl_matrix_alloc(m->size2,m->size1);
    gsl_matrix_fread(fid,n);
    gsl_matrix_transpose_memcpy(m,n);
    gsl_matrix_free(n);
}

void o_matrix_fwrite_col1(FILE *fid, gsl_matrix *var){
    gsl_matrix_view COL = gsl_matrix_submatrix(var,0,0,var->size1,1);
    gsl_matrix_fwrite(fid,&(COL.matrix));
}

void o_matrix_fread_char(FILE *fid, gsl_matrix *m){
    int i;
    char *N_ch = (char *)malloc(m->size1 * m->size2 * sizeof(char));
    fread(N_ch,m->size1 * m->size2,sizeof(char),fid);
    double *N = (double *)malloc(m->size1 * m->size2 * sizeof(double));
    for (i=0;i<m->size1 * m->size2;i++)
        N[i] = (double)N_ch[i];
    gsl_matrix_view N_view = gsl_matrix_view_array(N,m->size2,m->size1);
    gsl_matrix *m_trans = gsl_matrix_alloc(m->size2,m->size1);
    gsl_matrix_memcpy(m_trans,&(N_view.matrix));
    gsl_matrix_transpose_memcpy(m,m_trans);

    gsl_matrix_free(m_trans);
    free(N_ch);
    free(N);
}

void o_matrix_fwrite_char(FILE *fid, gsl_matrix *m){
    gsl_matrix *n = gsl_matrix_alloc(m->size2,m->size1);
    gsl_matrix_transpose_memcpy(n,m);
    char *N_ch = (char *)malloc(n->size1 * n->size2 * sizeof(char));
    int i;
    for (i=0;i<n->size1 * n->size2;i++)
        N_ch[i] = (char)(n->data[i]);
    fwrite(N_ch,n->size1 * n->size2, sizeof(char),fid);
    free(N_ch);
    gsl_matrix_free(n);
}
void o_matrix_fwrite(FILE *fid, gsl_matrix *m){
    gsl_matrix *n;
    n = gsl_matrix_alloc(m->size2,m->size1);
    gsl_matrix_transpose_memcpy(n,m);
    gsl_matrix_fwrite(fid,n);
    gsl_matrix_free(n);
}

void o_init_rng(gsl_rng *r){
    unsigned long int seed;
    seed = random_seed();
    gsl_rng_set(r,seed);
}

unsigned long int random_seed()
{

 unsigned int seed;
 FILE *devrandom;

 if ((devrandom = fopen("/dev/random","r")) == NULL) {
   fprintf(stderr,"Cannot open /dev/random, setting seed to 0\n");
   seed = 0;
 } else {
   fread(&seed,sizeof(seed),1,devrandom);
//   if(verbose == D_SEED) printf("Got seed %u from /dev/random\n",seed);
   fclose(devrandom);
 }

 return(seed);

}

void o_vector_add(gsl_vector *x, gsl_vector *y, gsl_vector *z)
{
    gsl_vector_memcpy(z,x);
    gsl_vector_add(z,y);
}

void o_matrix_invert(gsl_matrix *x, gsl_matrix *y){
    gsl_matrix *b;
    gsl_permutation *p;
    int sn;

    p = gsl_permutation_alloc(x->size1);
    b = gsl_matrix_alloc(x->size1,x->size2);

    gsl_matrix_memcpy(b,x);
    gsl_linalg_LU_decomp(b,p,&sn);
    gsl_linalg_LU_invert(b,p,y);

    gsl_matrix_free(b);
    gsl_permutation_free(p);
}

double o_vector_sum(gsl_vector *x){
    int i;
    double s = 0;
    for (i=0;i<x->size;i++)
        s+=gsl_vector_get(x,i);
    return(s);
}

void o_vector_cumsum(gsl_vector *dest, gsl_vector *src){
    int i;
    gsl_vector_set(dest,0,gsl_vector_get(src,0));
    for (i=1;i<dest->size;i++)
        gsl_vector_set(dest,i,gsl_vector_get(dest,i-1) + gsl_vector_get(src,i));
}

void o_randsample(gsl_rng *r, gsl_vector *dest, gsl_vector *w){
    double u;
    int flag,i,j,k;
    gsl_vector *g = gsl_vector_alloc(w->size);
    // assumes w's are already normalized

    o_vector_cumsum(g,w);

    for (i=0;i<w->size;i++){
      gsl_vector_set(dest,i,-1);
        u = gsl_rng_uniform(r);
        flag = 1;
        j = 0;
        while (flag){
            if (j==g->size){
                // if j exceeds the length of g for whatever reason (ie g filled with nans or infs) just pick an index at random
                flag = 0;
                j = floor(u*g->size);
                gsl_vector_set(dest,i,j);
                printf("randsample range exceeded!\n");
            }
            else if (u<gsl_vector_get(g,j)){
                flag = 0;
                gsl_vector_set(dest,i,j);
            }
            else
                j++;
        }
    }

    gsl_vector_free(g);
}

void o_matrix_replace_row(gsl_matrix *dest, int i, gsl_matrix *src, int j){
    gsl_vector *tmp = gsl_vector_alloc(dest->size2);
    gsl_matrix_get_row(tmp,src,j);
    gsl_matrix_set_row(dest,i,tmp);
    gsl_vector_free(tmp);
}

