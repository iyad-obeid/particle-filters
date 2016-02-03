#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "oMatrixPrimitives.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <dirent.h>
#include "createData.h"

//double hh = 0.01;
double nParticles = 100;
double alpha_est = 4;
double sigma_est = 15;
double blocksz = 25;
double dt = 2e-3;
double nNeurons;
char   *dataType;
double dataParam[5];
char   fNameIn[50];
char   fNameOut[50];
double tMax;
double gen;
double nSamples;
int blockIndex;
gsl_rng *r;
char outDir[] = "../results/";
int nFNDigits = 5;

void liu_step1(gV *prior, gM* mu_sim, gM *posn_sim, gM *mu_est, gM *posn_est, gM *mu_sMu, gM *posn_sMu){

  int i,j;
  double tmp;
//  double A = 1;
//  double B = 5; // test case
  
  double A = 0.05;
  double B = 10;
    
  gM *mu_prior   = gsl_matrix_alloc(nParticles,nNeurons);
  gM *posn_prior = gsl_matrix_alloc(nParticles,1);
  
  o_matrix_randn(r,0, A,mu_sim);
  o_matrix_randn(r,0, B,posn_sim);
  
  gsl_matrix_memcpy(mu_sMu   , mu_est);
  gsl_matrix_memcpy(posn_sMu , posn_est);
  
  gsl_matrix_add(mu_sMu   , mu_sim);
  gsl_matrix_add(posn_sMu , posn_sim);
  
  o_matrix_normpdf(mu_prior  , mu_sim   , A);
  o_matrix_normpdf(posn_prior, posn_sim , B);

  // Multiply all of the priors for each particle
  for (i=0;i<prior->size;i++){
    tmp = 1;
    for (j=0;j<mu_prior->size2;j++)
      tmp *= gsl_matrix_get(mu_prior,i,j);
    tmp *= gsl_matrix_get(posn_prior,i,0);
    gsl_vector_set(prior,i,tmp);
  }
  
  tmp = o_vector_sum(prior);
  if (tmp!=0)
    gsl_vector_scale(prior,1/tmp);
  else
    gsl_vector_set_all(prior,1/nParticles);

  gsl_matrix_free(mu_prior);
  gsl_matrix_free(posn_prior);
  
}

void liu_compute_likelihood(gV *wg, gM *N, gM *posn_sMu, gM *mu_sMu) {
    double subBlock_start = blockIndex*blocksz;
    int i,j,k;
    double lambda;

    for (i=0;i<nNeurons;i++)
        for (j=0;j<nParticles;j++)
            for (k=subBlock_start;k<subBlock_start + blocksz;k++) {
                lambda = dt * exp(alpha_est - pow( (gsl_matrix_get(posn_sMu,j,0) - gsl_matrix_get(mu_sMu,j,i) ) , 2) / (2*pow(sigma_est,2)));
                if (gsl_matrix_get(N,i,k)==0)
                    gsl_vector_set(wg,j,exp(-lambda) * gsl_vector_get(wg,j));
                else
                    gsl_vector_set(wg,j,lambda * exp(-lambda) * gsl_vector_get(wg,j));
            }
}

void liu_step2a(gV *wg, gM *N, gM *posn_sMu, gM *mu_sMu, gV *prior){

    double sum;

    liu_compute_likelihood(wg,N,posn_sMu,mu_sMu);
    gsl_vector_mul(wg,prior);

    sum = o_vector_sum(wg);
    if (sum!=0)
        gsl_vector_scale(wg,1/sum);
    else
        gsl_vector_set_all(wg,1/nParticles);

}

void liu_step2b(gV *wg, gM *mu_sMu, gM *mu_est, gM *rpfx_mu, gM *posn_sMu, gM *posn_est, gM *rpfx_posn){

    int i,j;
    gV *resample_indices = gsl_vector_alloc(wg->size);

    gM *posn_sMu_tmp    = gsl_matrix_alloc(posn_sMu->size1,posn_sMu->size2);
    gM *posn_est_tmp    = gsl_matrix_alloc(posn_sMu->size1,posn_sMu->size2);

    gM *mu_sMu_tmp      = gsl_matrix_alloc(mu_sMu->size1,mu_sMu->size2);
    gM *mu_est_tmp      = gsl_matrix_alloc(mu_sMu->size1,mu_sMu->size2);

    gM *rpfx_mu_tmp   = gsl_matrix_alloc(mu_sMu  ->size1 , mu_sMu  ->size2);
    gM *rpfx_posn_tmp = gsl_matrix_alloc(posn_sMu->size1 , posn_sMu->size2);

    gV *wg_tmp = gsl_vector_alloc(wg -> size);

    o_randsample(r,resample_indices,wg);

    for (i=0;i<nParticles;i++){
        j = (int)gsl_vector_get(resample_indices,i);

        o_matrix_replace_row(posn_sMu_tmp,    i, posn_sMu,    j);
        o_matrix_replace_row(posn_est_tmp,    i, posn_est,    j);

        o_matrix_replace_row(mu_sMu_tmp,    i, mu_sMu,    j);
        o_matrix_replace_row(mu_est_tmp,    i, mu_est,    j);

        o_matrix_replace_row(rpfx_mu_tmp,   i, rpfx_mu,j);
        o_matrix_replace_row(rpfx_posn_tmp, i, rpfx_posn,j);

	gsl_vector_set(wg_tmp, i, gsl_vector_get(wg,j));
    }

    gsl_matrix_memcpy(posn_sMu   , posn_sMu_tmp);
    gsl_matrix_memcpy(posn_est   , posn_est_tmp);

    gsl_matrix_memcpy(mu_sMu   , mu_sMu_tmp);
    gsl_matrix_memcpy(mu_est   , mu_est_tmp);

    gsl_matrix_memcpy(rpfx_mu, rpfx_mu_tmp);
    gsl_matrix_memcpy(rpfx_posn, rpfx_posn_tmp);

    gsl_vector_memcpy(wg, wg_tmp);

    gsl_vector_free(resample_indices);
    gsl_matrix_free(posn_sMu_tmp);
    gsl_matrix_free(posn_est_tmp);

    gsl_matrix_free(mu_sMu_tmp);
    gsl_matrix_free(mu_est_tmp);

    gsl_matrix_free(rpfx_mu_tmp);
    gsl_matrix_free(rpfx_posn_tmp);
    gsl_vector_free(wg_tmp);
}

void liu_assign_next_particle(gM *rpfx_posn, gM *posn_est_new, gM *rpfx_mu, gM *mu_est_new){
    int i,j;

    for (i=0;i<rpfx_posn->size1;i++)
        gsl_matrix_set(rpfx_posn,i,0,gsl_matrix_get(posn_est_new,i,0));

    for (j=0;j<rpfx_mu->size2;j++)
        for (i=0;i<rpfx_mu->size1;i++)
            gsl_matrix_set(rpfx_mu,i,j, gsl_matrix_get(mu_est_new,i,j) );

}

void liu_assign_next_estimate(gV *wg, gM *rpfx_posn, gV *xEstimate_posn, gM *rpfx_mu, gM *xEstimate_mu){
    int i,j;
    double s=0;
    double wgSum = o_vector_sum(wg);

    for (i=0;i<wg->size;i++)
        s += gsl_vector_get(wg,i) * gsl_matrix_get(rpfx_posn,i,0);
    gsl_vector_set(xEstimate_posn,blockIndex+1,s/wgSum);

    for (j=0;j<xEstimate_mu->size2;j++){
        s=0;
        for (i=0;i<wg->size;i++)
            s += gsl_vector_get(wg,i) * gsl_matrix_get(rpfx_mu,i,j);
        gsl_matrix_set(xEstimate_mu,blockIndex+1,j,s/wgSum);
    }
}

void loadData(double *PTRS_A[], gM *PTRS_B[]){

    FILE *fid = fopen(fNameIn,"rb");
    char foo[50];

    fread(foo,sizeof(char),50,fid);
    fread(&dt,sizeof(double),1,fid);
    fread(dataParam,sizeof(double),5,fid);
    fread(&nNeurons,sizeof(double),1,fid);
    fread(&dt,sizeof(double),1,fid);
    fread(&nSamples,sizeof(double),1,fid);

    PTRS_A[0] = (double *)malloc(nSamples*sizeof(double));

    PTRS_B[0] = gsl_matrix_alloc(nNeurons,nSamples);
    PTRS_B[1] = gsl_matrix_alloc(nNeurons,nSamples);

    fread(PTRS_A[0],sizeof(double),nSamples,fid); // posn

    o_matrix_fread_char(fid,PTRS_B[0]); // N
    o_matrix_fread(fid,PTRS_B[1]); // mu

    fclose(fid);
}

double getenvNumeric(char *env){
    char *str;
    double val;
    if ((str=getenv(env))==NULL)
        val = -1;
    else
        val = atof(str);
    return(val);
}

double getNextFilenum(){
    DIR *pdir = opendir(outDir);
    struct dirent *pent;
    int maxFN = 0;
    int currFN, newFN;
    char theNum[nFNDigits];
    char prefix[] = "d";
    int len = strlen(prefix);

    while ((pent=readdir(pdir))){
        if (strncmp(pent->d_name,prefix,len)==0){
            strncpy(theNum,pent->d_name+len,nFNDigits);
            currFN = atoi(theNum);
            if (currFN > maxFN)
                maxFN = currFN;
        }
    }

    closedir(pdir);
    newFN = maxFN+1;

    if (dataParam[4] != -1){
        newFN = dataParam[4];
        printf("overriding natural numbering sequence!\n");
    }

    return(newFN);
}

int main (int argc, char *argv[], char *envp[]){
    r = gsl_rng_alloc(gsl_rng_mt19937);

    o_init_rng(r); // reinstate this at some point
    //gsl_rng_set(r,(unsigned long int)-2004896911);

    // ver1 is the original version of the code
    // ver2 is where the initial guesses for rpfx have been modified to TRUE_VALUE +/- 5%

    // initialize variables
    int i,j,qq,filenum;
    char   codeVer[10] = "ver2";
    time_t rawtime;
    struct tm *timeinfo;
    char    timechar[50];
    double *PTRS_A[1];
    gM     *PTRS_B[4];
    double *posn;
    gM     *N, *mu_true;
    int SIM7FLAG = 0;

    // Grab environment variables and stash in global vars
    nNeurons     = getenvNumeric("O_nNeurons");
    dataType     = getenv("O_dataStyle");
    dataParam[0] = getenvNumeric("O_dataParam0");
    dataParam[1] = getenvNumeric("O_dataParam1");
    dataParam[2] = getenvNumeric("O_dataParam2");
    dataParam[3] = getenvNumeric("O_dataParam3");
    dataParam[4] = getenvNumeric("O_dataParam4");
    tMax         = getenvNumeric("O_tMax");

    if (strncmp(dataType,"confidence_interval",19)==0){
      SIM7FLAG = 1;
      alpha_est = 3.5;
    }

    // grab the date & time - use "asctime" to print them to file
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    sprintf(timechar,"%s",asctime(timeinfo));

    filenum = getNextFilenum();
    filenum = 9501;
    sprintf(fNameOut,"%sd%05ivpf.bin",outDir,filenum); // default
	
	filenum = dataParam[3];
    sprintf(fNameIn,"%sd%05iinp.bin",outDir,filenum);

    // load data
    loadData(PTRS_A,PTRS_B);
	
    posn    = PTRS_A[0];
    N       = PTRS_B[0];
    mu_true = PTRS_B[1];

    gen  = floor(nSamples/blocksz);

    // allocate memory for rpfx_posn
    gM *rpfx_posn = gsl_matrix_alloc(nParticles,1);

    // initialize rpfx_posn
    for (i=0;i<nParticles;i++)
        gsl_matrix_set(rpfx_posn,i,0, posn[0] +  30*(gsl_rng_uniform(r)-0.5) ); // add pm5% err

    // allocate memory for rpfx_mu
    gM *rpfx_mu = gsl_matrix_alloc(nParticles,nNeurons);

    // initialize rpfx_mu
    for (j=0;j<nNeurons;j++)
        for (i=0;i<nParticles;i++)
            gsl_matrix_set(rpfx_mu,i,j,gsl_matrix_get(mu_true,j,0) + 30*(gsl_rng_uniform(r)-0.5) );

    // allocate and initialize w
    gV *wg = gsl_vector_alloc(nParticles);
    gsl_vector_set_all(wg,1/nParticles);

    // allocate and initialize xEstimate
    gM *xEstimate_mu = gsl_matrix_alloc(gen,nNeurons);
    gV *xEstimate_posn = gsl_vector_alloc(gen);

    gsl_vector_set(xEstimate_posn, 0 ,o_matrix_mean_col(rpfx_posn,0));
    for (i=0;i<nNeurons;i++)
        gsl_matrix_set(xEstimate_mu,0,i,o_matrix_mean_col(rpfx_mu,i));

    gM *mu_est      = gsl_matrix_alloc(nParticles,nNeurons);
    gM *posn_est    = gsl_matrix_alloc(nParticles,1);
    gM *mu_est_new  = gsl_matrix_alloc(nParticles,nNeurons);
    gM *posn_est_new= gsl_matrix_alloc(nParticles,1);
    gM *mu_sMu      = gsl_matrix_alloc(nParticles,nNeurons);
    gM *posn_sMu    = gsl_matrix_alloc(nParticles,1);
    gM *mu_sim      = gsl_matrix_alloc(nParticles,nNeurons);
    gM *posn_sim    = gsl_matrix_alloc(nParticles,1);
    gV *prior       = gsl_vector_alloc(nParticles);

    printf("nParticles = %i\n",(int)nParticles);

    for (blockIndex=0; blockIndex<(gen-1); blockIndex++){
      if ((blockIndex%100)==0)
	printf("starting blockIndex %3i of %3i, file %s\n",blockIndex+1, (int)gen, fNameOut);

      gsl_matrix_memcpy(mu_est,rpfx_mu);
      gsl_matrix_memcpy(posn_est,rpfx_posn);

      liu_step1(prior, mu_sim, posn_sim, mu_est, posn_est, mu_sMu, posn_sMu);
      liu_step2a(wg,N,posn_sMu,mu_sMu,prior);

      liu_step2b(wg,mu_sMu,mu_est,rpfx_mu,posn_sMu,posn_est,rpfx_posn);

      liu_assign_next_particle(rpfx_posn,posn_sMu,rpfx_mu,mu_sMu);
      liu_assign_next_estimate(wg,rpfx_posn,xEstimate_posn,rpfx_mu,xEstimate_mu);

    }

    FILE *fid = fopen(fNameOut,"wb");

    fwrite(codeVer,sizeof(char),10,fid);
    fwrite(timechar,sizeof(char),50,fid);
    fwrite(fNameIn,sizeof(char),50,fid);
    fwrite(&nNeurons,sizeof(double),1,fid);
    fwrite(&nSamples,sizeof(double),1,fid);
    fwrite(&dt,sizeof(double),1,fid);
    fwrite(&blocksz,sizeof(double),1,fid);
    fwrite(&nParticles,sizeof(double),1,fid);
    fwrite(&gen,sizeof(double),1,fid);
    gsl_vector_fwrite(fid,xEstimate_posn);
    gsl_matrix_fwrite(fid,xEstimate_mu);
    fclose(fid);

    free(posn);
    gsl_rng_free(r);
    gsl_matrix_free(rpfx_posn);
    gsl_matrix_free(rpfx_mu);
    gsl_matrix_free(N);
    gsl_matrix_free(mu_true);
    gsl_matrix_free(xEstimate_mu);
    gsl_vector_free(xEstimate_posn);
    gsl_vector_free(wg);
    gsl_matrix_free(mu_est);
    gsl_matrix_free(posn_est);
    gsl_matrix_free(mu_est_new);
    gsl_matrix_free(posn_est_new);
    gsl_matrix_free(mu_sMu);
    gsl_matrix_free(posn_sMu);
    gsl_matrix_free(mu_sim);
    gsl_matrix_free(posn_sim);
    gsl_vector_free(prior);
    return(0);
}

