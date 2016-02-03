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

double hh = 0.01; // default
double nParticles = 100;
double alpha_est = 4;
double sigma_est = 15;
double blocksz = 25;
double nRepeats =  5;
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

void liu_step1(gM *mu_est, gM *posn_est, gM *mu_sMu, gM *mu_sTheta, gM *mu_m, gM *posn_sMu, gM *posn_sTheta, gM *posn_m){

    double a = sqrt(1-pow(hh,2));
    int i, NCOLS;
    double SD;

    gM *sTheta_local;
    gM *mLocal;
    gM *theta_local;
    gM *sMu_local;
    gV *tmp1;
    gM *tmp2;
    gM *x_local;

    for (i=0;i<2;i++){

    if (i==0){ // MU
        SD = 1;
        NCOLS = nNeurons;
        x_local = mu_est;
    }
    else { // POSN
        SD = 20; // default
        NCOLS = 1;
        x_local = posn_est;
    }

    sTheta_local = gsl_matrix_alloc(nParticles,NCOLS);
    mLocal       = gsl_matrix_alloc(nParticles,NCOLS);
    theta_local  = gsl_matrix_alloc(nParticles,NCOLS);
    sMu_local    = gsl_matrix_alloc(nParticles,NCOLS);
    tmp1         = gsl_vector_alloc(NCOLS);
    tmp2         = gsl_matrix_alloc(nParticles,NCOLS);

    o_matrix_randn(r,0,SD,sTheta_local);

    gsl_matrix_memcpy(mLocal,sTheta_local);
    gsl_matrix_scale(mLocal, a);
    o_matrix_mean1(tmp1,sTheta_local);
    gsl_vector_scale(tmp1,(1-a));
    o_matrix_rep_rows(tmp2,tmp1);
    gsl_matrix_add(mLocal,tmp2); // mLocal is defined

    o_matrix_std1(tmp1,sTheta_local);
    o_matrix_rep_rows(theta_local,tmp1);
    gsl_matrix_scale(theta_local,hh);
    o_matrix_randn(r,0,1,tmp2);
    gsl_matrix_mul_elements(theta_local,tmp2); // theta_local is defined

    gsl_matrix_memcpy(sMu_local,theta_local);
    gsl_matrix_add(sMu_local,x_local); // sMu_local is defined

    if (i==0){
        gsl_matrix_memcpy(mu_m,mLocal);
        gsl_matrix_memcpy(mu_sTheta,sTheta_local);
        gsl_matrix_memcpy(mu_sMu,sMu_local);
    }
    else {
        gsl_matrix_memcpy(posn_m,mLocal);
        gsl_matrix_memcpy(posn_sTheta,sTheta_local);
        gsl_matrix_memcpy(posn_sMu,sMu_local);
    }

    gsl_matrix_free(sTheta_local);
    gsl_matrix_free(mLocal);
    gsl_matrix_free(theta_local);
    gsl_matrix_free(sMu_local);
    gsl_vector_free(tmp1);
    gsl_matrix_free(tmp2);

    }
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

void liu_step2a(gV *wg, gM *N, gM *posn_sMu, gM *mu_sMu){

    double sum;

    liu_compute_likelihood(wg,N,posn_sMu,mu_sMu);

    sum = o_vector_sum(wg);
    if (sum!=0)
        gsl_vector_scale(wg,1/sum);
    else
        gsl_vector_set_all(wg,1/nParticles);
}

void liu_step2b(gV *wg, gM *posn_sTheta, gM *posn_sMu, gM *posn_est, gM *posn_m, gM *mu_sTheta, gM *mu_sMu, gM *mu_est, gM *mu_m, gM *rpfx_mu, gM *rpfx_posn){
    int i,j;
    gV *resample_indices = gsl_vector_alloc(wg->size);

    gM *posn_sTheta_tmp = gsl_matrix_alloc(posn_sTheta->size1,posn_sTheta->size2);
    gM *posn_sMu_tmp    = gsl_matrix_alloc(posn_sTheta->size1,posn_sTheta->size2);
    gM *posn_est_tmp    = gsl_matrix_alloc(posn_sTheta->size1,posn_sTheta->size2);
    gM *posn_m_tmp      = gsl_matrix_alloc(posn_sTheta->size1,posn_sTheta->size2);

    gM *mu_sTheta_tmp   = gsl_matrix_alloc(mu_sTheta->size1,mu_sTheta->size2);
    gM *mu_sMu_tmp      = gsl_matrix_alloc(mu_sTheta->size1,mu_sTheta->size2);
    gM *mu_est_tmp      = gsl_matrix_alloc(mu_sTheta->size1,mu_sTheta->size2);
    gM *mu_m_tmp        = gsl_matrix_alloc(mu_sTheta->size1,mu_sTheta->size2);

    gM *rpfx_mu_tmp= gsl_matrix_alloc(mu_sTheta->size1,mu_sTheta->size2);
    gM *rpfx_posn_tmp= gsl_matrix_alloc(posn_sTheta->size1,posn_sTheta->size2);

    o_randsample(r,resample_indices,wg);

    for (i=0;i<nParticles;i++){
        j = (int)gsl_vector_get(resample_indices,i);

        o_matrix_replace_row(posn_sTheta_tmp, i, posn_sTheta, j);
        o_matrix_replace_row(posn_sMu_tmp,    i, posn_sMu,    j);
        o_matrix_replace_row(posn_est_tmp,    i, posn_est,    j);
        o_matrix_replace_row(posn_m_tmp,      i, posn_m,      j);

        o_matrix_replace_row(mu_sTheta_tmp, i, mu_sTheta, j);
        o_matrix_replace_row(mu_sMu_tmp,    i, mu_sMu,    j);
        o_matrix_replace_row(mu_est_tmp,    i, mu_est,    j);
        o_matrix_replace_row(mu_m_tmp,      i, mu_m,      j);

        o_matrix_replace_row(rpfx_mu_tmp,   i, rpfx_mu,j);
        o_matrix_replace_row(rpfx_posn_tmp, i, rpfx_posn,j);
    }

    gsl_matrix_memcpy(posn_sTheta, posn_sTheta_tmp);
    gsl_matrix_memcpy(posn_sMu   , posn_sMu_tmp);
    gsl_matrix_memcpy(posn_est   , posn_est_tmp);
    gsl_matrix_memcpy(posn_m     , posn_m_tmp);

    gsl_matrix_memcpy(mu_sTheta, mu_sTheta_tmp);
    gsl_matrix_memcpy(mu_sMu   , mu_sMu_tmp);
    gsl_matrix_memcpy(mu_est   , mu_est_tmp);
    gsl_matrix_memcpy(mu_m     , mu_m_tmp);

    gsl_matrix_memcpy(rpfx_mu, rpfx_mu_tmp);
    gsl_matrix_memcpy(rpfx_posn, rpfx_posn_tmp);

    gsl_vector_free(resample_indices);
    gsl_matrix_free(posn_sTheta_tmp);
    gsl_matrix_free(posn_sMu_tmp);
    gsl_matrix_free(posn_est_tmp);
    gsl_matrix_free(posn_m_tmp);

    gsl_matrix_free(mu_sTheta_tmp);
    gsl_matrix_free(mu_sMu_tmp);
    gsl_matrix_free(mu_est_tmp);
    gsl_matrix_free(mu_m_tmp);

    gsl_matrix_free(rpfx_mu_tmp);
    gsl_matrix_free(rpfx_posn_tmp);
}

void liu_step3(gM *posn_sTheta,gM *mu_sTheta, gM *posn_m, gM *mu_m){
    int i,j;
    double sigma;

    sigma = sqrt(hh)*o_matrix_std_col(posn_sTheta,0);
    for (i=0;i<posn_m->size1;i++)
        gsl_matrix_set(posn_sTheta,i,0,gsl_matrix_get(posn_m,i,0)+gsl_ran_gaussian(r,sigma));

    for (j=0;j<mu_m->size2;j++){
        sigma = sqrt(hh)*o_matrix_std_col(mu_sTheta,j);
        for (i=0;i<mu_m->size1;i++)
            gsl_matrix_set(mu_sTheta,i,j,gsl_matrix_get(mu_m,i,j)+gsl_ran_gaussian(r,sigma));
    }

}

void liu_step4(gM *posn_est, gM *mu_est, gM *posn_sTheta, gM *mu_sTheta, gM *posn_est_new, gM *mu_est_new){
    int i,j;

    for (i=0;i<posn_est->size1;i++)
    // default
	gsl_matrix_set(posn_est_new,i,0,gsl_matrix_get(posn_est,i,0) + gsl_matrix_get(posn_sTheta,i,0) + gsl_ran_gaussian(r,20));

    for (j=0;j<mu_est->size2;j++)
        for (i=0;i<mu_est->size1;i++)
        // default
		gsl_matrix_set(mu_est_new,i,j, gsl_matrix_get(mu_est,i,j) + gsl_matrix_get(mu_sTheta,i,j) + gsl_ran_gaussian(r,1));
}

void liu_assign_next_particle(gM *rpfx_posn, gM *posn_est_new, gM *rpfx_mu, gM *mu_est_new){
    int i,j;

    for (i=0;i<rpfx_posn->size1;i++)
        gsl_matrix_set(rpfx_posn,i,0,gsl_matrix_get(posn_est_new,i,0));

    for (j=0;j<rpfx_mu->size2;j++)
        for (i=0;i<rpfx_mu->size1;i++)
            gsl_matrix_set(rpfx_mu,i,j, gsl_matrix_get(mu_est_new,i,j) );

}

void liu_step5(gV *wg, gM *N, gM *posn_sMu, gM *mu_sMu, gM *posn_est_new, gM *mu_est_new){
    int i;
    double s;

    gV *num = gsl_vector_alloc(nParticles);
    gV *den = gsl_vector_alloc(nParticles);

    gsl_vector_set_all(num,1);
    gsl_vector_set_all(den,1);

    liu_compute_likelihood(num, N, posn_est_new, mu_est_new);
    liu_compute_likelihood(den, N, posn_sMu, mu_sMu);

    gsl_vector_div(num,den);
    gsl_vector_memcpy(wg,num);

    // don't divide by zero!
    for (i=0;i<wg->size;i++)
        if (gsl_vector_get(den,i)==0)
            gsl_vector_set(wg,i,0);

    s = o_vector_sum(wg);
    if (s!=0)
        gsl_vector_scale(wg,1/s);
    else
        gsl_vector_set_all(wg,1/nParticles);

    gsl_vector_free(num);
    gsl_vector_free(den);
}

void liu_assign_next_estimate(gV *wg, gM *rpfx_posn, gV *xEstimate_posn, gM *rpfx_mu, gM *xEstimate_mu){
    int i,j;
    double s=0;

    for (i=0;i<wg->size;i++)
        s += gsl_vector_get(wg,i) * gsl_matrix_get(rpfx_posn,i,0);
    gsl_vector_set(xEstimate_posn,blockIndex+1,s);

    for (j=0;j<xEstimate_mu->size2;j++){
        s=0;
        for (i=0;i<wg->size;i++)
            s += gsl_vector_get(wg,i) * gsl_matrix_get(rpfx_mu,i,j);
        gsl_matrix_set(xEstimate_mu,blockIndex+1,j,s);
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
    o_init_rng(r);

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
    sprintf(fNameOut,"%sd%05iout.bin",outDir,filenum);
    if (SIM7FLAG){
      filenum = dataParam[3];
      printf("Overriding input file\n");
    }
    sprintf(fNameIn,"%sd%05iinp.bin",outDir,filenum);

    // create data file
    if (SIM7FLAG == 0)
      createData(fNameIn);

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
    gM *mu_sTheta   = gsl_matrix_alloc(nParticles,nNeurons);
    gM *mu_m        = gsl_matrix_alloc(nParticles,nNeurons);
    gM *posn_sMu    = gsl_matrix_alloc(nParticles,1);
    gM *posn_sTheta = gsl_matrix_alloc(nParticles,1);
    gM *posn_m      = gsl_matrix_alloc(nParticles,1);

    printf("nParticles = %i\n",(int)nParticles);

    for (blockIndex=0; blockIndex<(gen-1); blockIndex++){
        if ((blockIndex%10)==0)
            printf("starting blockIndex %3i of %3i, file %s\n",blockIndex+1, (int)gen, fNameOut);

        gsl_matrix_memcpy(mu_est,rpfx_mu);
        gsl_matrix_memcpy(posn_est,rpfx_posn);

        liu_step1(mu_est, posn_est, mu_sMu, mu_sTheta, mu_m, posn_sMu, posn_sTheta, posn_m);

    for (qq=0; qq<nRepeats; qq++){
        liu_step2a(wg,N,posn_sMu,mu_sMu);
        liu_step2b(wg,posn_sTheta,posn_sMu,posn_est,posn_m,mu_sTheta,mu_sMu,mu_est,mu_m,rpfx_mu,rpfx_posn);
        liu_step3(posn_sTheta,mu_sTheta,posn_m,mu_m);
        liu_step4(posn_est,mu_est,posn_sTheta,mu_sTheta,posn_est_new,mu_est_new);
        liu_assign_next_particle(rpfx_posn,posn_est_new,rpfx_mu,mu_est_new);
        liu_step5(wg,N,posn_sMu,mu_sMu,posn_est_new,mu_est_new);
        }

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
    gsl_matrix_free(mu_sTheta);
    gsl_matrix_free(mu_m);
    gsl_matrix_free(posn_sMu);
    gsl_matrix_free(posn_sTheta);
    gsl_matrix_free(posn_m);
    return(0);
}

