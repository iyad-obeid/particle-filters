#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <string.h>
#include "createData.h"
#include "oMatrixPrimitives.h"

double max_firing = 100;
extern double nNeurons;
extern double nSamples;
extern gsl_rng *r;
int SIMTYPE = 0;

void create_dormant_center(double *, double *);
void create_jump(double *, double *);
void create_sine(double *, double *);
void create_distractor(double *, double);
void stuff_alpha(gsl_matrix *);
void stuff_alpha_7(gsl_matrix *);
void stuff_sigma(gsl_matrix *);
void stuff_sigma_7(gsl_matrix *);
void stuff_mu   (gsl_matrix *, double *);

void write_out(double nNeurons, double nSamples, double *var, FILE *fid_out){
    int i, cnt = 0;
    for (i=0;i<nNeurons;i++)
        fwrite(var+(cnt+=nSamples) , sizeof(double),1,fid_out);
}

void createData(char *fNameIn)
{
    FILE *fid = fopen(fNameIn,"wb");
    int i,j;
    double dx = 0;
    extern double tMax;
    extern double dt;
    extern char *dataType;
    double *t, *posn, *dstr, tmp;
    double iSecondary = 0;
    double pctErr = 0, pctMissed = 0, pctFA = 0;
    double weightPrimary = 0, weightSecondary = 0;
    extern double dataParam[];

    nSamples   = tMax/dt;
    if (strcmp(dataType,"multineuron_sine")==0){
        SIMTYPE=0;
        iSecondary = INFINITY;
    }
    else if (strcmp(dataType,"mixture_model_type_1")==0){
        SIMTYPE=1;
        iSecondary = round(dataParam[0] * nNeurons); // in this case dataParam[0] tells us the pct of primary neurons
    }
    else if (strcmp(dataType,"multineuron_randwalk")==0){
        SIMTYPE=2;
        iSecondary = INFINITY;
        dx = dataParam[0];
    }
    else if (strcmp(dataType,"randwalk_sort_err")==0){
        SIMTYPE = 3;
        dx = 1.5; // hard wire in a particular value
        iSecondary = INFINITY;
        pctErr = dataParam[0];
    }
    else if (strcmp(dataType,"randwalk_missed_detections")==0){
        SIMTYPE = 4;
        dx = 1.5;
        iSecondary = INFINITY;
        pctMissed = dataParam[0];
        printf("pctMissed = %f\n",pctMissed);
    }
    else if (strcmp(dataType,"randwalk_false_alarms")==0){
        SIMTYPE = 5;
        dx = 1.5;
        iSecondary = INFINITY;
        pctFA = dataParam[0];
    }
    else if (strcmp(dataType,"mixture_model_type_2")==0){
        SIMTYPE = 6;
        dx = 1.5;
        iSecondary = INFINITY;
        weightPrimary = dataParam[0];
        weightSecondary = 1-weightPrimary;
    }
    else if (strcmp(dataType,"confidence_interval_1")==0){
        SIMTYPE = 7;
        // dx = 0.5; // default
        dx = 3; // test case
        iSecondary = INFINITY;
    }
    else if (strcmp(dataType,"confidence_interval_2")==0){
        SIMTYPE = 8;
        dx = 0.5;
        iSecondary = INFINITY;
    }
    else if (strcmp(dataType,"confidence_interval_3")==0){
        SIMTYPE = 9;
        dx = 0.5;
        iSecondary = INFINITY;
    }
    else if (strcmp(dataType,"confidence_interval_4")==0){
        SIMTYPE = 10;
        dx = 0.5;
        iSecondary = INFINITY;
    }

    t    = (double *)malloc(nSamples*sizeof(double));
    posn = (double *)malloc(nSamples*sizeof(double));
    dstr = (double *)malloc(nSamples*sizeof(double));

    gsl_matrix *alpha  = gsl_matrix_alloc(nNeurons,nSamples);
    gsl_matrix *sigma  = gsl_matrix_alloc(nNeurons,nSamples);
    gsl_matrix *mu     = gsl_matrix_alloc(nNeurons,nSamples);
    gsl_matrix *alpha2 = gsl_matrix_alloc(nNeurons,nSamples);
    gsl_matrix *sigma2 = gsl_matrix_alloc(nNeurons,nSamples);
    gsl_matrix *mu2    = gsl_matrix_alloc(nNeurons,nSamples);

    // Initialize t and posn vectors
    t[0] = dt;
    for (i=1;i<nSamples;i++)
        t[i] = t[i-1]+dt;

    if      (SIMTYPE==0) // standard pure sine wave
        create_sine(posn,t);
    else if (SIMTYPE==1) // some are sine wave some are random walk
        create_sine(posn,t);
    else if (SIMTYPE==2) // Random walk
        create_distractor(posn, dx);
    else if (SIMTYPE==3) // rand walk with sorting error
        create_distractor(posn,dx);
    else if (SIMTYPE==4) // rand walk with missed detections
        create_distractor(posn,dx);
    else if (SIMTYPE==5) // rand walk with false alarms
        create_distractor(posn,dx);
    else if (SIMTYPE==6){ // mixture model type 2 - each neuron tuned by both
        for (i=0;i<nSamples;i++)
            posn[i] = 150*sin(2*M_PI/5 * t[i]) + 150;
        create_distractor(dstr,dx);
    }
    else if (SIMTYPE==7) // confidence interval
    	create_distractor(posn,dx);
    else if (SIMTYPE==8) // confidence interval
    	create_sine(posn,t);
    else if (SIMTYPE==9) // confidence interval    	
    	create_dormant_center(posn,t);
    else if (SIMTYPE==10) // confidence interval
    	create_jump(posn,t);


    // stuff values into neuron parameter variables
      if (SIMTYPE<=6){
	stuff_alpha(alpha);
	stuff_sigma(sigma);
	stuff_mu(mu,t);
      }
      else{
	stuff_alpha_7(alpha);
	stuff_sigma_7(sigma);
	stuff_mu(mu,t);
      }

    if (SIMTYPE==6){ // mixture model type 2 - each neuron tuned twice simultaneously
        stuff_alpha(alpha2);
        stuff_sigma(sigma2);
        stuff_mu   (mu2, t);
    }

    // Define N
    gsl_matrix *N;
    N = gsl_matrix_alloc(nNeurons,nSamples);
    double lambda=0,lambda_1,lambda_2;
    for (i=0;i<nNeurons;i++){
        if (i>=iSecondary)
            create_distractor(dstr,1.5);
        for (j=0;j<nSamples;j++){
            if (SIMTYPE==1){
                if (i>=iSecondary)
                    lambda = exp(gsl_matrix_get(alpha,i,j) - pow(dstr[j]-gsl_matrix_get(mu,i,j) , 2) / (2*pow(gsl_matrix_get(sigma,i,j),2)) );
            }
            else if (SIMTYPE==6){
                lambda_1 = exp(gsl_matrix_get(alpha, i,j) - pow(posn[j]-gsl_matrix_get(mu, i,j) , 2) / (2*pow(gsl_matrix_get(sigma, i,j),2)) );
                lambda_2 = exp(gsl_matrix_get(alpha2,i,j) - pow(dstr[j]-gsl_matrix_get(mu2,i,j) , 2) / (2*pow(gsl_matrix_get(sigma2,i,j),2)) );
                lambda = lambda_1*weightPrimary + lambda_2*weightSecondary;
            }
            else
                lambda = exp(gsl_matrix_get(alpha,i,j) - pow(posn[j]-gsl_matrix_get(mu,i,j) , 2) / (2*pow(gsl_matrix_get(sigma,i,j),2)) );
            tmp = gsl_ran_poisson(r,lambda*dt);
            tmp = (tmp>=1) ? 1 : 0;
            gsl_matrix_set(N,i,j,tmp);
        }
    }

    double odd=0, even=0;
    printf("SIMTYPE = %i\n",SIMTYPE);

    int cntLocal = 0;
    if (SIMTYPE==3){
        for (i=0;i<nNeurons;i=i+2){
            for(j=0;j<nSamples;j++){
                odd  = gsl_matrix_get(N,i,j);
                even = gsl_matrix_get(N,i+1,j);
                // technically, you don't need to check if odd or even == 1 before generating rng and swapping
                if (((odd==1) | (even==1)) & (gsl_rng_uniform(r)<=pctErr)) {
                    gsl_matrix_set(N,i,j,even);
                    gsl_matrix_set(N,i+1,j,odd);
                    cntLocal++;
                }
            }
        }
        printf("There were a total of %i swaps.\n",cntLocal);
    }

    if (SIMTYPE==4){
        for (i=0;i<nNeurons;i++)
            for (j=0;j<nSamples;j++)
                    if ((gsl_matrix_get(N,i,j)==1) & (gsl_rng_uniform(r)<=pctMissed)){
                        gsl_matrix_set(N,i,j,0);
                        cntLocal++;
                    }
        printf("There were a total of %i missed detections.\n",cntLocal);
    }

    if (SIMTYPE==5){
        for (i=0;i<nNeurons;i++)
            for (j=0;j<nSamples;j++)
                    if ((gsl_matrix_get(N,i,j)==0) & (gsl_rng_uniform(r)<=pctFA)){
                        gsl_matrix_set(N,i,j,1);
                        cntLocal++;
                    }
        printf("There were a total of %i false alarms.\n",cntLocal);
    }


    // Output variables to file
    fwrite(dataType, sizeof(char), 50,fid);
    fwrite(&tMax,    sizeof(double),1,fid);
    fwrite(dataParam,sizeof(double),5,fid);
    fwrite(&nNeurons,sizeof(double),1,fid);
    fwrite(&dt,      sizeof(double),1,fid);
    fwrite(&nSamples,sizeof(double),1,fid);
    fwrite(posn,     sizeof(double),nSamples,fid);
    o_matrix_fwrite_char(fid,N);
    o_matrix_fwrite(fid,mu);
    if (SIMTYPE<7){
    o_matrix_fwrite_col1(fid,alpha);
    o_matrix_fwrite_col1(fid,sigma);
    }
    else{
      o_matrix_fwrite(fid,alpha);
      o_matrix_fwrite(fid,sigma);
    }

    if (SIMTYPE==6){
        fwrite(dstr,     sizeof(double),nSamples,fid);
        o_matrix_fwrite(fid,mu2);
        o_matrix_fwrite_col1(fid,alpha2);
        o_matrix_fwrite_col1(fid,sigma2);
    }

    gsl_matrix_free(alpha);
    gsl_matrix_free(sigma);
    gsl_matrix_free(mu);
    gsl_matrix_free(N);

    gsl_matrix_free(alpha2);
    gsl_matrix_free(sigma2);
    gsl_matrix_free(mu2);

    free(t);
    free(posn);
    free(dstr);
    fclose(fid);
}

void create_dormant_center(double *posn, double *t)
{
  extern double nSamples;
  int i;

  for (i=0;i<nSamples;i++)
    if (t[i] < 15)
      posn[i] = 100 + (350.0/15.0)*t[i];
    else
      posn[i] = 800 - (350.0/15.0)*t[i];
}

void create_jump(double *posn, double *t)
{
  extern double nSamples;
  int i;
  for (i=0;i<nSamples;i++){
    if (t[i] < 15)
      posn[i] =   50 + (200.0/15.0)*t[i];
    else
      posn[i] = -150 + (200.0/15.0)*t[i];
  }
}

void create_sine(double *posn, double *t)
{
    int i;
    extern double nSamples;
    for (i=0;i<nSamples;i++)
    posn[i] = 150*sin(2*M_PI/5 * t[i]) + 150;
}

void create_distractor(double *dstr, double dx)
{
      extern gsl_rng *r;
      extern double nSamples;
    int i,flag;
    flag = 1;
    i = 1;
    
    gsl_rng *rr;
    gsl_rng *s = gsl_rng_alloc(gsl_rng_mt19937);
    unsigned long int seed = 314159;
    gsl_rng_set(s,seed);
    
    if (SIMTYPE==7)
      rr = s;
    else
      rr = r;

    dstr[0] = 150 + gsl_ran_gaussian(rr,dx);

    while (flag==1){
        dstr[i] = dstr[i-1] + gsl_ran_gaussian(rr,dx);
        if ( (dstr[i]>300) | (dstr[i]<0) )
            i = 1;
        else
            i++;
        if (i==nSamples)
            flag = 0;
    }
}

void stuff_alpha(gsl_matrix *alpha){
    int i,j;
    double tmp;

    for (i=0;i<nNeurons;i++){
        tmp = log(max_firing)*gsl_rng_uniform(r);
        for(j=0;j<nSamples;j++)
            gsl_matrix_set(alpha,i,j,tmp);
    }
}

void stuff_alpha_7(gsl_matrix *alpha){
    int i,j;
    double tmp;
    int flag;
    double aStart,aEnd;
    double dy;

    for (i=0;i<nNeurons;i++){

      flag = 1;
      while (flag){
        aStart = 50*gsl_rng_uniform(r);
		aEnd   = 50*gsl_rng_uniform(r);
	if (fabs(aStart-aEnd)<5)
		flag = 0;
      }

      dy = (aEnd-aStart)/nSamples;
      tmp = aStart;
      for(j=0;j<nSamples;j++){
	gsl_matrix_set(alpha,i,j,log(tmp));
	tmp+=dy;
      }
    }

}

void stuff_sigma_7(gsl_matrix *sigma){
    int i,j;
    double tmp;
    int flag;
    double sStart,sEnd;
    double dy;

    for (i=0;i<nNeurons;i++){
      flag = 1;
      while (flag){
        sStart = 10*gsl_rng_uniform(r) + 10;
	sEnd   = 10*gsl_rng_uniform(r) + 10;
	if (abs(sStart-sEnd)<5)
		flag = 0;
      }

      dy = (sEnd-sStart)/nSamples;
      tmp = sStart;
      for(j=0;j<nSamples;j++)
	gsl_matrix_set(sigma,i,j,tmp+=dy);
    }
}

void stuff_sigma(gsl_matrix *sigma){
    int i,j;
    double tmp;
    for (i=0;i<nNeurons;i++){
        tmp = 10*gsl_rng_uniform(r) + 10;
        for(j=0;j<nSamples;j++)
            gsl_matrix_set(sigma,i,j,tmp);
    }
}

void stuff_mu(gsl_matrix *mu, double *t){
    int i,j;
    double tmp;
    double mu_slope, mu_initial;
    for (i=0;i<nNeurons;i++){
        mu_slope = (2*gsl_rng_uniform(r)-1)/2;
        mu_initial = 400*gsl_rng_uniform(r)-50;
        for(j=0;j<nSamples;j++){
            tmp = mu_initial + mu_slope*t[j];
            gsl_matrix_set(mu,i,j,tmp);
        }
    }
}
