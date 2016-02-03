#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../c_code_liu/oMatrixPrimitives.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <dirent.h>

double hh = 0.01;
double nParticles = 100;
double alpha_est = 4;
double sigma_est = 15;
double blocksz = 25;
double nRepeats = 5;
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

void loadData(double *PTRS_A[], gM *PTRS_B[]){

    FILE *fid = fopen(fNameIn,"rb");
    char foo[50];

    fread(foo,sizeof(char),50,fid);
    fread(&tMax,sizeof(double),1,fid);
    fread(dataParam,sizeof(double),5,fid);
    fread(&nNeurons,sizeof(double),1,fid);
    fread(&dt,sizeof(double),1,fid);
    fread(&nSamples,sizeof(double),1,fid);

    PTRS_A[0] = (double *)malloc(nSamples*sizeof(double));
    PTRS_A[1] = (double *)malloc(nSamples*sizeof(double));

    PTRS_B[0] = gsl_matrix_alloc(nNeurons,nSamples);
    PTRS_B[1] = gsl_matrix_alloc(nNeurons,nSamples);
    PTRS_B[2] = gsl_matrix_alloc(nNeurons,nSamples);
    PTRS_B[3] = gsl_matrix_alloc(nNeurons,nSamples);

    fread(PTRS_A[0],sizeof(double),nSamples,fid); // t
    fread(PTRS_A[1],sizeof(double),nSamples,fid); // posn

    o_matrix_fread(fid,PTRS_B[0]); // N
    o_matrix_fread(fid,PTRS_B[1]); // alpha
    o_matrix_fread(fid,PTRS_B[2]); // mu
    o_matrix_fread(fid,PTRS_B[3]); // sigma

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

void buildLagMatrix(gM *dst[], gM *src, int lag){
    int nTimes = src->size1;
    int nCols  = src->size2;
    int i,j,k,h;
    int rngStart;
    int dstSize1 = nTimes-lag+1;
    int dstSize2 = lag*nCols+1;

    *dst = gsl_matrix_calloc(dstSize1 , dstSize2);

    rngStart = 0;

    k = 1;
    for (i=0;i<lag;i++){
        for (j=0;j<nCols;j++){
            for(h=0;h<dstSize1;h++)
                gsl_matrix_set(*dst,h,k,gsl_matrix_get(src,rngStart+h,j));
            k++;
        }
        rngStart++;
    }

    for (h=0;h<dstSize1;h++)
        gsl_matrix_set(*dst,h,0,1);
}

int main (int argc, char *argv[], char *envp[]){
    r = gsl_rng_alloc(gsl_rng_mt19937);
    o_init_rng(r);
    double *PTRS_A[2];
    gM     *PTRS_B[4];
    double *t,*posn;
    gM *N, *mu_true;
    int i,j,k;
    int nLags=10;

    sprintf(fNameIn,"../results/d00200inp.bin");
    loadData(PTRS_A,PTRS_B);

    t       = PTRS_A[0];
    posn    = PTRS_A[1];
    N       = PTRS_B[0];
    mu_true = PTRS_B[2];

    double windowDur = 100e-3;
    int winWidth = floor(windowDur / dt);

    double nFRblocks = floor(nSamples / winWidth);
    gM *FR = gsl_matrix_calloc(nFRblocks,nNeurons);

    int sum = 0;
    for (i=0;i<nNeurons;i++)
        for (j=0;j<nFRblocks;j++){
            sum = 0;
            for (k=0;k<winWidth;k++)
                sum += gsl_matrix_get(N,i,(j*winWidth+k));
            gsl_matrix_set(FR,j,i,sum);
        }

    int rTrainLast = 50; // number of training BLOCKS

//    gM *frTrain = gsl_matrix_alloc(           rTrainLast -nLags+1 , nLags*nNeurons+1);
//    gM *frTest  = gsl_matrix_alloc((nFRblocks-rTrainLast)-nLags+1 , nLags*nNeurons+1);

    gM *frTrain[1];
    gM *frTest[1];

    gsl_matrix_view V1;
    gsl_matrix *M1;

    V1 = gsl_matrix_submatrix(FR,0,0,rTrainLast,FR->size2);
    M1 = &(V1.matrix);
    buildLagMatrix(frTrain,M1,nLags);

    V1 = gsl_matrix_submatrix(FR,rTrainLast,0,nFRblocks-rTrainLast,FR->size2);
    M1 = &(V1.matrix);
    buildLagMatrix(frTest,M1,nLags);

    gM *posnTrain = gsl_matrix_calloc(frTrain[0]->size1,1);
    gM *posnTest  = gsl_matrix_calloc(frTest[0] ->size1,1);

    int offset;
    offset = winWidth*(rTrainLast - posnTrain->size1); // ie 2500 - 2050 = 450
    for (i=0;i<posnTrain->size1;i++){
        sum = 0;
        for (j=0;j<winWidth;j++)
            sum += posn[offset + i*winWidth + winWidth];
        gsl_matrix_set(posnTrain,i,0,sum/winWidth);
    }

    offset += rTrainLast*winWidth;
    for (i=0;i<posnTest->size1;i++){
        sum = 0;
        for (j=0;j<winWidth;j++)
            sum += posn[offset + i*winWidth + winWidth];
        gsl_matrix_set(posnTest,i,0,sum/winWidth);
    }

    gV *A = gsl_vector_alloc(frTrain[0]->size2);

// A = inv(X' X)X' Y

    free(t);
    free(posn);
    gsl_matrix_free(N);
    gsl_matrix_free(mu_true);
    gsl_matrix_free(FR);
    gsl_matrix_free(frTrain[0]);
    gsl_matrix_free(frTest[0]);
    gsl_matrix_free(posnTrain);
    gsl_matrix_free(posnTest);
    gsl_vector_free(A);

    gM *m = gsl_matrix_alloc(3,4);
    double s = 1;
    for (i=0;i<3;i++)
        for (j=0;j<4;j++)
            gsl_matrix_set(m,i,j,s++);

    for (i=0;i<3;i++){
        for (j=0;j<4;j++)
            printf("%i\t", (int)gsl_matrix_get(m,i,j));
        printf("\n");
    }

    gsl_matrix_view M;
    M = gsl_matrix_submatrix(m,0,1,2,1);
    printf("nRows = %i\n", (int)(&(M.matrix))->size1);
    printf("nCols = %i\n", (int)(&(M.matrix))->size2);

    gsl_matrix *M2 = &(M.matrix);
    for (i=0;i< M2->size1;i++){
        for (j=0;j<M2->size2;j++)
            printf("%i\t", (int)gsl_matrix_get(M2,i,j));
        printf("\n");
    }

    gsl_matrix_free(m);

    return(0);
}

