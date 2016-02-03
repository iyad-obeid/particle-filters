#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "oMatrixPrimitives.h"
#include "createData.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <dirent.h>

double getenvNumeric(char *env){
  char *str;
  double val;
  if ((str=getenv(env))==NULL)
    val = -1;
  else
    val = atof(str);
  return(val);
}

double nSamples;
double nNeurons;
gsl_rng *r;
char   *dataType;
char outDir[] = "../results/";
char   fNameIn[50];
double dt = 2e-3;
double dataParam[5];
double tMax;

int main(){
    r = gsl_rng_alloc(gsl_rng_mt19937);
    o_init_rng(r);

  nNeurons = getenvNumeric("O_nNeurons");
  dataType = getenv("O_dataStyle");
    dataParam[0] = getenvNumeric("O_dataParam0");
    dataParam[1] = getenvNumeric("O_dataParam1");
    dataParam[2] = getenvNumeric("O_dataParam2");
    dataParam[3] = getenvNumeric("O_dataParam3");
    dataParam[4] = getenvNumeric("O_dataParam4");
  tMax         = getenvNumeric("O_tMax");

  double filenum      = getenvNumeric("O_dataParam3");
  sprintf(fNameIn,"%sd%05iinp.bin",outDir,(int)filenum);
  printf("%s\n",fNameIn);
  createData(fNameIn);

  return (0);
}
