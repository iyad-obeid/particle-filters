#include <cmath>
#include <iostream>
#include "particle.h"

double randn(){
	double u1,u2;
	u1 = (double)rand()/(double)RAND_MAX;
	u2 = (double)rand()/(double)RAND_MAX;
	return sqrt(-2.0 * log(u1) ) * cos(2.0*M_PI*u2);
}

double compute_likelihood(particle &part, unsigned int *N, int blockSize, int nNeurons, double dt){
	int j;
	double lambda, neuronal_likelihood, f;
	f = 1;
	
	for (j=0;j<nNeurons;j++){
		lambda = part.beta[j] + part.gamma[j]*cos(2*(part.s - part.theta_p[j]));
		neuronal_likelihood = pow((lambda*dt),N[j]) * exp(-lambda*dt*blockSize);
		f *= neuronal_likelihood;
	}
	return f;
}

void resample(int *resampleIndices, double *cumWeights, int nParticles){
	double randVal;
	int p,j;
	bool flag;
	
	for (p=0 ; p<nParticles ; p++){
		randVal = (double)rand()/(double)RAND_MAX;
		flag = true;
		j = 0;
		while (flag){
			if (randVal < cumWeights[j]){
				flag = false;
				resampleIndices[p] = j;
			}
			else
				j++;
		}
	}
}
