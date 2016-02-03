#include<iostream>
#include<fstream>
#include <cmath>
#include <time.h>
#include "particle.h"
#include "myFuncs.h"
using namespace std;

#define BYTES 8

int main(int argc, char *argv[]){
    int i,j,k,p,b;
    int nNeurons   = 79;
    int nParticles = 100; //default value
    int blockSize  = 100;
	double dt      = 0.001;
	int length;
	int nIterations;
	bool flag;
	double v;

	if (argc  > 1)
		nParticles = atoi(argv[1]);
			
	
	unsigned int *nFirings    = new unsigned int[nNeurons];
	bool         *nFiringsTmp = new bool[nNeurons];
	int *resampleIndices      = new int[nParticles];
	double *cumWeights        = new double[nParticles];
	double *xEstimate;
	double norm_sum;
	
	struct {double s1, s2, theta_p1, theta_p2;} sd;
	
	sd.s1       = 0.15;
	sd.s2       = 0.05;
	sd.theta_p1 = 0.001;
	sd.theta_p2 = 0.001;
	srand ( time(NULL) );


	// Create Particles
	particle *pNew = new particle[nParticles];
	particle *pOld = new particle[nParticles];
	particle *pTmp = new particle[nParticles];

	// Open Files
	ifstream fid_s("kinematics.bin",ios::binary);
	ifstream fid_N("spikes.bin"    ,ios::binary);
	ifstream fid_I("init.bin"      ,ios::binary);
	ofstream fid_O("out.bin"       ,ios::binary);
		
	// Define Length
	fid_s.seekg(0,ios::end);
	length = fid_s.tellg()/BYTES;
	nIterations = length / blockSize;
	xEstimate = new double[nIterations];	

	// Initialize Particles		
	for (i=0;i<nParticles;i++){
		pNew[i].init(nNeurons);
		pOld[i].init(nNeurons);
		pTmp[i].init(nNeurons);
	}

	fid_I.read((char *)pOld[0].theta_p , nNeurons*sizeof(double));
	fid_I.read((char *)pOld[0].gamma   , nNeurons*sizeof(double));
	fid_I.read((char *)pOld[0].beta    , nNeurons*sizeof(double));
	fid_I.read((char *)&(pOld[0].s)    ,          sizeof(double));

	pNew[0] = pOld[0];
	for (i=1;i<nParticles;i++){
		pOld[i] = pOld[0];
		pNew[i] = pOld[0];
	}
	
	///////////////////////////////////////////////////////
	///// MASTER LOOP
	///////////////////////////////////////////////////////
	for (k=0 ; k<nIterations ; k++){
		cout << "block number " << k << " of " << nIterations << endl;
	
	///////////////////////////////////////////////////////
	///// Step 1: Draw Random Samples
	///////////////////////////////////////////////////////
	for (p=0 ; p<nParticles ; p++){
		flag = true;
		while (flag) {
			v = pOld[p].s + sd.s1*randn();
			if (v>0 && v<M_PI){
				pNew[p].s = v;
				flag = false;
			}
		}

		for (j=0; j<nNeurons; j++){
		flag = true;
		while (flag) {
			v = pOld[p].theta_p[j] + sd.theta_p1*randn();
			if (v>0 && v<M_PI){
				pNew[p].theta_p[j] = v;
				flag = false;
			}
		}
		}

	}

	///////////////////////////////////////////////////////
	///// Step 2: Compute Likelihoods
	///////////////////////////////////////////////////////
	norm_sum = 0;
	for (j=0;j<nNeurons;j++)
		nFirings[j] = 0;

	for (b=0;b<blockSize;b++){
		fid_N.read((char *)nFiringsTmp , nNeurons*sizeof(bool));
		for (j=0;j<nNeurons;j++)
			nFirings[j] += nFiringsTmp[j];
	}
	
	
	for (p=0 ; p<nParticles ; p++){
		pNew[p].likelihood = compute_likelihood(pNew[p],nFirings,blockSize,nNeurons,dt);
		pNew[p].g = pOld[p].w * pNew[p].likelihood;
		norm_sum += pNew[p].g;
	}

	for (p=0 ; p<nParticles ; p++)
		pNew[p].g /= norm_sum;

	///////////////////////////////////////////////////////
	///// Step 3: Resample
	///////////////////////////////////////////////////////
	cumWeights[0] = pNew[0].g;
	for (p=1 ; p<nParticles ; p++)
		cumWeights[p] = cumWeights[p-1] + pNew[p].g;

	resample(resampleIndices, cumWeights, nParticles);
	
	for (p=0 ; p<nParticles ; p++)
		pTmp[p] = pNew[p];
		
	for (p=0 ; p<nParticles ; p++){
		j = resampleIndices[p];
		pNew[p].resample(pTmp[j]);
	}

	///////////////////////////////////////////////////////
	///// Step 4: Draw Random Samples
	///////////////////////////////////////////////////////
	for (p=0 ; p<nParticles ; p++){
		flag = true;
		while (flag) {
			v = pNew[p].s + sd.s2*randn();
			if (v>0 && v<M_PI){
				pNew[p].s = v;
				flag = false;
			}
		}

		for (j=0; j<nNeurons; j++){
		flag = true;
		while (flag) {
			v = pNew[p].theta_p[j] + sd.theta_p2*randn();
			if (v>0 && v<M_PI){
				pNew[p].theta_p[j] = v;
				flag = false;
			}
		}
		}

	}

	///////////////////////////////////////////////////////
	///// Step 5: Compute Resampled Likelihoods
	///////////////////////////////////////////////////////
	norm_sum = 0;
	for (p=0 ; p<nParticles ; p++){
		pNew[p].resampled_likelihood = compute_likelihood(pNew[p],nFirings,blockSize,nNeurons,dt);
		pNew[p].w                    = pNew[p].resampled_likelihood / pNew[p].likelihood;
		norm_sum += pNew[p].w;
	}	

	for (p=0 ; p<nParticles ; p++)
		pNew[p].w /= norm_sum;
		
	///////////////////////////////////////////////////////
	///// Step 6: Compute Weighted Estimates
	///////////////////////////////////////////////////////
	xEstimate[k] = 0;
	for (p=0 ; p<nParticles ; p++)
		xEstimate[k] += pNew[p].s * pNew[p].w;	

	for (p=0 ; p<nParticles ; p++)
		pOld[p] = pNew[p];	
	
	}

	fid_O.write((char *)xEstimate , nIterations*sizeof(double));

	// Clean Up
	fid_s.close();
	fid_N.close();
	fid_I.close();
	fid_O.close();
	delete [] pNew;
	delete [] pOld;
	delete [] pTmp;
	delete [] nFirings;
	delete [] nFiringsTmp;
	delete [] resampleIndices;
	delete [] cumWeights;
	delete [] xEstimate;
	return 0;
}

