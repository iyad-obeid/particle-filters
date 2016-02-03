#include "particle.h"
#include <iostream>
using namespace std;

particle& particle::operator=(const particle &op2){
  int i;
  if (nNeurons != op2.nNeurons)
    cout << "Number of Neurons must be the same!!" << endl;

  s = op2.s;
  g = op2.g;
  w = op2.w;

  likelihood           = op2.likelihood;
  resampled_likelihood = op2.resampled_likelihood;

  for (i=0;i<nNeurons;i++){
  	gamma[i]   = op2.gamma[i];
  	theta_p[i] = op2.theta_p[i];
  	beta[i]    = op2.beta[i];
  }
  return *this;
}

void particle::resample(const particle &op2){
	int i;
	s = op2.s;
	for (i=0 ; i<nNeurons; i++)
		theta_p[i] = op2.theta_p[i];
}

void particle::init(int n){
  	nNeurons = n;
	s = 0;
	g = 0;
	w = 1;

	likelihood = 0;
	resampled_likelihood = 0;

	gamma 	= new double[nNeurons];
	theta_p = new double[nNeurons];
	beta    = new double[nNeurons];
	for (int i = 0; i<nNeurons; i++){
		gamma  [i] = 0;
		theta_p[i] = 0;
		beta   [i] = 0;
	}
}

particle::~particle(){
  delete [] gamma;
  delete [] theta_p;
  delete [] beta;
}

void particle::query(){
	int i;
	for (i=0;i<nNeurons;i++)
		cout << "gamma[" << i << "]=" << gamma[i] << endl;
	cout << endl;
}
