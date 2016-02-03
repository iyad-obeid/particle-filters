double randn();
double compute_likelihood(particle &, unsigned int *, int, int, double);
void resample(int *resampleIndices, double *cumWeights, int nParticles);
