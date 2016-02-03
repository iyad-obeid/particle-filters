class particle {
private:
    int nNeurons;
public:
	double *gamma;
	double *theta_p;
	double *beta;
	double s;
	double g;
	double w;
	double likelihood;
	double resampled_likelihood;
	void resample(const particle&);
    particle(){};
	~particle();
    void init(int nNeurons);
    particle &operator=(const particle&);
	void query();
};
