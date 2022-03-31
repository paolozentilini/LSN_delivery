#include <cmath>
#include <vector>
#include "distribution_trial.h"

using namespace std;

Distribution_trial::Distribution_trial(double mu, double sigma){
	_sigma=sigma;
	_mu = mu;
}

Distribution_trial::~Distribution_trial(){}

double Distribution_trial::eval(double x){
	double psi = exp(-(x-_mu)*(x-_mu)/(2.*_sigma*_sigma)) +  exp(-(x+_mu)*(x+_mu)/(2.*_sigma*_sigma));
	return psi*psi;
}
