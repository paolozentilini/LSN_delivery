#include <cmath>
#include "integral.h"

Integral::Integral(Metropolis* mt, double mu, double sigma){
  _mt=mt;
  _mu = mu;
  _sigma = sigma;
}

Integral::~Integral(){}

double Integral::single_evaluation(){
  double x = _mt->single_step();
  double potential = x*x*(x*x - 2.5);
  double kinetic = -(x*x+ _mu*_mu -_sigma*_sigma -2*x*_mu*tanh((_mu*x)/(_sigma*_sigma)) )/(2.*pow(_sigma,4));
  return potential+kinetic;
}
