#include "discrete_wiener_process.h"
#include "wiener_process.h"
#include <cmath>

using namespace std;

DiscreteWienerProcess :: DiscreteWienerProcess(double S_0, double T, unsigned int n_steps, double mu, double sigma, Random* rnd){
  _S_0 = S_0;
  _T=T;
  _n_steps = n_steps;
  _mu=mu;
  _sigma=sigma;
  _rnd = rnd;
  _rnd->Init();
}

DiscreteWienerProcess ::~DiscreteWienerProcess(){
  _rnd->SaveSeed();
}

double DiscreteWienerProcess :: asset_price(){
  double delta_t = _T/(1.*_n_steps);
  _S_t = _S_0;
  for(unsigned int j=0; j<_n_steps; j++){
    _S_t = _S_t*exp((_mu - 0.5*_sigma*_sigma)*delta_t + _sigma*sqrt(delta_t)*_rnd->get_gaussian(0,1));
  }
  return _S_t;
}

double DiscreteWienerProcess :: get_delivery_time(){
  return _T;
}
