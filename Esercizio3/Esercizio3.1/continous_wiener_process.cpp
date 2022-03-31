#include "continous_wiener_process.h"
#include "wiener_process.h"
#include <cmath>

using namespace std;

ContinousWienerProcess :: ContinousWienerProcess(double S_0, double T, double mu, double sigma, Random* rnd){
  _S_0 = S_0;
  _T=T;
  _mu=mu;
  _sigma=sigma;
  _rnd = rnd;
  _rnd->Init();
}

ContinousWienerProcess ::~ContinousWienerProcess(){
  _rnd->SaveSeed();
}

double ContinousWienerProcess :: asset_price(){
  _S_t = _S_0*exp((_mu - 0.5*_sigma*_sigma)*_T + _sigma*sqrt(_T)*_rnd->get_gaussian(0,1));
  return _S_t;
}

double ContinousWienerProcess :: get_delivery_time(){
  return _T;
}
