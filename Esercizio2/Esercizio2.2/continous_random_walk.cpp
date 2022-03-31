#include "continous_random_walk.h"
#include <cmath>

using namespace std;

ContinousRW::ContinousRW(double a,Random* rnd){
  _a=a;
  _rnd = rnd;
  _rnd->Init();
}

ContinousRW::~ContinousRW(){
  _rnd->SaveSeed();
}

double ContinousRW::single_measurement(unsigned int n_steps){
  double x=0;
  double y=0;
  double z=0;
  double phi;
  double theta;
  for(unsigned int j=0; j<n_steps; j++){
    phi = _rnd->get_uniform(0,2*M_PI);			//Generating random angles phi in [0,2pi] and theta between 0 and pi
    theta = acos(1-2*_rnd->get_uniform());
    x+=_a*sin(theta)*cos(phi);              //Ad ogni step aggiorno la direzione dello spostamento unitario nello spazio 3D
    y+=_a*sin(theta)*sin(phi);
    z+=_a*cos(theta);
  }
  double squared_module = x*x + y*y + z*z;
  return squared_module;                   //Ritorna il modulo quadro
}
