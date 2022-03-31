#include "montecarlo.h"
#include "uniform_sampling.h"
#include "funzionebase.h"
#include "random.h"

using namespace std;

Uniform_Sampling:: Uniform_Sampling(double a, double b, FunzioneBase* f, Random* rnd ,unsigned int n){
  _a=a;
  _b=b;
  _f=f;
  _n=n;
  _rnd = rnd;
  _rnd->Init();
}

Uniform_Sampling::~Uniform_Sampling(){
  _rnd->SaveSeed();
}

//Uso il metodo della media per calcolare l'integrale tramite p.d.f. uniforme
double Uniform_Sampling::single_measurement(){
  double sum=0;
  for (unsigned int i=0; i<_n; i++) sum += _f->eval(_rnd->get_uniform(_a,_b));
  return (_b-_a)*sum/_n;
}
