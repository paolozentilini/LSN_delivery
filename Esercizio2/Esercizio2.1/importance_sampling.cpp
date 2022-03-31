#include "montecarlo.h"
#include "importance_sampling.h"
#include <cmath>

using namespace std;

Importance_Sampling :: Importance_Sampling(double a, double b, FunzioneBase* f, Random* rnd, unsigned int n, unsigned int label){
  _a=a;
  _b=b;
  _f=f;
  _n=n;
  _label=label;
  _rnd = rnd;
  _rnd->Init();
}

Importance_Sampling ::~Importance_Sampling(){
  _rnd->SaveSeed();
}

double Importance_Sampling::single_measurement(){
	double sum=0;
  double appo;
  if(_label==0){
	  for(unsigned int i=0; i<_n; i++){
      appo = _rnd->get_cos2();
      sum += _f->eval(appo);
    }
  }else{
    for(unsigned int i=0; i<_n; i++){
      appo = _rnd->get_cos1();
      sum += _f->eval(appo);
   }
  }
  return sum/_n;
}
/*
double Importance_Sampling::hit_or_miss(double m) {
	unsigned int n_hit = 0;
	for (unsigned int i=0; i< _n; i++) {
		double z = _rnd->get_uniform(a,b);
		double y = _rnd->get_uniform(0,m);
		if (f->Eval(z) > y) n_hit++;
	}
	return (b-a)*(m*n_hit)/_n;
}
*/
