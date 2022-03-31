#include "discrete_random_walk.h"
#include "random_walk.h"
#include <cmath>

using namespace std;

DiscreteRW :: DiscreteRW(double a, Random* rnd){
  _a=a;
  _rnd = rnd;
  _rnd->Init();
}

DiscreteRW ::~DiscreteRW(){
  _rnd->SaveSeed();
}

double DiscreteRW ::single_measurement(unsigned int n_steps){
  double x=0;
  double y=0;
  double z=0;
  for(unsigned int j=0; j<n_steps; j++){
    double s = _rnd->get_uniform();
    if(s<1./3){				             	//x coordinate
      if(s<1./6) x+= _a;
      else x+= -1*_a;
    }
    if(s>1./3 && s<2./3){			   		//y coordinate
      if(s<0.5) y+= _a;
      else y+= -1*_a;
    }
    if(s>2./3){				            	//z coordinate
      if(s<5./6) z+= _a;
      else z+= -1*_a;
    }
  }
  double squared_module = x*x + y*y + z*z;
return squared_module;              //Ritorna il modulo quadro dello spostamento
}


/*vector<double> DiscreteRW :: show_RW(unsigned int n_steps){
  double x=0;
  double y=0;
  double z=0;
  vector<double> distances(n_steps,0);

  for(unsigned int j=0; j<n_steps; j++){
    double s = _rnd->get_uniform();
    if(s<1./3){				             	//x coordinate
      if(s<1./6) x+= _a;
      else x+= -1*_a;
    }
    if(s>1./3 && s<2./3){			   		//y coordinate
      if(s<0.5) y+= _a;
      else y+= -1*_a;
    }
    if(s>2./3){				            	//z coordinate
      if(s<5./6) z+= _a;
      else z+= -1*_a;
    }
    double squared_module = x*x+y*y+z*z;
    distances.append(squared_module;)
  }
  return distances;
}*/
