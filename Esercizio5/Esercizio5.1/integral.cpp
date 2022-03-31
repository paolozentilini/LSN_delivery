#include <cmath>
#include "integral.h"

Integral::Integral(Metropolis* mt){
  _mt=mt;
}

Integral::~Integral(){}

double Integral::single_evaluation(){
  vector<double> random_pos = _mt->single_step();
  double sum_squared=0;
  for(unsigned int i=0; i < random_pos.size(); i++)
    sum_squared += random_pos.at(i)*random_pos.at(i);
  return sqrt(sum_squared);
}
