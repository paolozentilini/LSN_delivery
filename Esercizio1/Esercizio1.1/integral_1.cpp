#include <iostream>
#include "integral_1.h"

using namespace std;

Integral_1:: Integral_1(){}

Integral_1::~Integral_1(){}

double Integral_1::single_measurement(double random){
  return (random-0.5)*(random-0.5);
}
