#include "coseno.h"
#include <cmath>

Coseno::Coseno()  {}
Coseno::~Coseno() {}

double Coseno::eval(double x){
	return 0.5*M_PI*cos(0.5*M_PI*x);
}
