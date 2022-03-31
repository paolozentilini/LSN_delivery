#include "f_function.h"
#include <cmath>

F_Function::F_Function() {}
F_Function::~F_Function() {}

double F_Function::eval(double x){
	return M_PI*cos(M_PI*0.5*x)/(4*(1-x));
}
