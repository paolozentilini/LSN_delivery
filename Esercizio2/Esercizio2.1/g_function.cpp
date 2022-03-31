#include "g_function.h"
#include <cmath>

G_Function::G_Function() {}
G_Function::~G_Function() {}

double G_Function::eval(double x){
	return M_PI*0.5*cos(M_PI*0.5*x)/(M_PI*0.5*cos(M_PI*0.5*x));
}
