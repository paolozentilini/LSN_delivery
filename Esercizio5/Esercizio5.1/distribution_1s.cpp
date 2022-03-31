#include <cmath>
#include <vector>
#include "distribution_1s.h"

using namespace std;

Distribution_1s::Distribution_1s(double a0){
	_a0=a0;
}

Distribution_1s::~Distribution_1s(){}

double Distribution_1s::eval(vector<double> v){
	double r=sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	return  (1./(M_PI*pow(_a0,3)))*exp(-2*r/_a0);	//|Ψ1,0,0(x,y,z)|^2
}
