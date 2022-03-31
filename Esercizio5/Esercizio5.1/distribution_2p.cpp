#include <cmath>
#include <vector>
#include "distribution_2p.h"

using namespace std;

Distribution_2p::Distribution_2p(double a0){
	_a0=a0;
}

Distribution_2p::~Distribution_2p(){}

double Distribution_2p::eval(vector<double> v){
	double r=sqrt(v.at(0)*v.at(0) + v.at(1)*v.at(1) + v.at(2)*v.at(2));
	return (1./(32.*M_PI*pow(_a0,5)))*v.at(2)*v.at(2)*exp(-1*r/_a0); 	//|Ψ2,1,0(x,y,z)|^2
}
