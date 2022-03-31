#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "gaussian.h"

using namespace std;

Gaussian :: Gaussian(double mean, double sigma){
	_mean = mean;
	_sigma = sigma;
}

Gaussian :: Gaussian(){}

Gaussian :: ~Gaussian(){}

void Gaussian :: set_random(int * s, int p1, int p2){
	m1 = 502;
	m2 = 1521;
	m3 = 4071;
	m4 = 2107;
	l1 = s[0]%4096;
	l2 = s[1]%4096;
	l3 = s[2]%4096;
	l4 = s[3]%4096;
	l4 = 2*(l4/2)+1;
	n1 = 0;
	n2 = 0;
	n3 = p1;
	n4 = p2;
	return;
}

void Gaussian :: save_seed(){
	ofstream WriteSeed;
	WriteSeed.open("seed.out");
	if (WriteSeed.is_open()){
		WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
	} else cerr << "PROBLEM: Unable to open seed.out" << endl;
		WriteSeed.close();
	return;
}

void Gaussian :: init(){
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
      		Primes >> p1 >> p2 ;
   	} else cerr << "PROBLEM: Unable to open Primes" << endl;
  	 Primes.close();
  	 ifstream input("seed.in");
  	 string property;
  	 if (input.is_open()){
     	 while ( !input.eof() ){
         	input >> property;
        	 if( property == "RANDOMSEED" ){
          	  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
          	  set_random(seed,p1,p2);
        	 }
   	   }
     	 input.close();
   	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

void Gaussian::set_parameters(double mean, double sigma){
	_mean=mean;
	_sigma=sigma;
}

double Gaussian :: get_uniform(){
	const double twom12=0.000244140625;
  int i1,i2,i3,i4;
	double r;
	i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
	i2 = l2*m4 + l3*m3 + l4*m2 + n2;
	i3 = l3*m4 + l4*m3 + n3;
	i4 = l4*m4 + n4;
	l4 = i4%4096;
	i3 = i3 + i4/4096;
	l3 = i3%4096;
	i2 = i2 + i3/4096;
	l2 = i2%4096;
	l1 = (i1 + i2/4096)%4096;
	r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));
	return r;
}

double Gaussian :: get_random_number() {
	double s=get_uniform();
	double t=get_uniform();
	double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
	return _mean + x * _sigma;
}
