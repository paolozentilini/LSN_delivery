#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "uniform.h"

using namespace std;

Uniform :: Uniform(double min, double max){
	_min=min;
	_max=max;
}

Uniform :: Uniform(){}

Uniform :: ~Uniform(){}

void Uniform :: set_random(int * s, int p1, int p2){
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

void Uniform :: save_seed(){
	ofstream WriteSeed;
	WriteSeed.open("seed.out");
	if (WriteSeed.is_open()){
		WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
	} else cerr << "PROBLEM: Unable to open seed.out" << endl;
		WriteSeed.close();
	return;
}

void Uniform :: init(){
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

void Uniform::set_parameters(double center_of_interval, double delta){
	_min=center_of_interval-delta;
	_max=center_of_interval+delta;
}

double Uniform :: get_uniform(void){
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

double Uniform :: get_random_number(){
	return _min + (_max-_min)*get_uniform();
}
