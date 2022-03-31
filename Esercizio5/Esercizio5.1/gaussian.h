#ifndef __Gaussian__h_
#define __Gaussian__h_

#include "random.h"

class Gaussian: public Random {

	private:
		int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;
	protected:
    double _mean;
    double _sigma;
	public:
  	//constructors
 		Gaussian(double mean, double sigma);
		Gaussian();
  	//destructor
  	~Gaussian();
 		//methods
  	void set_random(int*, int, int);
  	void save_seed();
		void init();
		void set_parameters(double mean, double sigma);
  	double get_uniform(void);	//Restituisce un numero casuale uniformemente distribuito tra 0 e 1
  	double get_random_number();	//Restituisce un numero gaussiano con varianza sigma e media mean
};

#endif
