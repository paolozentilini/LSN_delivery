#ifndef __Uniform__h_
#define __Uniform__h_

#include "random.h"

class Uniform: public Random {

	private:
		int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;
	protected:
    double _min;
    double _max;
	public:
  	//constructors
 		Uniform(double min, double max);
		Uniform();
  	//destructor
  	~Uniform();
 		//methods
  	void set_random(int*, int, int);
  	void save_seed();
		void init();
		void set_parameters(double center_of_interval, double delta);
  	double get_uniform(void);	//Restituisce un numero casuale uniformemente distribuito tra 0 e 1
  	double get_random_number();	//Restituisce un numero uniforme tra min e max
};

#endif
