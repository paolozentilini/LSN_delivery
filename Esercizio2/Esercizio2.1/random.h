#ifndef __Random__
#define __Random__

class Random {

	private:
		int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;
	public:
  	//constructors
 		Random();
  	//destructor
  	~Random();
 		//methods
  	void SetRandom(int*, int, int);
  	void SaveSeed();
		void Init();
  	double get_uniform(void);	//Restituisce un numero casuale uniformemente distribuito tra 0 e 1
  	double get_uniform(double min, double max);	//Restituisce un numero casuale uniformemente distribuito tra min e max
  	double get_gaussian(double mean, double sigma);	//Restituisce un numero gaussiano con varianza sigma e media mean
		//Uso il metodo della trasformata per generare numeri casuali con distribuzione cos:
		double get_cos1(); //Restituisce una variabile casuale estratta da una distribuzione pi/2(1-x) tra 0 e 1
		double get_cos2(); //Restituisce una variabile casuale estratta da una distribuzione pi/2*cos(pi/2*x) tra 0 e 1

};

#endif // __Random__
