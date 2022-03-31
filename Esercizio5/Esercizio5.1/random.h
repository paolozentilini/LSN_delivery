#ifndef __Random__h_
#define __Random__h_

class Random {

	public:

  	virtual void set_random(int*, int, int)=0;
  	virtual void save_seed()=0;
		virtual void init()=0;
  	virtual double get_uniform()=0;	//Restituisce un numero casuale uniformemente distribuito tra 0 e 1 nelle classi concrete
  	virtual double get_random_number()=0;	//Restituisce un numero random della distribuzione oggetto corrispondente
		virtual void set_parameters(double,double)=0; //Setta i parametri della distribuzione 
};

#endif // __Random__
