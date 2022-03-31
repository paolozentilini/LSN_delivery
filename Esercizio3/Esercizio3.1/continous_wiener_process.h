#ifndef __ContinousWienerProcess__h
#define __ContinousWienerProcess__h

#include "wiener_process.h"
#include "random.h"

class ContinousWienerProcess : public WienerProcess{

  protected:
    double _S_0;            //asset price at t=0
    double _S_t;            //asset price at time t
    double _T;              //delivery time
    double _mu;             //risk-free interest rate
    double _sigma;          //volatility parameter
    Random* _rnd;

  public:
    //Costruttori e distruttore:
    ContinousWienerProcess( double S_0, double T, double mu, double sigma, Random* rnd );
    ~ContinousWienerProcess();
    //Metodo che calcola il prezzo dell'asset al tempo di delivery.
    double asset_price();
    //Metodo che restituisce il tempo di delivery
    double get_delivery_time();

};

#endif
