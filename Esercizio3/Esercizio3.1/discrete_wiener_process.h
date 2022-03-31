#ifndef __DiscreteWienerProcess__h
#define __DiscreteWienerProcess__h

#include "wiener_process.h"
#include "random.h"

class DiscreteWienerProcess : public WienerProcess{

  protected:
    double _S_0;            //asset price at t=0
    double _S_t;            //asset price at time t
    double _T;              //delivery time
    double _mu;             //risk-free interest rate
    double _sigma;          //volatility parameter
    unsigned int _n_steps;  //number of intervals
    Random* _rnd;

  public:
    //Costruttori e distruttore:
    DiscreteWienerProcess(double S_0, double T, unsigned int n_steps, double mu, double sigma, Random* rnd);
    ~DiscreteWienerProcess();
    //Metodo che calcola il prezzo dell'asset al tempo di delivery.
    double asset_price();
    //Metodo che restituisce il tempo di delivery
    double get_delivery_time();
};

#endif
