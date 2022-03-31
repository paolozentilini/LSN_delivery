#ifndef __Integrals__h
#define __Integrals__h

#include "metropolis.h"
#include <vector>

using namespace std;

class Integral{

  protected:
    Metropolis* _mt;
    double _sigma,_mu;
  public:
    //Costruttori e distruttore:
    Integral(Metropolis* mt, double mu, double sigma);
    ~Integral();

    double single_evaluation(); //Valutazione dell'integrale 
};

#endif
