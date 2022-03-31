#ifndef __Integrals__h
#define __Integrals__h

#include "metropolis.h"
#include <vector>

using namespace std;

class Integral{

  protected:
    Metropolis* _mt;
  public:
    //Costruttori e distruttore:
    Integral(Metropolis* mt);
    ~Integral();

    double single_evaluation(); //Valutazione dell'integrale <r> e scrittura su file
};

#endif
