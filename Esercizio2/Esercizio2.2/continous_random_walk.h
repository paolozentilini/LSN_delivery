#ifndef __ContinousRW__h
#define __ContinousRW__h

#include "random_walk.h"
#include "random.h"

class ContinousRW : public RandomWalk{

  protected:
    double _a;      //modulo della lunghezza di un singolo step 
    Random* _rnd;

  public:
    //Costruttori e distruttore:
    ContinousRW(double a, Random* rnd);
    ~ContinousRW();
    //Metodo ereditato dalla classe astratta RandomWalk e implementato qui.
    double single_measurement(unsigned int n_steps);

};

#endif
