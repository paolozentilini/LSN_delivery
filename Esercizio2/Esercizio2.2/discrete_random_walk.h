#ifndef __DiscreteRW__h
#define __DiscreteRW__h

#include "random_walk.h"
#include "random.h"

class DiscreteRW : public RandomWalk{

  protected:
    double _a;    //passo reticolare
    Random* _rnd;

  public:
    //Costruttori e distruttore:
    DiscreteRW(double a, Random* rnd);
    ~DiscreteRW();
    //Metodo ereditato dalla classe astratta RandomWalk e implementato qui.
    double single_measurement(unsigned int n_steps);

};

#endif
