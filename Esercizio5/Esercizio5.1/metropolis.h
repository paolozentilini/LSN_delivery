#ifndef __Metropolis__h
#define __Metropolis__h

#include "random.h"
#include "distribution.h"
#include <vector>

class Metropolis{

  protected:
    Distribution* _p;                        //Per semplicità in questa esercitazione verrà presa simmmetrica.
    Random* _rnd;
    std::vector<double> _current_position;  //Vettore di numeri random che viene aggiornato ogni step.
    double _delta;                          //Larghezza intervallo in cui viene estratto il n° random (caso Gauss: dev_std).
  public:
    //Costruttori e distruttore:
    Metropolis(Random* rnd, Distribution* p, std::vector<double> starting_point, double delta);
    ~Metropolis();
    //Metodi:
    std::vector<double> single_step(); //Singolo step Metropolis pensato per una T(x|x_n) simmetrica.
    double acceptance(); //Metodo che occorre al fine di equilibrare la larghezza dell'intervallo secondo la regola del 50%.
    void equilibrate(unsigned int n); //Metodo che implementa n step Metropolis a vuoto, per equilibrare il sistema.
};

#endif
