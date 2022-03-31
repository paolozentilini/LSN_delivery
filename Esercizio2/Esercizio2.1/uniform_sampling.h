#ifndef __Uniform_Sampling__h
#define __Uniform_Sampling__h

#include "montecarlo.h"
#include "funzionebase.h"
#include "random.h"

class Uniform_Sampling : public MonteCarlo{
  protected:
    double _a;        //Estremi dell'intervallo
    double _b;
    FunzioneBase* _f; //Funzione da valutare nell'integrale
    Random* _rnd;
    unsigned int _n;  //Numero di volte con cui calcolo il valore dell'integranda per usare il metodo della media

  public:
    //Costruttori e distruttore:
    Uniform_Sampling(double a, double b, FunzioneBase* f,Random* rnd,unsigned int n);
    ~Uniform_Sampling();
    //Metodo ereditato dalla classe astratta MonteCarlo e implementato qui.
    double single_measurement();

};

#endif
