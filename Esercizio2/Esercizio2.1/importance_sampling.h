#ifndef __Importance_Sampling__h
#define __Importance_Sampling__h

#include "montecarlo.h"
#include "funzionebase.h"
#include "random.h"

class Importance_Sampling : public MonteCarlo{

  protected:
    double _a;           //Estremi dell'intervallo
    double _b;
    FunzioneBase* _f;    //Funzione da valutare nell'integrale
    Random* _rnd;
    unsigned int _n;     //Numero di volte con cui calcolo il valore dell'integranda per usare il metodo della media
    unsigned int _label; //Flag utile per fare i due diversi calcoli per l'importance_sampling

  public:
    //Costruttori e distruttore:
    Importance_Sampling(double a, double b, FunzioneBase* f, Random* rnd, unsigned int n, unsigned int label);
    ~Importance_Sampling();
    //Metodo ereditato dalla classe astratta MonteCarlo e implementato qui.
    double single_measurement();

};

#endif
