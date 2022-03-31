#ifndef __Integral_1__h
#define __Integral_1__h

#include "experiment.h"

//Classe concreta figlia della classe virtuale Experiment. Questa classe "compie" una singola misura dell'integrale
// per il calcolo della varianza di r.

class Integral_1 : public Experiment{
  public:
    //Costruttore e distruttore:
    Integral_1();
    ~Integral_1();
    //Questa classe eredita da Experiment il metodo virtuale single_measurement e viene implementato qui al fine di rendere il codice
    //il più plastico possibile. Infatti è sufficiente richiamare nel data_blocking il metodo virtuale di Experiment senza dover ogni volta
    //ricostruire il blocking su misura per il singolo integrale.
    //In questo caso single_measurement restituisce la singola misura dell'integranda ( r - <r> )^2 dove <r>=0.5.
    double single_measurement(double);

};

#endif
