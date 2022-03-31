#ifndef __Integral__h
#define __Integral__h

#include "experiment.h"

//Classe concreta figlia della classe virtuale Experiment. Questa classe "compie" una singola misura dell'integrale
// per il calcolo della media di r.
class Integral : public Experiment{

  public:
    //Costruttori e distruttore:
    Integral();
    ~Integral();
    //Questa classe eredita da Experiment il metodo virtuale single_measurement e viene implementato qui al fine di rendere il codice
    //il più plastico possibile. Infatti è sufficiente richiamare nel data_blocking il metodo virtuale di Experiment senza dover ogni volta
    //ricostruire il blocking su misura per il singolo integrale.
    //In questo caso single_measurement restituisce la singola misura dell'integranda r.
    double single_measurement(double);

};

#endif
