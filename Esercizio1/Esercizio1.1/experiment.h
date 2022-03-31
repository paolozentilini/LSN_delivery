//Classe astratta per il calcolo dei due integrali richiesti nell'esercizio.
//Vedere le implementazioni delle classi figlie per lo scopo reale di questa classe virtuale.

#ifndef __Experiment__h_
#define __Experiment__h_

class Experiment{
    public:
        virtual double single_measurement(double random)=0;
};

#endif
