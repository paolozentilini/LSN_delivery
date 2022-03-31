//Classe astratta per il calcolo MonteCarlo

#ifndef __MonteCarlo__h_
#define __MonteCarlo__h_

class MonteCarlo{
    public:
        virtual double single_measurement()=0; //Metodo astratto che calcola il valore dell'integrale
};

#endif
