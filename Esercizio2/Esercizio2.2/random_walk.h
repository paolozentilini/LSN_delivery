#ifndef __RandomWalk__h
#define __RandomWalk__h

//Classe astratta che permette in Statistics di richiamare il metodo single_measurement senza specificare internamente il tipo di RandomWalk
class RandomWalk{ 
  public:
    virtual double single_measurement(unsigned int n_steps)=0;
};
#endif
