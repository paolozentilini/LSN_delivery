#ifndef __WienerProcess__h
#define __WienerProcess__h

//Classe astratta richiamata in statistic.h per consentire di generalizzare il pricing dell'opzione 
class WienerProcess {

  public:
    virtual double asset_price()=0;
    virtual double get_delivery_time()=0;
};
#endif
