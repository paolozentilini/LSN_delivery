#ifndef __PlainVanilla__h_
#define __PlainVanilla__h_

#include "wiener_process.h"
#include <string>

using namespace std;

class PlainVanilla {
  private:
    double _strike;         //strike price
    double _r;              //risk-free interest rate
    string _option_type;     //PlainVanilla option type (Put or Call)
    WienerProcess* _wiener_process; //Stochastic process of the asset price
  public:
    PlainVanilla(double strike, double r, WienerProcess* wiener_process, string option_type);
    ~PlainVanilla();
    double pay_off();
};

#endif
