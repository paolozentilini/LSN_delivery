#ifndef __Distribution_1s_h__
#define __Distribution_1s_h__

#include "distribution.h"
#include <vector>

class Distribution_1s: public Distribution{	//|Ψ1,0,0(x,y,z)|^2

  protected:
    double _a0;  //Raggio di Bohr
  public:
    Distribution_1s(double a0);
    ~Distribution_1s();
    double eval(std::vector<double>);
};

#endif
