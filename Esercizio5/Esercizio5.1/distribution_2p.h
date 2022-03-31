#ifndef __Distribution_2p_h__
#define __Distribution_2p_h__

#include "distribution.h"
#include <vector>

class Distribution_2p: public Distribution{	//|Ψ2,1,0(x,y,z)|^2

  protected:
    double _a0;  //Raggio di Bohr
  public:
    Distribution_2p(double a0);
    ~Distribution_2p();
    double eval(std::vector<double>);
};

#endif
