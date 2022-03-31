#ifndef __Distribution_trial_h__
#define __Distribution_trial_h__

#include "distribution.h"
#include <vector>

class Distribution_trial: public Distribution{	//|Ψ|^2

  protected:
    double _mu;
    double _sigma;
  public:
    Distribution_trial(double mu, double sigma);
    ~Distribution_trial();
    double eval(double x);
};

#endif
