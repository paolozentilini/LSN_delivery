#ifndef __Buffon_Experiment_
#define __Buffon_Experiment_

#include <vector>
#include "experiment.h"
#include "random.h"

using namespace std;

class Buffon_Experiment : public Experiment{

  protected:
    Random* _rnd;
    double _length;
    double _pi;
    double _distance;
    unsigned int _n_throws;

  public:

    Buffon_Experiment(Random* rnd, double length, double distance, unsigned int n_throws);
    ~Buffon_Experiment();

    vector<double> needle();
    double single_measurement();
};

#endif
