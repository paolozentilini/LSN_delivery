#ifndef __Statistics__
#define __Statistics__

#include <vector>
#include "random.h"
#include "experiment.h"

class Statistics {

  protected:

    unsigned int _n_blocks;
    unsigned int _n_measures;
    unsigned int _n_per_blocks;
    Random *_random;
    Experiment *_experiment;

  public:

    Statistics(unsigned int n_blocks, unsigned int n_measures, Random*, Experiment*);
    Statistics(Random*);
    ~Statistics();

    double mean_value(std::vector<double>, unsigned int);
    double dev_std(std::vector<double>, unsigned int);
    double chi_squared(std::vector<int>, double);

    void chi_squared_test(unsigned int n, unsigned int throws_per_block, const char*);
    void data_blocking(const char*);

};

#endif
