#ifndef _SimulatedAnnealing_h__
#define _SimulatedAnnealing_h__

#include <cmath>
#include <iostream>
#include <vector>
#include "random.h"
#include "data_structures.h"

using namespace std;

class SimulatedAnnealing{

  protected:

    Random* _rnd;
    unsigned int _n_cities;
    Chromosome* _chromosome1;
    Chromosome* _chromosome2;
    double _mutation_probability;

  public:

    SimulatedAnnealing(Random* rnd, unsigned int n_cities);
    ~SimulatedAnnealing();

    void start_on_circumference();
    void start_on_plane();
    void mutation1();
    void mutation2();
    void mutation3();
    void mutation4();
    void evolution_step(double beta);

    double cost();
    Chromosome get_chromosome();

};

  vector<Gene> swap(vector<Gene>, int, int);
  vector<Gene> permutation(vector<Gene>, unsigned int);
  vector<Gene> permutation(vector<Gene>, unsigned int, unsigned int);
  vector<Gene> inversion(vector<Gene>, unsigned int, unsigned int);

#endif
