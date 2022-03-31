#ifndef _GeneticAlgorithm_h__
#define _GeneticAlgorithm_h__

#include <cmath>
#include <iostream>
#include <vector>
#include "random.h"
#include "data_structures.h"

using namespace std;

class GeneticAlgorithm{

  protected:

    Random* _rnd;
    unsigned int _n_cities;
    unsigned int _n_chromosomes;
    vector<Chromosome> _population, _new_population;
    Chromosome _chromosome1, _chromosome2;
    int _position_chromosome1, _position_chromosome2;
    vector<Gene> _offspring1, _offspring2;
    double _mutation_probability, _crossover_probability;

  public:

    GeneticAlgorithm(Random* rnd, unsigned int n_cities, unsigned int n_chromosomes, double p_m, double p_c);
    ~GeneticAlgorithm();

    unsigned int get_n_chromosomes(){return _n_chromosomes; }
    vector<Chromosome> get_population(){ return _population; }
    vector<int> get_best_chromosome_sequence();
    void set_best_chromosome(vector<int>);
    double best_fit();
    double mean_fit();
    void clear_population();

    void start(vector<double>, vector<double>, vector<int>);
    void select();
    void crossover();
    void mutation1();
    void mutation2();
    void mutation3();
    void mutation4();
    void recreate_population();
    void population_update();

    void evolution();

};

bool my_cmp( Chromosome&, Chromosome&);
bool in( const Gene, vector<Gene>);
vector<Gene> swap(vector<Gene>, int, int);
vector<Gene> permutation(vector<Gene>, unsigned int);
vector<Gene> permutation(vector<Gene>, unsigned int, unsigned int);
vector<Gene> inversion(vector<Gene>, unsigned int, unsigned int);

#endif
