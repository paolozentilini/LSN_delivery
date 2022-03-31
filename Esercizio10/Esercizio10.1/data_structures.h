#ifndef __DataStructures__h_
#define __DataStructures__h_

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "random.h"

struct Gene{
  double _x;
  double _y;
  int _identity;
};

class Chromosome{
  protected:
    unsigned int _number_of_genes;
    std::vector<Gene> _sequence;
    double _fitness;
  public:
    Chromosome(std::vector<Gene>);
    Chromosome();
    ~Chromosome();
    //Methods:
    Gene get_gene(unsigned int i);           //Get the i th gene of the sequence
    void set_gene(unsigned int i, Gene);      //Set the i th gene of the sequence
    unsigned int get_number_of_genes(){ return _number_of_genes;}
    void set_sequence(std::vector<Gene>);
    std::vector<Gene> get_sequence(){ return _sequence;}

    double fitness();                        //Fitness calcolous
    double get_fitness(){ return _fitness;}
    void set_fitness(double fitness);

    void view();
};


#endif
