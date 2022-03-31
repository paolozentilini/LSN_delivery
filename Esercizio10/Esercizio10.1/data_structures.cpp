#include "data_structures.h"

Chromosome::Chromosome(std::vector<Gene> sequence){
  _sequence = sequence;
  _number_of_genes = _sequence.size();
  random_shuffle( _sequence.begin()+1, _sequence.end());    //Comando che "misclela" casualmente gli elementi di sequence
}

Chromosome::Chromosome(){}

Chromosome::~Chromosome(){}

Gene Chromosome::get_gene(unsigned int i){
  return _sequence.at(i);
}

void Chromosome::set_gene(unsigned int i, Gene gene){
  _sequence.at(i) = gene;
}

void Chromosome::set_sequence(std::vector<Gene> sequence){
  _sequence = sequence;
}

double Chromosome::fitness(){
  std::vector<Gene>::iterator it, init, final;
  double distance=0;
  init = _sequence.begin();
  final =  _sequence.end();
  for(it = init; it < final-1; it++)
    distance += sqrt( pow((it->_x)-((it+1)->_x),2) + pow((it->_y)-((it+1)->_y),2) );
  distance += sqrt( pow(((final-1)->_x)-(init->_x),2) + pow(((final-1)->_y)-(init->_y),2) );
  return distance;
}

void Chromosome::set_fitness(double fitness){
  _fitness = fitness;
}

void Chromosome::view(){
  Gene gene;
  for(unsigned int i=0; i<_number_of_genes; i++){
    gene = _sequence.at(i);
    std::cout << gene._identity << '\t' << gene._x << "\t" << gene._y << "\n";
  }
  std::cout << "Fitness:" << fitness() << "\n";
}
