#ifndef __DataStructures__h_
#define __DataStructures__h_

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "random.h"

//La città è caratterizzata da una posizione e da un'etichetta (numero progressivo di generazione)
struct Gene{
  double _x;
  double _y;
  int _identity;
};
//Classe che modifica/costruisce il vettore di città
class Chromosome{
  protected:
    unsigned int _number_of_genes;
    std::vector<Gene> _sequence;
    double _fitness;
  public:
    //Costruttori e distruttore:
    Chromosome(std::vector<Gene>);
    Chromosome();
    ~Chromosome();
    //Metodi:
    Gene get_gene(unsigned int i);                                //Restituisce il gene i esimo
    void set_gene(unsigned int i, Gene);                          //Assegna al gene i la struct gene
    unsigned int get_number_of_genes(){ return _number_of_genes;} //Restituisce il numero di geni
    void set_sequence(std::vector<Gene>);                         //Assegna alla sequenza data membro la sequenza di geni pr argomento
    std::vector<Gene> get_sequence(){ return _sequence;}          //Restituisce la sequenza di geni

    double fitness();                                             //Calcolo della fitness
    double get_fitness(){ return _fitness;}                       //Funzione che restituisce la fitness
    void set_fitness(double fitness);                             //Funzione che assegna la fitness

    void view();                                                  //Metodo di controllo(visualizzo il cromosoma a schermo)
};


#endif
