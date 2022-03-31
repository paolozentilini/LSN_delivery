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
    //Costruttore e distruttore:
    GeneticAlgorithm(Random* rnd, unsigned int n_cities, unsigned int n_chromosomes, double p_m, double p_c);
    ~GeneticAlgorithm();
    //Metodi per ottenere il numero di cromosomi e il vettore che rappresenta la popolazione
    unsigned int get_n_chromosomes(){return _n_chromosomes; }
    vector<Chromosome> get_population(){ return _population; }
    //Metodi per il calcolo della media della fitness e il best fitness
    double best_fit();
    double mean_fit();

    void start_on_circumference();  //Metodi che inizializzano le città sulla circonferenza o sul piano
    void start_on_plane();
    void select();                 //Operatore di selezione
    void crossover();              //Operatore di crossover
    void mutation1();              //Operatore di mutazione numero 1
    void mutation2();              //Operatore di mutazione numero 2
    void mutation3();              //Operatore di mutazione numero 3
    void mutation4();              //Operatore di mutazione numero 4
    void recreate_population();    //Funzione che sostituisce i cromosomi figli al posto dei padri
    void population_update();      //Caricamento della nuova popolazione
    void clear_population();       //Funzione che pulisce il vettore _population
    void evolution();              //Costruisce un'intera generazione

};
//Alcune funzioni utili:
bool my_cmp( Chromosome&, Chromosome&);  //Funzione booleana che confronta in base alla fitness
bool in( const Gene, vector<Gene>);   //Funzione che restituisce vero se il gene si trova nella sequenza di geni
//Funzioni che implementano inversioni e diverse permutazioni utili ai metodi Mutazione:
vector<Gene> swap(vector<Gene>, int, int);
vector<Gene> permutation(vector<Gene>, unsigned int);
vector<Gene> permutation(vector<Gene>, unsigned int, unsigned int);
vector<Gene> inversion(vector<Gene>, unsigned int, unsigned int);

#endif
