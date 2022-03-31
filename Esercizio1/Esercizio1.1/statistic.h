#ifndef __Statistics__
#define __Statistics__

#include "random.h"
#include "experiment.h"

class Statistics {

  protected:

    unsigned int _n_blocks;                   //Numero di blocchi
    unsigned int _n_measures;                 //Numero di misure
    unsigned int _n_per_blocks;               //Numero di misure per blocco
    Random *_random;                          //Puntatore alla classe Random
    Experiment *_experiment;                  //Puntatore alla classe Experiment

  public:
    //Costruttori e distruttore della classe Statistics (Entrambi i costruttori inizializzano il costruttore di Random):
    Statistics(unsigned int n_blocks, unsigned int n_measures, Random*, Experiment*);
    Statistics(Random*);
    ~Statistics();
    //Metodi:
    double mean_value(std::vector<double>, unsigned int); //Da un vector di dati calcola il valor medio dei dati da 0 a unsigned int
    double dev_std(std::vector<double>, unsigned int);    //Da un vector di dati calcola la deviazione standard dei dati da 0 a unsigned int
    double chi_squared(std::vector<int>, double); //Da un vector di dati calcola il chi quadro dei dati secondo uno specifico valore di aspettazione
    void chi_squared_test(unsigned int n, unsigned int throws_per_block, const char*); //Metodo che riproduce il test del chi quadro
    void data_blocking(const char*);  //Metodo che implementa il blocking e restituisce direttamente i risultati su file

};

#endif
