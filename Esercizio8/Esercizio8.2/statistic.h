#ifndef __Statistics__h_
#define __Statistics__h_

#include "integral.h"

class Statistics {

  protected:
    unsigned int _n_blocks;
    unsigned int _n_measures;
    unsigned int _n_per_blocks;
    Integral* _integral;
    double _sum;
    double _sum_squared;

  public:
    //Costruttori e distruttore
    Statistics(unsigned int n_blocks, unsigned int n_measures, Integral*);
    ~Statistics();
    //Metodi per il calcolo della media e della varianza
    double mean_value(std::vector<double>, unsigned int);
    double dev_std(std::vector<double>, unsigned int);
    double mean_value(double sum, unsigned int);
    double dev_std(double sum, double sum_squared, unsigned int);
    //Metodo che stampa su file i risultati dell'analisi statistica del random walk
    void data_blocking(const char* filename);
    //Metodo che restituisce la media finale del data_blocking e la deviazione standard:
    std::vector<double> final_mean_sigma();

};

#endif
