#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "random.h"
#include "data_structures.h"
#include "genetic_algorithm.h"

using namespace std;

int main(){

  unsigned int n_cities = 32;
  unsigned int n_chromosomes = 1000;
  unsigned int n_evolutions = 200;
  ofstream outf1,outf2,outf3;
  Random* rnd = new Random();
  GeneticAlgorithm* ga = new GeneticAlgorithm(rnd, n_cities, n_chromosomes, 0.05, 0.75);

  cout << "------MY GENETIC ALGORITHM------" << endl;
  cout << "-     Number of cities: " << n_cities << endl;
  cout << "-     Number of chromosomes: " << n_chromosomes << endl;
  cout << "-     Number of generations: " << n_evolutions << endl;
  cout << "--------------------------------"<< endl;
//Cities are placed on a circumference of unitary radius:
  cout << "1) Cities randomly placed on a circumference." << endl;
  ga->start_on_circumference();
  //Salvo su file il percorso iniziale
  vector<Chromosome> population = ga->get_population();
  Chromosome ch = population.at(n_chromosomes-1);
  Gene city;
  outf3.open("../../Esercizio9/Risultati/cities_positions_init.dat");
  for(unsigned int i=0; i<n_cities; i++){
    city = ch.get_gene(i);
    outf3 << city._x << "\t" << city._y << "\t" << city._identity << endl;
  }
  outf3.close();
  //Faccio girare l'algoritmo e contemporaneamente registro su file il fitness medio e il best fitness
  outf1.open("../../Esercizio9/Risultati/best_fit.dat");
  outf2.open("../../Esercizio9/Risultati/mean_fit.dat");
  for(unsigned int i=0; i<n_evolutions; i++){
    ga->evolution();
    if(i%10==0 && i!=0) cout << "Generation " << i << endl;
    outf1 << i << "\t" << ga->best_fit() << endl;
    outf2 << i << "\t" << ga->mean_fit() << endl;
  }
  outf1.close();
  outf2.close();
  //Salvo su file il percorso finale
  population = ga->get_population();
  ch = population.at(n_chromosomes-1);
  outf3.open("../../Esercizio9/Risultati/cities_positions.dat");
  for(unsigned int i=0; i<n_cities; i++){
    city = ch.get_gene(i);
    outf3 << city._x << "\t" << city._y << "\t" << city._identity << endl;
  }
  outf3.close();
//////////////////////////////////////////////////////////////////////////////////////////////
//Cities are positioned randomly into the plane [-1,1]x[-1,1]:
  cout << "2) Cities randomly placed on a plane." << endl;
  rnd->save_seed();
  ga->clear_population();
  ga->start_on_plane();
  //Salvo su file il percorso iniziale
  population = ga->get_population();
  ch = population.at(n_chromosomes-1);
  outf3.open("../../Esercizio9/Risultati/cities_positions1_init.dat");
  for(unsigned int i=0; i<n_cities; i++){
    city = ch.get_gene(i);
    outf3 << city._x << "\t" << city._y << "\t" << city._identity << endl;
  }
  outf3.close();
  //Faccio girare l'algoritmo e contemporaneamente registro su file il fitness medio e il best fitness
  outf1.open("../../Esercizio9/Risultati/best_fit1.dat");
  outf2.open("../../Esercizio9/Risultati/mean_fit1.dat");
  for(unsigned int i=0; i<n_evolutions; i++){
    ga->evolution();
    if(i%10==0 && i!=0) cout << "Generation " << i << endl;
    outf1 << i << "\t" << ga->best_fit() << endl;
    outf2 << i << "\t" << ga->mean_fit() << endl;
  }
  outf1.close();
  outf2.close();
  //Salvo su file il percorso finale:
  vector<Chromosome> population1 = ga->get_population();
  Chromosome chr = population1.at(n_chromosomes-1);
  outf3.open("../../Esercizio9/Risultati/cities_positions1.dat");
  for(unsigned int i=0; i<n_cities; i++){
    city = chr.get_gene(i);
    outf3 << city._x << "\t" << city._y << "\t" << city._identity << endl;
  }
  outf3.close();


return 0;
}
