#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "random.h"
#include "data_structures.h"
#include "simulated_annealing.h"

using namespace std;

int main(){

  unsigned int n_cities = 32;
  //unsigned int n_evolutions = 100;

  ofstream outf1;
  Random* rnd = new Random();
  SimulatedAnnealing* ga = new SimulatedAnnealing(rnd, n_cities);

//Cities are placed on a circumference of unitary radius:

  ga->start_on_circumference();
  outf1.open("../../Esercizio10/Risultati/cost.dat");
  int j=1;
  for(double beta=0.; beta<150; beta+=1){
    for(int i=0; i<3000; i++){
      ga->evolution_step(beta);
    }
    if(int(beta)%10==0) cout << beta << "\t" << ga->cost() << endl;
    outf1 << beta << "\t" << ga->cost() << endl;
    j++;
  }
  outf1.close();

  Chromosome ch = ga->get_chromosome();
  Gene city;
  outf1.open("../../Esercizio10/Risultati/cities_positions.dat");
  for(unsigned int i=0; i<n_cities; i++){
    city = ch.get_gene(i);
    outf1 << city._x << "\t" << city._y << "\t" << city._identity << endl;
  }
  outf1.close();

//Cities are placed on a plane [-1,1]x[-1,1]:
  rnd->save_seed();
  ga->start_on_plane();
  outf1.open("../../Esercizio10/Risultati/cost1.dat");
  j=1;
  for(double beta=0.; beta<150; beta+=1){
    for(int i=0; i<3000; i++){
      ga->evolution_step(beta);
    }
    if(int(beta)%10==0) cout << beta << "\t" << ga->cost() << endl;
    outf1 << beta << "\t" << ga->cost() << endl;
    j++;
  }
  outf1.close();

  ch = ga->get_chromosome();
  outf1.open("../../Esercizio10/Risultati/cities_positions1.dat");
  for(unsigned int i=0; i<n_cities; i++){
    city = ch.get_gene(i);
    outf1 << city._x << "\t" << city._y << "\t" << city._identity << endl;
  }
  outf1.close();


return 0;
}
