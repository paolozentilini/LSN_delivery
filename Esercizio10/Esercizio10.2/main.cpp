#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "mpi.h"
#include "random.h"
#include "data_structures.h"
#include "genetic_algorithm.h"

using namespace std;

int main(int argc, char* argv[]){

//~~MPI initialisation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  int size, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status stat;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  unsigned int n_cities = 32;
  unsigned int n_chromosomes = 1000;
  unsigned int n_evolutions = 200;
  ofstream outf1,outf2,outf3;
  Random* rnd = new Random();
  GeneticAlgorithm* ga = new GeneticAlgorithm(rnd, n_cities, n_chromosomes, 0.05, 0.75);

//~~~~~Let's start with generating cities and sharing them across all processes~~~~~~~~~~

  vector<double> city_x(n_cities),city_y(n_cities);
  vector<int> city_identity(n_cities);
  if(rank==0){
    city_x.clear();
    city_y.clear();
    city_identity.clear();
    for(unsigned int i=0; i<n_cities; i++){
      double r1 = rnd->get_uniform(-1,1);
      double r2 = rnd->get_uniform(-1,1);
      city_x.push_back(r1);
      city_y.push_back(r2);
      city_identity.push_back(i+1);
    }
  }
  MPI_Bcast(&city_x.front(),n_cities,MPI_DOUBLE,0, MPI_COMM_WORLD);
  MPI_Bcast(&city_y.front(),n_cities,MPI_DOUBLE,0, MPI_COMM_WORLD);
  MPI_Bcast(&city_identity.front(),n_cities,MPI_INT,0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  vector<int> city_identity1(n_cities);
  vector<int> city_identity2(n_cities);
  vector<int> process_id = {0, 1, 2, 3};

  ga->start(city_x,city_y,city_identity);
  outf1.open("../../Esercizio10/Risultati/best_cost"+to_string(rank)+".dat");
  outf2.open("../../Esercizio10/Risultati/mean_cost"+to_string(rank)+".dat");
  for(unsigned int i=0; i<n_evolutions; i++){
    ga->evolution();
    if(i%20==0 && rank==0) cout << "Generazione " << i << endl;
    if(i%10==0){
      for(unsigned int j=0; j<2; j++){
        int tag1 = 2*j;
        int tag2 = 2*j+1;
        if(rank==process_id[tag1]){
          city_identity1 = ga->get_best_chromosome_sequence();
          MPI_Send(&city_identity1.front(), n_cities, MPI_INT, process_id[2*j+1], tag1, MPI_COMM_WORLD);
          MPI_Recv(&city_identity2.front(), n_cities, MPI_INT, process_id[2*j+1], tag2, MPI_COMM_WORLD, &stat);
          ga->set_best_chromosome(city_identity2);
        }
        if(rank==process_id[tag2]){
          city_identity2 = ga->get_best_chromosome_sequence();
          MPI_Recv(&city_identity1.front(), n_cities, MPI_INT, process_id[2*j], tag1, MPI_COMM_WORLD, &stat);
          MPI_Send(&city_identity2.front(), n_cities, MPI_INT, process_id[2*j], tag2, MPI_COMM_WORLD);
          ga->set_best_chromosome(city_identity1);
        }
      }
    }
    outf1 << i << "\t" << ga->best_fit() << endl;
    outf2 << i << "\t" << ga->mean_fit() << endl;
  }
  outf1.close();
  outf2.close();

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  vector<Chromosome> population = ga->get_population();
  Chromosome ch = population.at(n_chromosomes-1);
  Gene city;
  outf3.open("../../Esercizio10/Risultati/cities_positions_"+to_string(rank)+".dat");
  for(unsigned int i=0; i<n_cities; i++){
    city = ch.get_gene(i);
    outf3 << city._x << "\t" << city._y << "\t" << city._identity << endl;
  }
  outf3.close();

  MPI_Finalize();

return 0;
}
