#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "statistic.h"
#include "random.h"
#include "continous_random_walk.h"
#include "discrete_random_walk.h"
#include "random_walk.h"

using namespace std;

int main(){

	//Dichiarazione dei costruttori Random e DiscreteRW + calcolo del random walk discreto e stampa su file
	Random *rnd = new Random();
	RandomWalk *rw = new DiscreteRW(1,rnd);
	Statistics stat(100,5e5,rw);
	stat.random_walk_data("../../Esercizio2/Risultati/discrete_random_walk.dat");
	//Dichiarazione dei costruttori Random e ContinousRW + calcolo del random walk continuo e stampa su file
	rnd = new Random();
	rw = new ContinousRW(1,rnd);
	Statistics stat1(100,5e5,rw);
	stat.random_walk_data("../../Esercizio2/Risultati/continous_random_walk.dat");


return 0;
}
