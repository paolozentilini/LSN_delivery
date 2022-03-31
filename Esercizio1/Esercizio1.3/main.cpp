#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "statistic.h"
#include "buffon_experiment.h"
#include "experiment.h"
#include <cmath>

using namespace std;

int main(){

	double length = 0.05;
	double distance = 0.1;

	Random* rnd = new Random();
	Experiment* buffon_exp = new Buffon_Experiment(rnd,length,distance, 1e3);
	Statistics stat(100, 1e4, rnd, buffon_exp);

	stat.data_blocking("../../Esercizio1/Risultati/Buffon_Experiment.dat");

	return 0;
}
