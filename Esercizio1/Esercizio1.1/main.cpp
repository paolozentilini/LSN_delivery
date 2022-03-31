#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "statistic.h"
#include "random.h"
#include "experiment.h"
#include "integral.h"
#include "integral_1.h"

using namespace std;

int main(){
	//Esercizio1.1:
	Random *rand = new Random();									//Creo dinamicamente il costruttore di Random
	Experiment *integral = new Integral();				//Creo dinamicamente il costruttore di Integral partendo dalla classe virtuale madre Experiment
	Statistics stat(100, 1e6, rand, integral);		//Creo staticamente il costruttore della classe Statistics con i parametri necessari
	stat.data_blocking("../../Esercizio1/Risultati/measurements.dat");	//Evoco il metodo data_blocking e salvo i risultati
	//Esercizio1.2:
	Experiment *integral_1 = new Integral_1();		//Creo dinamicamente il costruttore di Integral_1 partendo dalla classe virtuale madre Experiment
	Statistics stat_1(100, 1e6, rand, integral_1);//Creo staticamente il costruttore della classe Statistics con i nuovi parametri
	stat_1.data_blocking("../../Esercizio1/Risultati/measurements1.dat"); //Evoco il metodo data_blocking e salvo i risultati
	//Esercizio1.3:
	Statistics stat_2(rand); 											//Creo staticamente il secondo costruttore della classe Statistics
	stat_2.chi_squared_test(100, 1e4, "../../Esercizio1/Risultati/chi_squared_test10000.dat"); //Svolgo il test del chi quadro e salvo i risultati

return 0;
}
