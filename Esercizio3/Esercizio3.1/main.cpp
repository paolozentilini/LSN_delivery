#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "statistic.h"
#include "random.h"
#include "continous_wiener_process.h"
#include "discrete_wiener_process.h"
#include "wiener_process.h"
#include "plain_vanilla.h"

using namespace std;

int main(){

//Parametri che verranno passati agli oggetti per il calcolo statistico dei prezzi delle opzioni PlainVanilla

	double S_0=100;
 	unsigned int n_steps=100;
	double T=1;
	double mu=0.1;
	double sigma=0.25;
	double strike=110;
	unsigned int n_blocks=100;
	unsigned int n_measures=1e5;
	Random *rnd = new Random();

///////////////////////Processo di Wiener discreto//////////////////////////////

	WienerProcess *wp = new DiscreteWienerProcess(S_0,T,n_steps,mu,sigma,rnd);

	PlainVanilla* option = new PlainVanilla(strike,mu,wp,"Call");
	Statistics stat(n_blocks,n_measures,option);
	stat.data_blocking("../../Esercizio3/Risultati/discrete_call.dat");

	option = new PlainVanilla(strike,mu,wp,"Put");
	Statistics stat1(n_blocks,n_measures,option);
	stat1.data_blocking("../../Esercizio3/Risultati/discrete_put.dat");

///////////////////////Processo di Wiener continuo//////////////////////////////

	wp = new ContinousWienerProcess(S_0, T, mu, sigma, rnd);

	option = new PlainVanilla(strike,mu,wp,"Call");
	Statistics stat2(n_blocks,n_measures,option);
	stat2.data_blocking("../../Esercizio3/Risultati/continous_call.dat");

	option = new PlainVanilla(strike,mu,wp,"Put");
	Statistics stat3(n_blocks,n_measures,option);
	stat3.data_blocking("../../Esercizio3/Risultati/continous_put.dat");

return 0;
}
