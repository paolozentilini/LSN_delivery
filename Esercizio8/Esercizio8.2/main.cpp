#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>

#include "statistic.h"
#include "random.h"
#include "uniform.h"
#include "distribution.h"
#include "distribution_trial.h"
#include "metropolis.h"
#include "integral.h"

using namespace std;

int main(){

	unsigned int n_blocks=100;																//Numero di blocchi
	unsigned int n_measures=5e6;															//Numero di misure totali
	double starting_point = 1.;																//Punto di partenza per il campionamento della catena di Markov
	ofstream outfile;

	Random *rnd;
	Distribution *distribution;
	Metropolis *mt;
	Integral *integral;

/*
	//Studio <H> al variare di sigma e mu
	outfile.open("../../Esercizio8/Risultati/energy_sigma_mu.dat");
	int i=0;
	vector<double> v;
	for(double mu=0.7; mu<=0.95; mu+=0.01 ){
		for(double sigma=0.5; sigma<=0.75; sigma+=0.01){
			rnd = new Uniform();
			distribution = new Distribution_trial(mu,sigma);
			mt = new Metropolis(rnd, distribution, starting_point, 2.67);
			integral= new Integral(mt,mu,sigma);
			Statistics stat(n_blocks,n_measures,integral);
			if(i%20==0) cout << "Ciclo " << i << endl;
			i++;
			v = stat.final_mean_sigma();
			outfile << mu << '\t' << sigma << '\t' << v.at(0) <<  '\t' << v.at(1) << endl;
		}
	}
	outfile.close();
*/

	//Valori ottimali di sigma e mu trovati:
	double sigma=0.61;
	double mu=0.8;

	//Calcolo statistico(tramite data_blocking) dell'energia che più si avvicina a E_GroundState
	rnd = new Uniform();
	distribution = new Distribution_trial(mu,sigma);
	mt = new Metropolis(rnd, distribution, starting_point, 2.67);
	integral= new Integral(mt,mu,sigma);
	Statistics stat(n_blocks,n_measures,integral);
	stat.data_blocking("../../Esercizio8/Risultati/energy.dat");

	//Generazione dei punti per l'istogramma
	unsigned int n_hist = 1e6;
	outfile.open("../../Esercizio8/Risultati/hist.dat");
	rnd = new Uniform();
	distribution = new Distribution_trial(mu,sigma);
	mt = new Metropolis(rnd, distribution, starting_point, 2.67);
	for(unsigned int i=0; i<n_hist; i++)
		outfile <<  mt->single_step() << endl;
	outfile.close();

 	//OTTIMIZZAZIONE DI A(x|y)
	/*
	rnd = new Uniform();
	distribution = new Distribution_trial(mu,sigma);
	double sum=0;
	int N=1e4;
	outfile.open("acceptance.dat");
	for(double delta=0.1; delta<5.; delta+=0.01){
			sum=0;
			mt = new Metropolis(rnd, distribution, starting_point, delta);
			for(int i=0; i<N; i++){
				sum += mt->acceptance();
			}
			outfile << delta << "		" << sum/N << endl;
	}
	outfile.close();
	*/

return 0;
}
