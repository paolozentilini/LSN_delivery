#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "statistic.h"
#include "random.h"
#include "gaussian.h"
#include "uniform.h"
#include "distribution.h"
#include "distribution_2p.h"
#include "distribution_1s.h"
#include "metropolis.h"
#include "integral.h"

using namespace std;

int main(){

	unsigned int n_blocks=100;																//Numero di blocchi
	unsigned int n_measures=1e6;															//Numero di misure totali
	double a0=1;																							//Raggio di Bohr
	vector <double> starting_point = {a0, a0, a0};						//Punto di partenza per il campionamento della catena di Markov
	vector <double> starting_point1 = {2*a0, 2*a0, 3*a0};

	Random *rnd;
	Distribution *distribution;
	Metropolis *mt;
	Integral *integral;

	//Studio <r> per |Ψ1,0,0(x,y,z)|^2

	rnd = new Uniform();
	distribution = new Distribution_1s(a0);
	mt = new Metropolis(rnd, distribution, starting_point, 1.22);
	integral= new Integral(mt);
	Statistics stat(n_blocks,n_measures,integral);
	stat.data_blocking("../../Esercizio5/Risultati/1s_uniform.dat");

	rnd = new Gaussian();
	distribution = new Distribution_1s(a0);
	mt = new Metropolis(rnd, distribution, starting_point, 0.76);
	integral= new Integral(mt);
	Statistics stat1(n_blocks,n_measures,integral);
	stat1.data_blocking("../../Esercizio5/Risultati/1s_gaussian.dat");

	//Studio <r> per |Ψ2,1,0(x,y,z)|^2

	rnd = new Uniform();
	distribution = new Distribution_2p(a0);
	mt = new Metropolis(rnd, distribution, starting_point1, 2.95);
	integral= new Integral(mt);
	Statistics stat2(n_blocks,n_measures,integral);
	stat2.data_blocking("../../Esercizio5/Risultati/2p_uniform.dat");

  rnd = new Gaussian();
	distribution = new Distribution_2p(a0);
	mt = new Metropolis(rnd, distribution, starting_point1, 1.89);
	integral= new Integral(mt);
	Statistics stat3(n_blocks,n_measures,integral);
	stat3.data_blocking("../../Esercizio5/Risultati/2p_gaussian.dat");
/*
 //OTTIMIZZAZIONE DI A(x|y)
	//rnd = new Uniform();
	rnd = new Gaussian();
	//distribution = new Distribution_1s(a0);
	distribution = new Distribution_2p(a0);
	double sum=0;
	vector <double> starting_point1 = {2*a0, 2*a0, 3*a0};
	ofstream out;
	int N=1e4;
	out.open("acceptance.dat");
	for(double delta=0.1; delta<4.; delta+=0.01){
			sum=0;
			mt = new Metropolis(rnd, distribution, starting_point1, delta);
			for(int i=0; i<N; i++){
				sum += mt->acceptance();
			}
			out << delta << "		" << sum/N << endl;
	}
	out.close();
*/
return 0;
}
