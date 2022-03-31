#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "statistic.h"
#include "random.h"
#include "uniform_sampling.h"
#include "coseno.h"
#include "montecarlo.h"
#include "importance_sampling.h"
#include "g_function.h"
#include "f_function.h"

using namespace std;

int main(){
	//Dichiaro prima rnd, funzionebase e montecarlo per poi usare il metodo data_blocking della classe Statistics.
	//Uniform_Sampling:
	Random *rnd = new Random();
	FunzioneBase *f = new Coseno();
	MonteCarlo *mc = new Uniform_Sampling(0,1,f,rnd,100);
	Statistics stat(100,1e4,mc);
	stat.data_blocking("../../Esercizio2/Risultati/uniform_sampling.dat");
	//Importance_Sampling esatto
	rnd = new Random();
	f = new G_Function();
	mc = new Importance_Sampling(0,1,f,rnd,100,0);
	Statistics stat1(100,1e4,mc);
	stat1.data_blocking("../../Esercizio2/Risultati/importance_sampling1.dat");
	//Importance_Sampling con distribuzione approssimata:
	rnd = new Random();
	f = new F_Function();
	mc = new Importance_Sampling(0,1,f,rnd,100,1);
	Statistics stat2(100,1e4,mc);
	stat2.data_blocking("../../Esercizio2/Risultati/importance_sampling2.dat");

return 0;
}
