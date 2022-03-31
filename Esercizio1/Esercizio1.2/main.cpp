#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "random.h"

using namespace std;

double mean_value(vector<double>); //Restituisce la media dei dati contenuti in un vector

int main(){

	Random *rnd = new Random();										//Costruttore dinamico di Random
	ofstream myfile,myfile1,myfile2,myfile3;			//Dichiarazione deglli operatori su file
	unsigned int M=1e4;														//Numero di dati per file
	vector<double> data,data1,data2;							//Contenitori dati
	//Inizializzo il generatore e apro i file sui quali salvare i dati.
	rnd->Init();
	myfile.open("../../Esercizio1/Risultati/tlc.dat");
	myfile1.open("../../Esercizio1/Risultati/tlc1.dat");
	myfile2.open("../../Esercizio1/Risultati/tlc2.dat");
	myfile3.open("../../Esercizio1/Risultati/tlc3.dat");
	//Ciclo da 0 a 10^4. Per ogni ciclo estraggo in successione N= 1,2,10,100 numeri random da tre diverse distribuzioni:
	//uniforme, esponenziale e Cauchy-Lorentz.
	//Per ogni N fissato faccio la media degli N numeri random e questo è il dato da salvare su file.
	for(unsigned int j=0; j<M; j++){
		myfile << rnd->get_uniform() << "\t" << rnd->get_esponential(1) << "\t" <<  rnd->get_cauchy_lorentz(0,1) << endl;
		for(unsigned int k=0; k<2; k++){
			data.push_back(rnd->get_uniform());
			data1.push_back(rnd->get_esponential(1));
			data2.push_back(rnd->get_cauchy_lorentz(0,1));
		}
		myfile1 << mean_value(data) << "\t" << mean_value(data1) << "\t" << mean_value(data2) << endl;
		data.clear();
		data1.clear();
		data2.clear();
		for(unsigned int k=0; k<10; k++){
			data.push_back(rnd->get_uniform());
			data1.push_back(rnd->get_esponential(1));
			data2.push_back(rnd->get_cauchy_lorentz(0,1));
		}
		myfile2 << mean_value(data) << "\t" << mean_value(data1) << "\t" << mean_value(data2) << endl;
		data.clear();
		data1.clear();
		data2.clear();
		for(unsigned int k=0; k<100; k++){
			data.push_back(rnd->get_uniform());
			data1.push_back(rnd->get_esponential(1));
			data2.push_back(rnd->get_cauchy_lorentz(0,1));
		}
		myfile3 << mean_value(data) << "\t" << mean_value(data1) << "\t" << mean_value(data2) << endl;
		data.clear();
		data1.clear();
		data2.clear();
	}
	myfile.close();
	myfile1.close();
	myfile2.close();
	myfile3.close();
	rnd->SaveSeed();

return 0;
}

double mean_value(vector<double> v){
	double appo=0;
	unsigned int N= v.size();
	for(unsigned int i=0; i<N; i++)
		appo+=v[i];
	return appo/N;
}
