#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "statistic.h"
#include "random.h"
#include "experiment.h"

using namespace std;

Statistics :: Statistics(unsigned int n_blocks, unsigned int n_measures, Random* random_generator, Experiment* experiment) {
	_n_blocks = n_blocks;
	_n_measures = n_measures;
	_n_per_blocks = n_measures/n_blocks;
	_random = random_generator;
	_experiment = experiment;
	_random -> Init();
}

Statistics :: Statistics(Random* random_generator) {
	_random = random_generator;
	_random -> Init();
}

Statistics :: ~Statistics() {}

double Statistics :: mean_value(vector<double> v, unsigned int counter){
	double mean_value = 0.;
	for (unsigned int i = 0; i < counter; i++)
		mean_value += v.at(i);
	mean_value /= counter;
	return mean_value;
}

double Statistics :: dev_std(vector<double> v, unsigned int counter){
	double sum_squared=0;
	double sum=0;
	for (unsigned int i=0; i<counter; i++) {
		sum += v.at(i);
		sum_squared += v.at(i)*v.at(i);
	}
	double dev_std = sqrt( (1./(counter-1))*( (1./counter)*sum_squared - ((1./counter)*sum)*((1./counter)*sum) ) );
	return dev_std;
}

double Statistics :: chi_squared(vector<int> counter, double expected_value){
	double chi_squared=0;
	for (unsigned int i = 0; i < counter.size(); i++)
		chi_squared += (( counter.at(i) - expected_value )*( counter.at(i) - expected_value))/expected_value;			//chi squared calculation
	return chi_squared;
}

void Statistics :: chi_squared_test(unsigned int n, unsigned int throws_per_block, const char* results_file){
	vector<double> chi_squared_measurements;
	vector<int> counter;
	double expected_value = throws_per_block/n;
	double appo;
	ofstream outfile;
	outfile.open(results_file);
	if (outfile.is_open()){
		for (unsigned int j = 0; j < n; j++) {
			counter.assign(n,0);
			for(unsigned int k=0; k<throws_per_block; k++){
				appo = _random -> get_uniform();
				for(unsigned int i=0; i<n; i++){
					if ( (double(i)/double(n) < appo) && (appo < double(i+1)/double(n))){ counter[i] += 1; }
				}
			}
			chi_squared_measurements.push_back(chi_squared(counter,expected_value));
		}
		_random -> SaveSeed();
		for (unsigned int j = 0 ; j < n; j++) {
			outfile << chi_squared_measurements.at(j) << endl;
		}
	} else cerr << "WARNING: not able to open " << results_file << " !" << endl;
	outfile.close();
}

void Statistics :: data_blocking(const char* results_file){
	vector<double> measurements_vector;
	ofstream outfile;
	outfile.open(results_file);
	if (outfile.is_open()){
		for(unsigned int i=0; i<_n_blocks; i++){
			double single_measure=0;
			for(unsigned int j=0; j<_n_per_blocks; j++){
				single_measure += _experiment -> single_measurement();
			}
			single_measure /= _n_per_blocks;
			measurements_vector.push_back(single_measure);
		}
		_random -> SaveSeed();
		for (unsigned int i = 2; i <= measurements_vector.size(); i++) {
			outfile << mean_value(measurements_vector,i) << "\t\t" << dev_std(measurements_vector,i) << endl;
		}
	} else cerr << "WARNING: not able to open " << results_file << " !" << endl;
	outfile.close();
}
