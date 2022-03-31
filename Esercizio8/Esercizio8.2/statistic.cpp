#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "statistic.h"

using namespace std;

Statistics :: Statistics(unsigned int n_blocks, unsigned int n_measures, Integral* integral) {
	_n_blocks = n_blocks;
	_n_measures = n_measures;
	_n_per_blocks = n_measures/n_blocks;
	_integral = integral;
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

double Statistics :: mean_value(double mean_value, unsigned int counter){
	return mean_value/counter;
}

double Statistics :: dev_std(double sum, double sum_squared, unsigned int counter){
	return sqrt( (1./(counter-1))*( (1./counter)*sum_squared - ((1./counter)*sum)*((1./counter)*sum) ) );
}

void Statistics :: data_blocking(const char* results_file){
	ofstream outfile;
	outfile.open(results_file);
	_sum=0;
	_sum_squared=0;
	if (outfile.is_open()){
		for(unsigned int i=0; i<_n_blocks; i++){
			double single_measure=0;
			for(unsigned int j=0; j<_n_per_blocks; j++){
				single_measure += _integral -> single_evaluation();
			}
			single_measure /= _n_per_blocks;
			_sum += single_measure;
			_sum_squared += single_measure*single_measure;
			if(i!=0){	outfile << mean_value(_sum,i+1) << "\t\t" << dev_std(_sum,_sum_squared,i+1) << endl;}
		}
	} else cerr << "WARNING: not able to open " << results_file << " !" << endl;
	outfile.close();
}

vector<double> Statistics :: final_mean_sigma(){
	_sum=0;
	_sum_squared=0;
	vector<double> v;
	for(unsigned int i=0; i<_n_blocks; i++){
		double single_measure=0;
		for(unsigned int j=0; j<_n_per_blocks; j++){
			single_measure += _integral -> single_evaluation();
		}
		single_measure /= _n_per_blocks;
		_sum += single_measure;
		_sum_squared += single_measure*single_measure;
	}
	v.push_back(mean_value(_sum,_n_blocks));
	v.push_back(dev_std(_sum,_sum_squared,_n_blocks));
	return v;
}
