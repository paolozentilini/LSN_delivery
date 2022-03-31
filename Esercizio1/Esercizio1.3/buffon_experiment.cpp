#include "buffon_experiment.h"
#include "experiment.h"
#include "random.h"
#include <cmath>

using namespace std;

Buffon_Experiment :: Buffon_Experiment(Random* rnd,double length, double distance, unsigned int n_throws){
	_rnd=rnd;
	_rnd->Init();
	_pi=_rnd->piGreco(100000);
	_length=length;
	_distance = distance;
	_n_throws= n_throws;
}

Buffon_Experiment :: ~Buffon_Experiment(){}


vector<double> Buffon_Experiment :: needle(){	//Vector representation of a random needle
	vector<double> needle;
	vector<double> position = _rnd->position2D();
	double theta = _rnd->get_uniform(0,2*_pi);
	//needle[0]=position1[0];
	needle.push_back(position.at(1));
	//needle[2]=position1[0]+length*cos(theta);
	needle.push_back(position.at(1)+_length*sin(theta));
	return needle;
}

double Buffon_Experiment :: single_measurement(){
	unsigned int n_hit=0;
	for(unsigned int i=0; i<_n_throws; i++){
		vector<double> Needle = needle();
		for(double y=0; y<=1; y+=_distance){
			if(Needle.at(0)<y && Needle.at(1)>y) n_hit++;
			if(Needle.at(0)>y && Needle.at(1)<y) n_hit++;
		}
	}
	return (2*_length*_n_throws)/(_distance*n_hit);
}
