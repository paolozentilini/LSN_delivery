#include <cmath>
#include "metropolis.h"

using namespace std;

Metropolis::Metropolis(Random* rnd, Distribution* p, vector<double> starting_point, double delta){
  _rnd=rnd;
  _p=p;
  _rnd->init();
  _delta=delta;
  _current_position = starting_point;
}

Metropolis::~Metropolis(){
  _rnd->save_seed();
}

vector<double> Metropolis::single_step(){
  vector<double> proposed_new_position;
  for(unsigned int i=0; i<_current_position.size(); i++){
    _rnd->set_parameters(_current_position.at(i),_delta);
    proposed_new_position.push_back(_rnd->get_random_number());
  }
  double alpha = min( 1. , _p->eval(proposed_new_position)/_p->eval(_current_position) );
  double r=_rnd->get_uniform();
  if(r<=alpha) _current_position = proposed_new_position;
  return _current_position;
}

double Metropolis::acceptance(){
  vector<double> proposed_new_position;
  for(unsigned int i=0; i<_current_position.size(); i++){
    _rnd->set_parameters(_current_position.at(i),_delta);
    proposed_new_position.push_back(_rnd->get_random_number());
  }
  double appo = _p->eval(proposed_new_position)/_p->eval(_current_position);
  double alpha = min( 1. , appo );
  return alpha;
}

void Metropolis::equilibrate(unsigned int n){
  vector<double> proposed_new_position;
  for(unsigned int k=0; k<n; k++){
    for(unsigned int i=0; i<_current_position.size(); i++){
      _rnd->set_parameters(_current_position.at(i),_delta);
      proposed_new_position.push_back(_rnd->get_random_number());
    }
    double appo = _p->eval(proposed_new_position)/_p->eval(_current_position);
    double alpha = min( 1. , appo );
    double r=_rnd->get_uniform();
    if(r<=alpha){
      for(unsigned int j=0; j<_current_position.size(); j++) _current_position[j]=proposed_new_position.at(j);
    }
  }
}
