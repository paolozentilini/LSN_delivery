#include <cmath>
#include "metropolis.h"

using namespace std;

Metropolis::Metropolis(Random* rnd, Distribution* p, double starting_point, double delta){
  _rnd=rnd;
  _p=p;
  _rnd->init();
  _delta=delta;
  _current_position = starting_point;
}

Metropolis::~Metropolis(){
  _rnd->save_seed();
}

double Metropolis::single_step(){
  double proposed_new_position;
  _rnd->set_parameters(_current_position,_delta);
  proposed_new_position = _rnd->get_random_number();
  double alpha = min( 1. , _p->eval(proposed_new_position)/_p->eval(_current_position) );
  double r=_rnd->get_uniform();
  if(r<=alpha) _current_position = proposed_new_position;
  return _current_position;
}

double Metropolis::acceptance(){
  double proposed_new_position;
  _rnd->set_parameters(_current_position,_delta);
  proposed_new_position = _rnd->get_random_number();
  double alpha = min( 1. , _p->eval(proposed_new_position)/_p->eval(_current_position) );
  double r=_rnd->get_uniform();
  if(r<=alpha) _current_position = proposed_new_position;
  return alpha;
}

void Metropolis::equilibrate(unsigned int n){
  double proposed_new_position;
  for(unsigned int k=0; k<n; k++){
    _rnd->set_parameters(_current_position,_delta);
    proposed_new_position = _rnd->get_random_number();
    double alpha = min( 1. , _p->eval(proposed_new_position)/_p->eval(_current_position) );
    double r=_rnd->get_uniform();
    if(r<=alpha) _current_position = proposed_new_position;
  }
}
