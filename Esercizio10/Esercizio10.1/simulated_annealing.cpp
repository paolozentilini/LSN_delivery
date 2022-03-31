#include "simulated_annealing.h"


SimulatedAnnealing::SimulatedAnnealing(Random* rnd, unsigned int n_cities){
  _rnd = rnd;
  _n_cities = n_cities;
}

SimulatedAnnealing::~SimulatedAnnealing(){}

void SimulatedAnnealing::start_on_circumference(){
  vector<Gene> ch;
  Gene city;
  _rnd->init();
  for(unsigned int i=0; i<_n_cities; i++){
    double theta = _rnd->get_uniform(0,2*M_PI);
    city._x = cos(theta);
    city._y = sin(theta);
    city._identity = i+1;
    ch.push_back(city);
  }
  _chromosome1 = new Chromosome(ch);
  _chromosome2 = new Chromosome();
}

void SimulatedAnnealing::start_on_plane(){
  vector<Gene> ch;
  Gene city;
  _rnd->init();
  for(unsigned int i=0; i<_n_cities; i++){
    city._x = _rnd->get_uniform(-1,1);
    city._y = _rnd->get_uniform(-1,1);
    city._identity = i+1;
    ch.push_back(city);
  }
  _chromosome1 = new Chromosome(ch);
  _chromosome2 = new Chromosome();
}

void SimulatedAnnealing::mutation1(){
  vector<Gene> chromo = _chromosome1->get_sequence();
  unsigned int pos1 = int((_n_cities-1)*_rnd->get_uniform())+1;
  unsigned int pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
  while(pos1==pos2){
    pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
  }
  chromo = swap(chromo,pos1,pos2);
  _chromosome2->set_sequence(chromo);
}

void SimulatedAnnealing::mutation2(){
  vector<Gene> chromo = _chromosome1->get_sequence();
  unsigned int pos1 = int((_n_cities-1)*_rnd->get_uniform())+1;
  chromo = permutation(chromo,pos1);
  _chromosome2->set_sequence(chromo);
}

void SimulatedAnnealing::mutation3(){
  vector<Gene> chromo = _chromosome1->get_sequence();
  unsigned int pos1 = int((_n_cities-1)*_rnd->get_uniform())+1;
  unsigned int pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
  while(pos1==pos2){
    pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
  }
  chromo = permutation(chromo,pos1,pos2);
  _chromosome2->set_sequence(chromo);
}

void SimulatedAnnealing::mutation4(){
  vector<Gene> chromo = _chromosome1->get_sequence();
  unsigned int pos1 = int((_n_cities-1)*_rnd->get_uniform())+1;
  unsigned int pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
  while(pos1==pos2){
    pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
  }
  chromo = inversion(chromo,pos1,pos2);
  _chromosome2->set_sequence(chromo);
}


void SimulatedAnnealing::evolution_step(double beta){
    double r1=_rnd->get_uniform();
    if(r1<=0.25) mutation1();
    if(r1>=0.25 && r1<=0.5) mutation2();
    if(r1>=0.5 && r1<=0.75) mutation3();
    if(r1>=0.75) mutation4();
    double l1 = _chromosome1->fitness();
    double l2 = _chromosome2->fitness();
    double alpha;
    if(l2>=l1) alpha = exp(-beta*(l2-l1));        //Metropolis usando la fitness
    else alpha = 1;
    double r2 = _rnd->get_uniform();
    if(r2<=alpha){
      _chromosome1->set_sequence(_chromosome2->get_sequence());
    }
}

double SimulatedAnnealing::cost(){
  return _chromosome1->fitness();
}

Chromosome SimulatedAnnealing::get_chromosome(){
  return *_chromosome1;
}


vector<Gene> swap(vector<Gene> v, int pos1, int pos2){
  Gene temp1 = v.at(pos1);
  Gene temp2 = v.at(pos2);
  v.erase(v.begin() + pos1);
  v.insert(v.begin() + pos1, temp2);
  v.erase(v.begin() + pos2);
  v.insert(v.begin() + pos2, temp1);
  return v;
}

vector<Gene> permutation(vector<Gene> v,unsigned int pos){
  Gene temp1;
  vector<Gene> temp2,temp3;
  unsigned int size = v.size();
  temp1 = v.at(0);
  for(unsigned int i=1; i<pos; i++)   temp2.push_back(v.at(i));
  for(unsigned int k=pos; k<size; k++)   temp3.push_back(v.at(k));
  v.clear();
  v.push_back(temp1);
  for(unsigned int k=0; k<size-pos; k++)   v.push_back(temp3.at(k));
  for(unsigned int i=0; i<pos-1; i++)   v.push_back(temp2.at(i));
  return v;
}

vector<Gene> permutation(vector<Gene> v,unsigned int a, unsigned int b){
  vector<Gene> temp1,temp2,temp3;
  unsigned int size = v.size();
  unsigned int pos1, pos2;
  if(a<=b){
    pos1 = a;
    pos2 = b;
  }else{
    pos1 = b;
    pos2 = a;
  }
  for(unsigned int i=0; i<pos1; i++)      temp1.push_back(v.at(i));
  for(unsigned int j=pos1; j<pos2; j++)   temp2.push_back(v.at(j));
  for(unsigned int k=pos2; k<size; k++)   temp3.push_back(v.at(k));
  v.clear();
  for(unsigned int j=0; j<pos1; j++)      v.push_back(temp1.at(j));
  for(unsigned int k=0; k<size-pos2; k++) v.push_back(temp3.at(k));
  for(unsigned int i=0; i<pos2-pos1; i++) v.push_back(temp2.at(i));
  return v;
}

vector<Gene> inversion(vector<Gene> v, unsigned int a, unsigned int b){
  vector<Gene> temp1,temp2,temp3;
  unsigned int size = v.size();
  unsigned int pos1, pos2;
  if(a<=b){
    pos1 = a;
    pos2 = b;
  }else{
    pos1 = b;
    pos2 = a;
  }
  for(unsigned int i=0; i<pos1; i++)      temp1.push_back(v.at(i));
  for(unsigned int j=pos1; j<pos2; j++)   temp2.push_back(v.at(j));
  for(unsigned int k=pos2; k<size; k++)   temp3.push_back(v.at(k));
  v.clear();
  reverse(temp2.begin(),temp2.end());
  for(unsigned int j=0; j<pos1; j++)      v.push_back(temp1.at(j));
  for(unsigned int i=0; i<pos2-pos1; i++) v.push_back(temp2.at(i));
  for(unsigned int k=0; k<size-pos2; k++) v.push_back(temp3.at(k));
  return v;
}
