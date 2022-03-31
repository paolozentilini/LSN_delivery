#include "genetic_algorithm.h"


GeneticAlgorithm::GeneticAlgorithm(Random* rnd, unsigned int n_cities, unsigned int n_chromosomes,double p_m, double p_c){
  _rnd = rnd;
  _rnd->init();
  _n_cities = n_cities;
  _n_chromosomes = n_chromosomes;
  _mutation_probability = p_m;
  _crossover_probability = p_c;
}

GeneticAlgorithm::~GeneticAlgorithm(){}

void GeneticAlgorithm::start(vector<double> x, vector<double> y, vector<int> id){
  vector<Gene> ch;
  Gene city;
  for(unsigned int i=0; i<_n_cities; i++){
    city._x = x.at(i);
    city._y = y.at(i);
    city._identity = id.at(i);
    ch.push_back(city);
  }
  Chromosome* chromo;
  for(unsigned int j=0; j<_n_chromosomes; j++){
    chromo = new Chromosome(ch);
    _population.push_back(*chromo);
  }
}

void GeneticAlgorithm::select(){
  for(unsigned int i=0; i<_n_chromosomes; i++)
    _population.at(i).set_fitness(_population.at(i).fitness());       //Calcolo del fitness per ciascun cromosoma
  sort(_population.begin(), _population.end(), my_cmp);               //Ordinamento della popolazione secondo il fitness
  double r1 = _rnd->get_uniform();                                    //L'ultimo elemento ha il fitness migliore (più basso)
  double r2 = _rnd->get_uniform();
  double p = 0.2;                                                     //valore scelto per l'esponente conveniente
  _position_chromosome1 = int(_n_chromosomes*pow(r1,p));              //Scelta dei due cromosomi genitori
  _position_chromosome2 = int(_n_chromosomes*pow(r2,p));
  _chromosome1 = _population.at(_position_chromosome1);
  _chromosome2 = _population.at(_position_chromosome2);
}

void GeneticAlgorithm::crossover(){
  double r1 = _rnd->get_uniform();
  if(r1 <= _crossover_probability){
      unsigned int cut = int((_n_cities-2)*_rnd->get_uniform())+1;
      for(unsigned int k=0; k<cut; k++){
      _offspring1.push_back(_chromosome1.get_gene(k));
      _offspring2.push_back(_chromosome2.get_gene(k));
    }
    for(unsigned int k=0; k<_chromosome1.get_number_of_genes(); k++){
      if(in(_chromosome2.get_gene(k), _offspring1)==false) _offspring1.push_back(_chromosome2.get_gene(k));
      if(in(_chromosome1.get_gene(k), _offspring2)==false) _offspring2.push_back(_chromosome1.get_gene(k));
    }
  }else{
    _offspring1 = _chromosome1.get_sequence();
    _offspring2 = _chromosome2.get_sequence();
  }
}

void GeneticAlgorithm::mutation1(){
  double r = _rnd->get_uniform();
  if(r <= _mutation_probability){
    unsigned int pos1 = int((_n_cities-1)*_rnd->get_uniform())+1;
    unsigned int pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
    while(pos1==pos2){
      pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
    }
    _offspring1 = swap(_offspring1,pos1,pos2);
    _offspring2 = swap(_offspring2,pos1,pos2);
  }
}

void GeneticAlgorithm::mutation2(){
  double r = _rnd->get_uniform();
  if(r <= _mutation_probability){
    unsigned int pos1 = int((_n_cities-1)*_rnd->get_uniform())+1;
    unsigned int pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
    _offspring1 = permutation(_offspring1,pos1);
    _offspring2 = permutation(_offspring2,pos2);
  }
}

void GeneticAlgorithm::mutation3(){
  double r = _rnd->get_uniform();
  if(r <= _mutation_probability){
    unsigned int pos1 = int((_n_cities-1)*_rnd->get_uniform())+1;
    unsigned int pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
    while(pos1==pos2){
      pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
    }
    _offspring1 = permutation(_offspring1,pos1,pos2);
    _offspring2 = permutation(_offspring2,pos1,pos2);
  }
}

void GeneticAlgorithm::mutation4(){
  double r = _rnd->get_uniform();
  if(r <= _mutation_probability){
    unsigned int pos1 = int((_n_cities-1)*_rnd->get_uniform())+1;
    unsigned int pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
    while(pos1==pos2){
      pos2 = int((_n_cities-1)*_rnd->get_uniform())+1;
    }
    _offspring1 = inversion(_offspring1,pos1,pos2);
    _offspring2 = inversion(_offspring2,pos1,pos2);
  }
}

void GeneticAlgorithm::recreate_population(){
  _chromosome1.set_sequence(_offspring1);
  _chromosome2.set_sequence(_offspring2);
  _new_population.push_back(_chromosome1);
  _new_population.push_back(_chromosome2);
  _offspring1.clear();
  _offspring2.clear();
}

void GeneticAlgorithm::population_update(){
  _population = _new_population;
  _new_population.clear();
}

void GeneticAlgorithm::clear_population(){
  _population.clear();
}

void GeneticAlgorithm::evolution(){
  for(unsigned int i=0; i<_n_chromosomes/2; i++){
    select();
    crossover();
    mutation1();
    mutation2();
    mutation3();
    mutation4();
    recreate_population();
  }
  population_update();
}

vector<int> GeneticAlgorithm::get_best_chromosome_sequence(){
  sort(_population.begin(), _population.end(), my_cmp);
  Chromosome chromo =_population.at(_population.size()-1);
  vector<int> id_vector;
  Gene gene;
  for(unsigned int i=0; i<_n_cities; i++){
    gene = chromo.get_gene(i);
    id_vector.push_back(gene._identity);
  }
  return id_vector;
}

void GeneticAlgorithm::set_best_chromosome(vector<int> label_sequence){
  Chromosome chromo =_population.at(_population.size()-1);
  Gene city, new_city;
  vector<Gene> sequence;
  for(unsigned int i=0; i<_n_cities; i++){
    for(unsigned int j=0; j<_n_cities; j++){
      city = chromo.get_gene(j);
      if(label_sequence.at(i)==city._identity){
        new_city._x = city._x;
        new_city._y = city._y;
        new_city._identity = city._identity;
        sequence.push_back(new_city);
      }
    }
  }
  chromo.set_sequence(sequence);
  _population.pop_back();
  _population.push_back(chromo);
}

double GeneticAlgorithm::best_fit(){
  sort(_population.begin(), _population.end(), my_cmp);
  return _population.at(_population.size()-1).get_fitness();
}

double GeneticAlgorithm::mean_fit(){
  double sum =0;
  unsigned int n = get_n_chromosomes();
  for(unsigned int i = 0; i < n; i++) {
    sum += _population.at(i).get_fitness();
  }
  return sum/n;
}


bool my_cmp(Chromosome& a, Chromosome& b){
  // smallest comes last
  return a.get_fitness() > b.get_fitness(); //>
}

bool in( const Gene element, vector<Gene> chromosome){
  for(unsigned int i=0; i<chromosome.size(); i++)
    if(element._identity == chromosome.at(i)._identity) return true;
  return false;
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
