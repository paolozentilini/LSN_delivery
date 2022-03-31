/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  for(temp=2.0; temp>0.499; temp-=0.01){
    Input(); //Inizialization
    for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation

      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep){

        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk);   //Print results for current block
    }
    ConfFinal(); //Write final configuration
    cout << "----------------------------" << endl << endl;
  }
  return 0;
}

void Input(void){

  ifstream ReadInput;
//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

//Read input informations
  ReadInput.open("input.dat");

  //ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> init;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

  if(init ==1){
    ofstream out;
    out.open("config.0");
    //initial configuration
    for (int i=0; i<nspin; ++i){
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
      out << s[i] << endl;
    }
    out.close();
    //Evaluate energy etc. of the initial configuration
    Measure();
    //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
  }else{
    ifstream infile;
    infile.open("config.final");
    for (int i=0; i<nspin; ++i)
      infile >> s[i];
    infile.close();
    //Evaluate energy etc. of the configuration
    Measure();
    //Print values for the potential energy and virial
    cout << "Energy = " << walker[iu]/(double)nspin << endl;
  }
}

void Reset(int iblk){ //Reset block averages
   if(iblk == 1){
       for(int i=0; i<n_props; ++i){
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }
   for(int i=0; i<n_props; ++i){
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Move(int metro){

  int o;
  double p,s_proposed, alpha, r, delta_E;

  for(int i=0; i<nspin; ++i){
    o = (int)(rnd.Rannyu()*nspin);    //Select randomly a spin (for C++ syntax, 0 <= o <= nspin-1)
    attempted++;
    if(metro==1){ //Metropolis
      if(rnd.Rannyu() >= 0.5) s_proposed = 1;
      else s_proposed = -1;
      delta_E = Boltzmann(s_proposed,o)-Boltzmann(s[o],o);
      if(delta_E < 0){
         s[o] = s_proposed;
         accepted++;
      }else{
        alpha = exp(-delta_E/temp);
        r = rnd.Rannyu();
        if(r<alpha){
           s[o] = s_proposed;
           accepted++;
        }
      }
    }else{ //Gibbs sampling
      accepted++;                     //Tutte le mosse sono accettate
      double r = rnd.Rannyu();
      if(s[o]==1){
        p= 1./(1. + exp(2*Boltzmann(s[o], o)/temp));  //Probabilità congiunta per s = +1
        if(r >= p) s[o] =-1;
      }else{
        p= 1./(1. + exp(-2*Boltzmann(s[o], o)/temp)); //Probabilità congiunta per s = -1
        if(r <= p) s[o] = 1;
      }
    }
  }
}

double Boltzmann(double sm, int ip){     //Restituisce l'energia legata allo spin sm di indice ip
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure(){

  double m = 0.0, u = 0.0;
  double spin;
  //cycle over spins
  for (int i=0; i<nspin; ++i){
    spin = s[i];
    u += -J * spin * s[Pbc(i+1)] - 0.5 * h * (spin + s[Pbc(i+1)]);
    m += spin;
  }
  walker[iu] = u;
  walker[ic] = u * u;
  walker[im] = m;
  walker[ix] = beta * m * m;
}


void Accumulate(void){ //Update block averages
   for(int i=0; i<n_props; ++i){
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){ //Print results for current block

   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;

    //cout << "Block number " << iblk << endl;
    //cout << "Acceptance rate " << accepted/attempted << endl << endl;

    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    //Ene.open("../../Esercizio6/Risultati/equilibration_en_gibbs.dat",ios::app);
    //Ene << iblk <<  setw(wd) << glob_av[iu]/(double)iblk << setw(wd) <<  err_u << endl;
    //Ene.close();
    if(iblk == nblk){
      cout << "Energy: " << endl;
      cout << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      Ene.open("../../Esercizio6/Risultati/ene_"+to_string(metro)+".0",ios::app);
      Ene << temp << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      Ene.close();
    }

    //stima_u_squared = blk_av[ic]/blk_norm/(double)nspin; //Energy squared
    stima_heat = beta*beta*(blk_av[ic]/blk_norm - (blk_av[iu]/blk_norm)*(blk_av[iu]/blk_norm))/(double)nspin; ;
    glob_av[ic]  += stima_heat;
    glob_av2[ic] += stima_heat*stima_heat;
    err_heat=Error(glob_av[ic],glob_av2[ic],iblk);
    if(iblk == nblk){
      cout << "Heat capacity:" << endl;
      cout << setw(wd) << iblk <<  setw(wd) << stima_heat << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_heat << endl;
      Heat.open("../../Esercizio6/Risultati/heat_"+to_string(metro)+".0",ios::app);
      Heat << temp << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_heat << endl;
      Heat.close();
    }
/*
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetisation
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);

    if(iblk == nblk){
      cout << "Magnetisation:" << endl;
      cout << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      Mag.open("../../Esercizio6/Risultati/m_"+to_string(metro)+".0",ios::app);
      Mag << temp << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      Mag.close();
    }
*/
    stima_x = blk_av[ix]/blk_norm/(double)nspin; // Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    if(iblk==nblk){
      cout << "Magnetic susceptibility:" << endl;
      cout << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      Chi.open("../../Esercizio6/Risultati/chi_"+to_string(metro)+".0",ios::app);
      Chi << temp << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      Chi.close();
    }
}


void ConfFinal(void){
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i){
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i){  //Algorithm for periodic boundary conditions
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk){
    if(iblk==1) return 0.0;
    else{
      double appo = sum2/(double)iblk - pow(sum/(double)iblk,2);
      if(appo >= 0) return sqrt(appo/(double)(iblk-1));
      else return sqrt(-appo/(double)(iblk-1));
    }
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
