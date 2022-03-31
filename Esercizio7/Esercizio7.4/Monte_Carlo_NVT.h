/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __NVT__
#define __NVT__
#include "random.h"
#include <string>

using namespace std;

// Random numbers
int seed[4];
Random rnd;

// parameters, observables
const int m_props=1000;
int n_props, iv, iw, igofr;
double vtail, ptail, bin_size, nbins, sd;
double walker[m_props];

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_pres, err_pot, err_press, err_gdir;
double stima_accettazione;

//configuration
const int m_part=108;
double x[m_part], y[m_part], z[m_part];

// thermodynamical state
int npart;
double beta, temp, vol, rho, box, rcut;
string phase;

// simulation
int nstep, nblk, restart;
double delta;

// pigreco
const double pi=3.1415927;

// functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
