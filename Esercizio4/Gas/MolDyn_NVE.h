#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;

// averages
double acc,att;

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_kin,stima_etot,stima_temp;
double err_pot,err_kin,err_etot,err_temp;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;
bool equilibration;
int n_equilibration,rescaling_step;
std::vector<double> sumv_sqr;
int nblk,n_per_blk;
double walker[m_props];

//functions
void Input(void);
void ConfInit(void);
void Equilibration(void);
void Move(void);
void ConfFinal(void);
void Measure_on_blocks();
void Measure_Equilibration();
void Accumulate();
void Averages(int);
void Reset(int);
double Force(int, int);
double Pbc(double);
double Error(double, double, int);
