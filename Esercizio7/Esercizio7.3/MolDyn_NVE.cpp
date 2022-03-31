#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include  <iomanip>	// setprecision
#include "MolDyn_NVE.h"

using namespace std;

int main(){
  Input(); // Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk){ // Simulation
	   Reset(iblk); // Reset block averages
	    for(int istep=1; istep <= nstep; ++istep){
	       Move(); // Move particles with Verlet algorithm
	       Measure();      // Properties measurement
	       Accumulate();   // Update block averages
	    }
	    Averages(iblk); // Print results for current block
  }
  ConfFinal();         // Write final configuration to restart
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator

  ReadInput.open("input.dat"); //Read input
  ReadInput >> temp;
  ReadInput >> npart;
  ReadInput >> rho;
  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> restart;
  ReadInput >> phase;
  ReadInput.close();

  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);

  cout << "Number of particles = " << npart << endl;
  cout << "Density of particles = " << rho << endl;
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;
  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Phase = " << phase << endl << endl;

  //Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
  stima_accettazione=0;
  //measurement of g(r)
  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

  //Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
  // Read old configuration
  if(restart==1){
    cout << "Read initial configuration from file config.final " << endl;
    ReadConf.open("config.final");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();
	  cout << "Read old configuration from file old.final " << endl << endl;
	  ReadConf.open("old.final");
	  for (int i=0; i<npart; ++i){
	    ReadConf >> xold[i] >> yold[i] >> zold[i];
	    xold[i] = xold[i] * box;
	    yold[i] = yold[i] * box;
	    zold[i] = zold[i] * box;
	  }
	  ReadConf.close();
  }
  if(restart == 0){
	   //Prepare initial velocities
	  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	  double sumv[3] = {0.0, 0.0, 0.0};
	  for (int i=0; i<npart; ++i){
	    vx[i] = rand()/double(RAND_MAX) - 0.5;
	    vy[i] = rand()/double(RAND_MAX) - 0.5;
	    vz[i] = rand()/double(RAND_MAX) - 0.5;
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
	  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	  double sumv2 = 0.0, fs;
	  for (int i=0; i<npart; ++i){
	    vx[i] = vx[i] - sumv[0];
	    vy[i] = vy[i] - sumv[1];
	    vz[i] = vz[i] - sumv[2];
	    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	  }
	  sumv2 /= (double)npart;

	  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
	  for (int i=0; i<npart; ++i){
	    vx[i] *= fs;
	    vy[i] *= fs;
	    vz[i] *= fs;

	    xold[i] = Pbc(x[i] - vx[i] * delta);
	    yold[i] = Pbc(y[i] - vy[i] * delta);
	    zold[i] = Pbc(z[i] - vz[i] * delta);
	  }
   }else{
	   Move(); // r(t + dt), r(t), v(t)
	   double sumv2 = 0., fs;
	   double v_x, v_y, v_z;
	   for(int i=0; i<npart; i++){
		   v_x = Pbc(x[i] - xold[i])/delta;
		   v_y = Pbc(y[i] - yold[i])/delta;
		   v_z = Pbc(z[i] - zold[i])/delta;
		   sumv2 += v_x*v_x + v_y*v_y + v_z*v_z;
	   }
	   sumv2 /= (double)npart;
	   fs = sqrt(3*temp / sumv2);
	   for (int i=0; i<npart; ++i){
	     // v(t) -> v_s(t)
	     vx[i] *= fs;
	     vy[i] *= fs;
	     vz[i] *= fs;
	     // r_new(t)
	     xold[i] = Pbc(x[i] - vx[i] * delta);
	     yold[i] = Pbc(y[i] - vy[i] * delta);
	     zold[i] = Pbc(z[i] - vz[i] * delta);
	   }
   }
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
}


double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  return f;
}


void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;
/*
  Epot.open("output.epot_ist.0",ios::app);
  Ekin.open("output.ekin_ist.0",ios::app);
  Temp.open("output.temp_ist.0",ios::app);
  Etot.open("output.etot_ist.0",ios::app);
*/
  v = 0.0; //reset observables
  t = 0.0;
  //reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;
  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){
     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );
     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     //histogram of g(r)
     bin = int(dr/bin_size);
     if(bin < nbins) walker[igofr + bin] += 2;
     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       //Potential energy
       v += vij;
     }
    }
  }
/*
//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    walker[iv] = stima_pot;
    walker[ik] = stima_kin;
    walker[ie] = stima_etot;
    walker[it] = stima_temp;

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();*/
}


void Reset(int iblk){
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


void Accumulate(void){
   for(int i=0; i<n_props; ++i) {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){
  double r, gdir;
  ofstream Gofr, Gave, Epot, Pres;
  const int wd=12;

  stima_accettazione += accepted/attempted/nblk;

  //Epot.open("../../Esercizio7/Risultati/output_"+phase+".epot.0",ios::app);
  //Pres.open("../../Esercizio7/Risultati/output_"+phase+".pres.0",ios::app);
  Gofr.open("../../Esercizio7/Risultati/molDyn_"+phase+".gofr.0",ios::app);
  Gave.open("../../Esercizio7/Risultati/molDyn_"+phase+".gave.0",ios::app);
/*
  stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot=Error(glob_av[iv],glob_av2[iv],iblk);

  stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
  glob_av[iw] += stima_pres;
  glob_av2[iw] += stima_pres*stima_pres;
  err_press=Error(glob_av[iw],glob_av2[iw],iblk);
  //Potential energy per particle
  Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
  //Pressure
  Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_press << endl;
*/  //g(r)
  for(int i=0; i<nbins; i++){
    r = bin_size*i;
    double r_dr = r + bin_size;
    blk_av[igofr+i] /= rho * npart * 4.*M_PI/3.*(r_dr*r_dr*r_dr - r*r*r);
    gdir = blk_av[igofr + i] / blk_norm;
    glob_av[igofr + i] += gdir;
    glob_av2[igofr + i] += gdir * gdir;
    err_gdir = Error(glob_av[igofr + i], glob_av2[igofr + i], iblk);
    Gofr<< setw(wd) << iblk << setw(wd) << r <<setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_gdir << endl;
    if(iblk == nblk) Gave << r << setw(wd) << glob_av[igofr + i]/double(nblk)<< setw(wd) << err_gdir << endl;
  }
  cout << "----------------------------" << endl << endl;

  //Epot.close();
  //Pres.close();
  Gofr.close();
  Gave.close();
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf, WriteOld;

  cout << endl << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for(int i=0; i<npart; ++i)  WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  WriteConf.close();

  cout << "Print final old configuration at to file old.final " << endl << endl;
  WriteOld.open("old.final");
  for(int i=0; i<npart; ++i)  WriteOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl; // aggiunto io
  WriteOld.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk){
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
