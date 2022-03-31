#include "MolDyn_NVE.h"

using namespace std;

int main(){
  Input();                                     //Inizialization
  ConfInit();                                  //Configuration inizialization
  Equilibration();                             //Equilibration
  for(int iblk=1; iblk<=nblk; iblk++){
    Reset(iblk);                               //Reset block averages
    for(int j=0; j<n_per_blk; j++){
      Move();                                  //Move particles with Verlet algorithm
      if(j%10 == 0){
        Measure_on_blocks();                        //Properties measurement
        Accumulate();                          //Update block averages
      }
    }
    Averages(iblk);                            //Print results for current block
    if(iblk%10==0) cout << "Block " << iblk << endl;
  }
  ConfFinal();                                 //Write final configuration to restart
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput;

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
  ReadInput >> n_per_blk;
  ReadInput >> n_equilibration;
  ReadInput >> rescaling_step;

  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);

  cout << "Number of particles = " << npart << endl;
  cout << "Density of particles = " << rho << endl;
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;
  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of measurements per blocks = " << n_per_blk << endl;
  ReadInput.close();

  //Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ConfInit(){
  //Read initial configuration
  ifstream ReadConf;
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
  //Prepare initial velocities
  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
  double sumv[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<npart; ++i){   //Non voglio moti di deriva per il mio sistema!!
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
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VelocitiesRescaling(){
  double sumv2_mean=0,fs;
  for(unsigned int i=0; i<sumv_sqr.size(); i++) sumv2_mean += sumv_sqr.at(i);
  sumv2_mean /= sumv_sqr.size();
  fs = sqrt(3 * temp / sumv2_mean);
  for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;
     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
  }
  sumv_sqr.clear();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Equilibration(){
  cout << "The system is being balanced.." << endl;
  temp = 1.0;
  for(int istep=1; istep<=n_equilibration; ++istep){
      Move();
      double sumv2 = 0.,v_x,v_y,v_z;
      for(int i=0; i<npart; i++){
        v_x = Pbc(x[i] - xold[i])/delta;
        v_y = Pbc(y[i] - yold[i])/delta;
        v_z = Pbc(z[i] - zold[i])/delta;
        sumv2 += v_x*v_x + v_y*v_y + v_z*v_z;
      }
      sumv2 /= (double)npart;
      sumv_sqr.push_back(sumv2);
      if(istep%rescaling_step==0) VelocitiesRescaling();
      if(istep%20) Measure_Equilibration();
      if(istep%2000==0 && temp!=1.2 && temp<1.2) temp +=0.02;
  }
  for(int istep=1; istep<=20*n_equilibration; ++istep){
    Move();
    if(istep%20) Measure_Equilibration();
  }
  cout << "The equilibration of the system is concluded." << endl;
  cout << "Let's start the simulation.." << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Measure_Equilibration(){ //Properties measurement
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("../../Esercizio4/Risultati/gas_equilibration_epot.dat",ios::app);
  Ekin.open("../../Esercizio4/Risultati/gas_equilibration_ekin.dat",ios::app);
  Temp.open("../../Esercizio4/Risultati/gas_equilibration_temp.dat",ios::app);
  Etot.open("../../Esercizio4/Risultati/gas_equilibration_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){
     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)
     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6); //Potential energy
       v += vij;
     }
    }
  }

  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);//Kinetic energy
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle

  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Measure_on_blocks(){ //Properties measurement
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  v = 0.0; //reset observables
  t = 0.0;
  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){
     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)
     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6); //Potential energy
       v += vij;
     }
    }
  }

  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);//Kinetic energy
  walker[iv] = v/(double)npart; //Potential energy per particle;
  walker[ik] = t/(double)npart; //Kinetic energy per particle
  walker[ie] = (t+v)/(double)npart; //Total energy per particle
  walker[it] = (2.0 / 3.0) * t/(double)npart; //Temperature
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Accumulate(void){ //Update block averages
   for(int i=0; i<n_props; ++i){
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Averages(int iblk){
  ofstream Epot, Ekin, Etot, Temp; //stima_pot, stima_kin, stima_etot, stima_temp;
  const int wd=12;

  Epot.open("../../Esercizio4/Risultati/gas_epot.dat",ios::app);
  Ekin.open("../../Esercizio4/Risultati/gas_ekin.dat",ios::app);
  Temp.open("../../Esercizio4/Risultati/gas_temp.dat",ios::app);
  Etot.open("../../Esercizio4/Risultati/gas_etot.dat",ios::app);

  stima_pot = blk_av[iv]/blk_norm; //Potential energy
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
  Epot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;

  stima_kin = blk_av[ik]/blk_norm; //Kinetic energy
  glob_av[ik] += stima_kin;
  glob_av2[ik] += stima_kin*stima_kin;
  err_kin=Error(glob_av[ik],glob_av2[ik],iblk);
  Ekin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin  << endl;

  stima_etot = blk_av[ie]/blk_norm; //Total energy
  glob_av[ie] += stima_etot;
  glob_av2[ie] += stima_etot*stima_etot;
  err_etot=Error(glob_av[ie],glob_av2[ie],iblk);
  Etot <<  "   " << setprecision(8) << glob_av[ie]/(double)iblk <<  "   " <<  setprecision(8) << err_etot << endl;

  stima_temp = blk_av[it]/blk_norm; //Temperature
  glob_av[it] += stima_temp;
  glob_av2[it] += stima_temp*stima_temp;
  err_temp=Error(glob_av[it],glob_av2[it],iblk);
  Temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<npart; ++i)  WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  WriteConf.close();
  cout << "Print old configuration to file config.old " << endl << endl;
  WriteConf.open("config.old");
  for (int i=0; i<npart; ++i)  WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  WriteConf.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Error(double sum, double sum2, int iblk){
    if(iblk==1) return 0.0;
    else{
      double appo = sum2/(double)iblk - pow(sum/(double)iblk,2);
      if(appo >= 0) return sqrt(appo/(double)(iblk-1));
      else return sqrt(-appo/(double)(iblk-1));
    }
}
