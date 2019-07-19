/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __vmc_
#define __vmc_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters
double mu, sigma;
double x0, x;
const int m_props=1000;
int n_props, ip, ie;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_ene,err_ene;

//histogram of |Psi|^2
int nbins;
double histo_xmax, histo_xmin, bin_size, N;

// simulation
int nstep, nblk;
double delta;
int restart, equilibration;

//functions
void Input(void);
double Psi_Squared(double y);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Equilibrate(void);
void Measure(void);
double Error(double,double,int);
void Histo(double y);
void Normalize_Histo(void);

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
