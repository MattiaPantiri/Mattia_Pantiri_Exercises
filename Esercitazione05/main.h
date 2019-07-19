/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <vector>

using namespace std;

//variables, parameters
const double a0 = 0.0529E-9;
double delta_1s_uniform, delta_1s_gaussian, delta_2p_uniform, delta_2p_gaussian;
double xold, yold, zold, xnew, ynew, znew;
double prob_old, prob_new;
int n_accepted;
double probability;
double acceptance;
int equilibration;

//block average
int n_throws, n_blocks, block_length;
double prog_sum, prog_blkaverage, prog_squared_blkaverage;
double blkaverage, squared_blkaverage;
double average, squared_average;
double error;

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//functions
void Input();
void Initialize(double x0, double y0, double z0);
void Accumulate(double x);
void Block_Average(int iblock);
void Simulate_1s(string mode, double delta);
void Simulate_2p(string mode, double delta);
double Equilibrate_1s(string mode, double delta);
double Equilibrate_2p(string mode, double delta);
void Metropolis();
double r (double x, double y, double z);
double theta (double x, double y, double z);
double gs (double x, double y, double z);
double es (double x, double y, double z);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
