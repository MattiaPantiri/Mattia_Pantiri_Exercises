/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

using namespace std;

//variables, parameters
int n_throws, n_blocks, block_length;
double r;
int n;
double prog_sum, prog_blkaverage, prog_squared_blkaverage;
double blkaverage, squared_blkaverage;
double average, squared_average;
double error;
double S0, T, K, mu, sigma;
double S_T;
double t[101];

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//output stream
ofstream output;

//functions
void Input();
void Initialize();
void Block_Average(int i);
void Print();
void Direct_Call();
void Direct_Put();
void Discrete_Call();
void Discrete_Put();

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
