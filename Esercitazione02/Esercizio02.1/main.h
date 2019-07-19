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

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//output stream
ofstream output;

//functions
void Input();
void Initialize();
void Accumulate(double x);
void Block_Average(int i);
void Print();

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
