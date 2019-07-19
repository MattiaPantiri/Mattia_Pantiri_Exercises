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
int n_throws, n_bins;
double r;
double progSum;
int N[4] = {1, 2, 10, 100};
vector<double> sum;

//random numbers
#include "random.h"
int seed[4];
Random rnd;

//output stream
ofstream output;

//functions
void Input();
void Histo_Throw(string mode, int k);
void Histo_Bin(int bins, double xmin, double xmax);


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
