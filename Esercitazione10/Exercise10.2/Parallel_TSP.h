/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Parallel_
#define __Parallel_

#include <vector>

using namespace std;

class Individual {

private:
	vector<int> _path;
	double _length;
	
public:
	
	Individual();
	Individual(int size);
	
	vector<int> GetPath() {return _path;}
	void SetPathSize(int size) {_path.resize(size);}
	int GetCity(int a) {return _path[a];}
	void SetCity(int city, int position) {_path[position] = city;}
	void Shuffle() {random_shuffle(_path.begin(), _path.end());}
	void Print() {for (unsigned int i=0; i<_path.size(); i++) cout << _path[i] << "  " ; cout << endl;}
	void SetLength();
	double GetLength() {return _length;}
	
	//mutation functions
	void Pair_Permutation();
	void Partial_Shift();
	void Permutation();
	void Inversion();
	
};


//Random numbers
#include "random.h"
int seed[4];
Random rnd[4], s;

//Parallel computing
int Size, Rank;

//variables
Individual TSP;
Individual New_TSP;
int n_cities;
double Temp, Temp_in, beta, R, L, theta, Scale_Factor, n_temp;
double prob;
vector<double> x_circ;
vector<double> y_circ;
vector<double> x_square;
vector<double> y_square;
vector<double> lengths;
vector<vector<double>> distances;
double r;
int N_accept[4], N_tot[4];
double acc[4];

ofstream output;

//parallel
int best_rank;

//functions
void Input(int a);
void Initialize(string mode, int a);
void Metropolis();
void Print (string mode);

int Pbc(int n, int size);

double min_value (double* v, int size);
int pos_min (double*v, int size);









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
