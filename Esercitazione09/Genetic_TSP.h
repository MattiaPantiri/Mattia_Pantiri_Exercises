/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __TSP_
#define __TSP_

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
	void Shift();
	void Partial_Shift();
	void Permutation();
	void Inversion();
	
};
	

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//variables
int n_cities, n_elements, n_generations;
vector<Individual> Population;
vector<Individual> New_Population;
Individual Father;
Individual Mother;
Individual Son;
Individual Daughter;
double R, theta, L, ave;
vector<double> x;
vector<double> y;
vector<vector<double>> distances;

ofstream output;

//functions
void Input();
void Initialize(string mode);
bool Check();
bool mySort(Individual a, Individual b) {return a.GetLength() < b.GetLength();}
void SortPopulation();
Individual Select();
void Generate();
void GenerateNewPopulation();
void PrintBestPath(string filename);
double Average();
int Pbc(int n, int size);




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
