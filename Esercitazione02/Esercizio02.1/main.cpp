/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "main.h"


using namespace std;
 
int main (int argc, char *argv[]){
	
	Input();
	
	Initialize();
	output.open("uniform_sampling.txt");
	for (int i=0; i<n_blocks; i++) {
		prog_sum = 0.;
		for (int j=0; j<block_length; j++) {
			r = rnd.Rannyu();
			Accumulate(M_PI/2. * cos(M_PI/2. * r));
		}
		Block_Average(i);
		if (n == n_blocks) cout << "1) Uniform Sampling:" << endl;
		Print();
	}
	output.close();
	
	Initialize();
	output.open("importance_sampling.txt");
	for (int i=0; i<n_blocks; i++) {
		prog_sum = 0.;
		for (int j=0; j<block_length; j++) {
			r = rnd.ImportanceSampling();
			Accumulate(M_PI/2. * cos(M_PI/2. * r) / (-2*(r-1)));
		}
		Block_Average(i);
		if (n == n_blocks) cout << "2) Importance Sampling:" << endl;
		Print();
	}
	output.close();
	
	//rnd.SaveSeed();
	return 0;
}

void Input() {
	
	ifstream ReadInput;
	
	cout << endl << "Evaluation of the integral between 0 and 1 of pi/2 * cos(pi/2 * x) dx" << endl << endl;
	
	//Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
	
	//Read variables from input.dat and print
	ReadInput.open("input.dat");
	
	ReadInput >> n_throws;
	cout << "--> Number of throws = " << n_throws << endl;
	
	ReadInput >> n_blocks;  
	cout << "--> Number of blocks = " << n_blocks << endl << endl;
	
	block_length = n_throws/n_blocks;	

	ReadInput.close();
	
}


void Initialize() {
		
	prog_sum = 0.;
	prog_blkaverage = 0.;
	prog_squared_blkaverage = 0.;

}
	

void Accumulate(double x) {
	
	prog_sum += x;

}

void Block_Average(int i) {

	n = i+1;

	blkaverage = prog_sum/(double)block_length;
	squared_blkaverage = blkaverage*blkaverage;
	
	prog_blkaverage += blkaverage;
	prog_squared_blkaverage += squared_blkaverage;
	
	average = 1./(double)n * prog_blkaverage;
	squared_average = 1/(double)n * prog_squared_blkaverage;
	
	if (n == 1) error = 0.;
	else error = sqrt(1/(double)(n-1) * (squared_average - average*average));
	
}

void Print() {
	
	output << n << "	" << average << "	" << error << endl;
	if (n == n_blocks) cout << fixed << setprecision(6) << "	" << average << " with uncertainty " << error << endl << endl;

}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
