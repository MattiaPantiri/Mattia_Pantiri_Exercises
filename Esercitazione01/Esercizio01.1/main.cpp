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
#include "main.h"

using namespace std;
 
int main (int argc, char *argv[]){
	
	Input();
	
	Initialize();
	output.open("mean.txt");
	for (int i=0; i<n_blocks; i++) {
		prog_sum = 0.;
		for (int j=0; j<block_length; j++) {
			r = rnd.Rannyu();
			Accumulate(r);
		}
		Block_Average(i);
		if (n == n_blocks) cout << "1) Mean value of random number generator at the end of the test:" << endl;
		Print();
	}
	output.close();
	
	Initialize();
	output.open("variance.txt");
	for (int i=0; i<n_blocks; i++) {
		prog_sum = 0.;
		for (int j=0; j<block_length; j++) {
			r = rnd.Rannyu();
			Accumulate((r-0.5)*(r-0.5));
		}
		Block_Average(i);
		if (n == n_blocks) cout << "2) Variance of random number generator at the end of the test:" << endl;
		Print();
	}
	output.close();
	
	output.open("chi.txt");
	prog_sum = 0.;	
	for (int i=0; i<n_chitest; i++) {
		Chi_Initialize();
		Chi_Fill();
		for (int j=0; j<n_subintervals; j++) chi += (double)(pow(chi_counter[j] - n_chithrows/n_subintervals, 2))/(n_chithrows/n_subintervals);
		Accumulate(chi);
		output << i+1 << "	" << chi << endl;
		if (i == n_chitest-1) {
			chi_mean = prog_sum/(double)n_chitest;
			cout << "3) Mean value of chi-square between all chi-square tests: " << chi_mean << endl << endl;
		}
	}
	output.close();
	
	return 0;
}


void Input() {
	
	ifstream ReadInput;
	
	cout << endl << "Random Number Generator Test" << endl << endl;
	
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
	
	cout << "--> Seed used: " << seed[0] << "  " << seed[1] << "  " << seed[2] << "  " << seed[3] << endl;
	
	ReadInput >> n_throws;
	cout << "--> Number of throws = " << n_throws << endl;
	
	ReadInput >> n_blocks;  
	cout << "--> Number of blocks = " << n_blocks << endl;
	
	ReadInput >> n_chitest;
	cout << "--> Number of chi-square tests = " << n_chitest << endl;
	
	ReadInput >> n_subintervals;
	cout << "--> Number of sub-intervals for chi-square test = " << n_subintervals << endl;
	
	ReadInput >> n_chithrows;
	cout << "--> Number of throws for each chi-square test = " << n_chithrows << endl << endl;
	
	block_length = n_throws/n_blocks;
	chi_bin = 1./(double)n_subintervals;	

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

void Chi_Initialize() {
	
	chi_counter.resize(n_blocks);		
	for (int j=0; j<n_blocks; j++) chi_counter[j] = 0;
	chi = 0.;
	
}	

void Chi_Fill() {
	
	for (int j=0; j<n_chithrows; j++) { 
		r = rnd.Rannyu();
		for (int k=0; k<n_subintervals; k++) {
			if (r > k*chi_bin && r < (k+1)*chi_bin) chi_counter[k]++;
		}			
	}

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
