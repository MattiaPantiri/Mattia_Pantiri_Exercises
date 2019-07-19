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
#include "random.h"
#include "main.h"

using namespace std;
 
int main (int argc, char *argv[]){
	
	Input();
		
	for (int k=0; k<4; k++) {
		if (k==0) output.open("standard1.txt");
		if (k==1) output.open("standard2.txt");
		if (k==2) output.open("standard10.txt");
		if (k==3) output.open("standard100.txt");
		Histo_Throw("standard", k);
		Histo_Bin(n_bins, 0., 1.);
		output.close();		
	}
	
	
	for (int k=0; k<4; k++) {
		if (k==0) output.open("exponential1.txt");
		if (k==1) output.open("exponential2.txt");
		if (k==2) output.open("exponential10.txt");
		if (k==3) output.open("exponential100.txt");
		Histo_Throw("exponential", k);
		Histo_Bin(n_bins, 0., 6.);
		output.close();		
	}
	
	
	for (int k=0; k<4; k++) {
		if (k==0) output.open("cauchy_lorentz1.txt");
		if (k==1) output.open("cauchy_lorentz2.txt");
		if (k==2) output.open("cauchy_lorentz10.txt");
		if (k==3) output.open("cauchy_lorentz100.txt");
		Histo_Throw("cauchy_lorentz", k);
		Histo_Bin(n_bins, -4., 4.);
		output.close();	
	}
		
	//rnd.SaveSeed();
	return 0;
}


void Input() {
	
	ifstream ReadInput;
	
	cout << endl << "Random Number Generator extension & Central Limit Theorem test" << endl;
	
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
	cout << "--> Number of throws per histogram = " << n_throws << endl;
	
	ReadInput >> n_bins;
	cout << "--> Number of bins per histogram = " << n_bins << endl;
	
	sum.resize(n_throws);
	
	ReadInput.close();

}

void Histo_Throw(string mode, int k) {

	for (int i=0; i<n_throws; i++) {
		progSum = 0.;
		for (int j=0; j<N[k]; j++) {
			if (mode == "standard") r = rnd.Rannyu();
			if (mode == "exponential") r = rnd.Exponential(1.);
			if (mode == "cauchy_lorentz") r = rnd.CauchyLorentz(0., 1.);
			progSum += r;
		}
		sum[i] = 1./(double)N[k] * progSum;
	}
		
}

void Histo_Bin(int bins, double xmin, double xmax) {

	vector<int> bin_vector;
	bin_vector.resize(bins);
	for (int i=0; i<bins; i++) bin_vector[i] = 0;
	double bin_size = (xmax-xmin)/(double)bins;
	
	for (int i=0; i<n_throws; i++) {
		for (int k=0; k<bins; k++) {
			if (sum[i] >= xmin+k*bin_size && sum[i] < xmin+(k+1)*bin_size) bin_vector[k]++;
		}
	}
	
	for (int i=0; i<bins; i++) {
		output << xmin+i*bin_size << "	" << bin_vector[i] << endl;
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
