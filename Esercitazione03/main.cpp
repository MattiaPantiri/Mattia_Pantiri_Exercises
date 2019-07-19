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
#include <algorithm>
#include <iomanip>
#include "random.h"
#include "main.h"

using namespace std;

int main (int argc, char *argv[]){
	
	Input();
	
	Initialize();
	output.open("call_direct.txt");
	for (int i=0; i<n_blocks; i++) {
		prog_sum = 0.;
		for (int j=0; j<block_length; j++) Direct_Call();
		Block_Average(i);
		if (n == n_blocks) cout << "1) Direct Call value at the end of the simulation:" << endl;
		Print();
	}
	output.close();
	
	Initialize();
	output.open("put_direct.txt");
	for (int i=0; i<n_blocks; i++) {
		prog_sum = 0.;
		for (int j=0; j<block_length; j++) Direct_Put();
		Block_Average(i);
		if (n == n_blocks) cout << "2) Direct Put value at the end of the simulation:" << endl;
		Print();
	}
	output.close();
	
	Initialize();
	output.open("call_discrete.txt");
	for (int i=0; i<n_blocks; i++) {
		prog_sum = 0.;
		for (int j=0; j<block_length; j++) Discrete_Call();
		Block_Average(i);
		if (n == n_blocks) cout << "3) Discrete Call value at the end of the simulation:" << endl;
		Print();
	}
	output.close();
	
	Initialize();
	output.open("put_discrete.txt");
	for (int i=0; i<n_blocks; i++) {
		prog_sum = 0.;
		for (int j=0; j<block_length; j++) Discrete_Put();
		Block_Average(i);
		if (n == n_blocks) cout << "4) Discrete Put value at the end of the simulation:" << endl;
		Print();
	}
	output.close();
	
	rnd.SaveSeed();
	return 0;
}

void Input() {
	
	ifstream ReadInput;
	
	cout << endl << "Call-Put option prices" << endl << endl;
	
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
	cout << "--> Number of blocks = " << n_blocks << endl;
	
	block_length = n_throws/n_blocks;
	
	S0 = 100.;
	T = 1.;
	K = 100.;
	mu = 0.1;
	sigma = 0.25;
	
	for (int i=0; i<=100; i++) t[i] = i*1./100.;
	
	
	ReadInput.close();
	
}

void Initialize() {
		
	prog_sum = 0.;
	prog_blkaverage = 0.;
	prog_squared_blkaverage = 0.;

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

void Direct_Call() {
	
	r = rnd.Gauss(0., T);
	S_T = S0*exp((mu-0.5*sigma*sigma)*T + sigma*r);
	prog_sum += exp(-mu*T)*max(S_T - K, 0.);
	
}

void Direct_Put() {
	
	r = rnd.Gauss(0., T);
	S_T = S0*exp((mu-0.5*sigma*sigma)*T + sigma*r);
	prog_sum += exp(-mu*T)*max(K - S_T, 0.);

}

void Discrete_Call() {

	S_T = S0;
	for (int i=0; i<100; i++) {
		r = rnd.Gauss(0., 1.);
		S_T = S_T*exp((mu-0.5*sigma*sigma)*(t[i+1]-t[i]) + sigma*r*sqrt(t[i+1]-t[i]));
		}
	prog_sum += exp(-mu*T)*max(S_T - K, 0.);

}

void Discrete_Put() {

	S_T = S0;
	for (int i=0; i<100; i++) {
		r = rnd.Gauss(0., 1.);
		S_T = S_T*exp((mu-0.5*sigma*sigma)*(t[i+1]-t[i]) + sigma*r*sqrt(t[i+1]-t[i]));
	}
	prog_sum += exp(-mu*T)*max(K - S_T, 0.);

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
