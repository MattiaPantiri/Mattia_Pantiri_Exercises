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
	output.open("discreteRW.txt");
	for (int j=1; j<=n_simulations; j++) {
		RW* randw = new RW();
		Simulate_DiscreteRW(randw);
		if (j == block_length*h) {
			Block_Average(h);
			for (int i=0; i<n_steps; i++) prog_sum[i] = 0.;
			h++;
		}
		delete randw; 
	}
	cout << "1) Discrete random walk with step length 1: " << endl;
	Print_FinalRW();	
	output.close();
	
	Initialize();	
	output.open("continuosRW.txt");
	for (int j=1; j<=n_simulations; j++) {
		RW* randw = new RW();
		Simulate_ContinuosRW(randw);
		if (j == block_length*h) {
			Block_Average(h);
			for (int i=0; i<n_steps; i++) prog_sum[i] = 0.;
			h++;
		}
		delete randw;
	}
	cout << "2) Continuos random walk with step length 1: " << endl;
	Print_FinalRW();	
	output.close();
	
	rnd.SaveSeed();
	return 0;
}

void Input() {
	
	ifstream ReadInput;
	
	cout << endl << "Random Walk simulation, starting from the origin" << endl << endl;
	
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
	
	ReadInput >> n_simulations;
	cout << "--> Number of simulations = " << n_simulations << endl;
	
	ReadInput >> n_steps;  
	cout << "--> Number of steps for each random walk = " << n_steps << endl << endl;	
	
	ReadInput >> n_blocks;
	
	block_length = n_simulations/n_blocks;
	
	prog_sum.resize(n_steps);
	prog_blkaverage.resize(n_steps);
	prog_squared_blkaverage.resize(n_steps);
	blkaverage.resize(n_steps);
	squared_blkaverage.resize(n_steps);
	average.resize(n_steps);
	squared_average.resize(n_steps);
	error.resize(n_steps);

	ReadInput.close();
	
}

void Initialize() {
	
	h = 1;	
	for (int i=0; i<n_steps; i++) {
		prog_sum[i] = 0.;
		prog_blkaverage[i] = 0.;
		prog_squared_blkaverage[i] = 0.;
	}

}

void Block_Average(int h) {
	
	for (int k=0; k<n_steps; k++) {
	
		blkaverage[k] = prog_sum[k]/(double)block_length;
		squared_blkaverage[k] = blkaverage[k]*blkaverage[k];
	
		prog_blkaverage[k] += blkaverage[k];
		prog_squared_blkaverage[k] += squared_blkaverage[k];
	
		average[k] = 1./(double)h * prog_blkaverage[k];
		squared_average[k] = 1/(double)h * prog_squared_blkaverage[k];
	
		if (h == 1) error[k] = 0.;
		else error[k] = sqrt(1/(double)(h-1) * (squared_average[k] - average[k]*average[k]));	
		
	}	
	
}

void Print_FinalRW() {
	
	cout << "Step	r	Error" << endl;
	for (int i=0; i<n_steps; i++) {
		output << i+1 << "	" << sqrt(average[i]) << "	" << (1./(2*sqrt(average[i])))*error[i] << endl;
		cout << i+1 << "	" << sqrt(average[i]) << "	" << (1./(2*sqrt(average[i])))*error[i] << endl;
	}
	cout << endl;

}

void Simulate_DiscreteRW(RW* randw) {
	
	for (int i=0; i<n_steps; i++) {
		r = rnd.Rannyu(0., 6.);
		if((int)r == 0) {randw -> incrX();}
		if((int)r == 1) {randw -> decrX();}
		if((int)r == 2) {randw -> incrY();}
		if((int)r == 3) {randw -> decrY();}
		if((int)r == 4) {randw -> incrZ();}
		if((int)r == 5) {randw -> decrZ();}
		prog_sum[i] += randw -> Square_Modulus();	
	}

}

void Simulate_ContinuosRW(RW* randw) {

	for (int i=0; i<n_steps; i++) {
		theta = rnd.Theta();
		phi = rnd.Rannyu(0., 2*M_PI);
		randw -> Spherical(1., theta, phi);
		prog_sum[i] += randw -> Square_Modulus();	
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
