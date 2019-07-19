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
#include <cstdlib>
#include "random.h"
#include "main.h"

using namespace std;
 
int main (int argc, char *argv[]) {

	Input();
	ofstream coordinates, averages;
	
	Initialize(0., 0., 0.);
	//1s equilibration, uniform mode
	if (equilibration == 1) { 
		cout << "Equilibration for 1s orbital (uniform sampling)" << endl;
		delta_1s_uniform = Equilibrate_1s("uniform", delta_1s_uniform);
		cout << "Equilibrated! Use delta_1s_uniform = " << delta_1s_uniform << endl << endl;
	}
	//1s simulation, uniform mode
	if (equilibration == 0) { 
		cout << "Simulation of 1s Hydorgen Atom Orbital (uniform sampling)" << endl;
		coordinates.open("coordinates_1s_uniform.txt");
		averages.open("averages_1s_uniform.txt");
		for (int i=1; i<=n_blocks; i++) {
			prog_sum = 0.;
			for (int j=0; j<block_length; j++) {
				coordinates << xold << "	" << yold << "	" << zold << endl;
				Simulate_1s("uniform", delta_1s_uniform);
				Metropolis();
				Accumulate(r(xnew,ynew,znew));
			}
			Block_Average(i);
			averages << i << "	" << average << "	" << error << endl;		
		}
		coordinates.close();
		averages.close();
	}
	
	Initialize(0., 0., 0.);
	//1s equilibration, gaussian mode
	if (equilibration == 1) { 
		cout << "Equilibration for 1s orbital (gaussian sampling)" << endl;
		delta_1s_gaussian = Equilibrate_1s("gaussian", delta_1s_gaussian);
		cout << "Equilibrated! Use delta_1s_gaussian = " << delta_1s_gaussian << endl << endl;
	}
	//1s simulation, gaussian mode
	if (equilibration == 0) { 
		cout << "Simulation of 1s Hydorgen Atom Orbital (gaussian sampling)" << endl;
		coordinates.open("coordinates_1s_gaussian.txt");
		averages.open("averages_1s_gaussian.txt");
		for (int i=1; i<=n_blocks; i++) {
			prog_sum = 0.;
			for (int j=0; j<block_length; j++) {
				coordinates << xold << "	" << yold << "	" << zold << endl;
				Simulate_1s("gaussian", delta_1s_gaussian);
				Metropolis();
				Accumulate(r(xnew,ynew,znew));
			}
			Block_Average(i);
			averages << i << "	" << average << "	" << error << endl;		
		}
		coordinates.close();
		averages.close();
	}

	Initialize(0., 0., 1.);
	//2p equilibration, uniform mode
	if (equilibration == 1) { 
		cout << "Equilibration for 2p orbital (uniform sampling)" << endl;
		delta_2p_uniform = Equilibrate_2p("uniform", delta_2p_uniform);
		cout << "Equilibrated! Use delta_2p_uniform = " << delta_2p_uniform << endl << endl;
	}
	//2p simulation, uniform mode
	if (equilibration == 0) { 
		cout << "Simulation of 2p Hydorgen Atom Orbital (uniform sampling)" << endl;
		coordinates.open("coordinates_2p_uniform.txt");
		averages.open("averages_2p_uniform.txt");
		for (int i=1; i<=n_blocks; i++) {
			prog_sum = 0.;
			for (int j=0; j<block_length; j++) {
				coordinates << xold << "	" << yold << "	" << zold << endl;
				Simulate_2p("uniform", delta_2p_uniform);
				Metropolis();
				Accumulate(r(xnew,ynew,znew));
			}
			Block_Average(i);
			averages << i << "	" << average << "	" << error << endl;		
		}
		coordinates.close();
		averages.close();
	}
		
	Initialize(0., 0., 1.);
	//2p equilibration, gaussian mode
	if (equilibration == 1) { 
		cout << "Equilibration for 2p orbital (gaussian sampling)" << endl;
		delta_2p_gaussian = Equilibrate_2p("gaussian", delta_2p_gaussian);
		cout << "Equilibrated! Use delta_2p_gaussian = " << delta_2p_gaussian << endl << endl;
	}
	//2p simulation, gaussian mode
	if (equilibration == 0) { 
		cout << "Simulation of 2p Hydorgen Atom Orbital (gaussian sampling)" << endl;
		coordinates.open("coordinates_2p_gaussian.txt");
		averages.open("averages_2p_gaussian.txt");
		for (int i=1; i<=n_blocks; i++) {
			prog_sum = 0.;
			for (int j=0; j<block_length; j++) {
				coordinates << xold << "	" << yold << "	" << zold << endl;
				Simulate_2p("gaussian", delta_2p_gaussian);
				Metropolis();
				Accumulate(r(xnew,ynew,znew));
			}
			Block_Average(i);
			averages << i << "	" << average << "	" << error << endl;		
		}
		coordinates.close();
		averages.close();
	}
	
	rnd.SaveSeed();
	return 0;
}

void Input() {
	
	ifstream ReadInput;
	
	cout << endl << "Hydrogen Atom Orbitals" << endl << endl;
	
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
	ReadInput >> n_blocks;
	ReadInput >> delta_1s_uniform;
	ReadInput >> delta_1s_gaussian;
	ReadInput >> delta_2p_uniform;
	ReadInput >> delta_2p_gaussian;
	ReadInput >> equilibration;
	
	block_length = n_throws/n_blocks;

}

void Initialize(double x0, double y0, double z0) {
		
	prog_sum = 0.;
	prog_blkaverage = 0.;
	prog_squared_blkaverage = 0.;
	xold = x0;
	yold = y0;
	zold = z0;
	n_accepted = 0;	

}
	

void Accumulate(double x) {
	
	prog_sum += x;

}

void Block_Average(int iblock) {

	blkaverage = prog_sum/(double)block_length;
	squared_blkaverage = blkaverage*blkaverage;
	
	prog_blkaverage += blkaverage;
	prog_squared_blkaverage += squared_blkaverage;
	
	average = 1./(double)(iblock) * prog_blkaverage;
	squared_average = 1./(double)(iblock) * prog_squared_blkaverage;
	
	if (iblock == 1) error = 0.;
	else error = sqrt(1./(double)(iblock-1) * (squared_average - average*average));
	
}

void Simulate_1s(string mode, double delta) {

	prob_old = gs(xold,yold,zold);
	if (mode == "uniform") {
		xnew = rnd.Rannyu(xold-delta, xold+delta);
		ynew = rnd.Rannyu(yold-delta, yold+delta);
		znew = rnd.Rannyu(zold-delta, zold+delta);
	}
	if (mode == "gaussian") {
		xnew = rnd.Gauss(xold, delta);
		ynew = rnd.Gauss(yold, delta);
		znew = rnd.Gauss(zold, delta);
	}
	prob_new = gs(xnew,ynew,znew);

}

void Simulate_2p(string mode, double delta) {

	prob_old = es(xold,yold,zold);
	if (mode == "uniform") {
		xnew = rnd.Rannyu(xold-delta, xold+delta);
		ynew = rnd.Rannyu(yold-delta, yold+delta);
		znew = rnd.Rannyu(zold-delta, zold+delta);
	}
	if (mode == "gaussian") {
		xnew = rnd.Gauss(xold, delta);
		ynew = rnd.Gauss(yold, delta);
		znew = rnd.Gauss(zold, delta);
	}
	prob_new = es(xnew,ynew,znew);

}

void Metropolis() {

	probability = min(1., prob_new/prob_old);
	double s = rnd.Rannyu();
	if (s<probability) {
		xold = xnew;
		yold = ynew;
		zold = znew;
		n_accepted++;
	}
	else {
		xnew = xold;
		ynew = yold;
		znew = zold;
	}
	
}

double Equilibrate_1s(string mode, double delta) {
	
	do {
		cout << "Equilibrating..." << endl;
		delta += 0.1;
		n_accepted = 0;
		for (int i=0; i<n_throws; i++) {
			Simulate_1s(mode, delta);
			Metropolis();
		}
		acceptance = (double)(n_accepted)/(double)(n_throws);}
	while (acceptance < 0.47 or acceptance > 0.53);
	
	return delta;

}

double Equilibrate_2p(string mode, double delta) {
	
	do {
		cout << "Equilibrating..." << endl;
		delta += 0.1;
		n_accepted = 0;
		for (int i=0; i<n_throws; i++) {
			Simulate_2p(mode, delta);
			Metropolis();
		}
		acceptance = (double)(n_accepted)/(double)(n_throws);}
	while (acceptance < 0.47 or acceptance > 0.53);
	
	return delta;
}


double r (double x, double y, double z) {
	return sqrt(x*x + y*y + z*z);
}


double theta (double x, double y, double z) {
	return acos(z/sqrt(x*x + y*y + z*z));
}


double gs (double x, double y, double z) {
	double wf = exp(-r(x,y,z))/sqrt(M_PI);
	return wf*wf;
}


double es (double x, double y, double z) {
	double wf = 1./8.*sqrt(2./M_PI)*r(x,y,z)*exp(-r(x,y,z)/2.)*cos(theta(x,y,z));
	return wf*wf;
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
