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
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Variational_MC.h"

using namespace std;

int main() { 
	Input(); //Inizialization
	if (equilibration == 1) Equilibrate();
	else {
		for(int iblk=1; iblk <= nblk; ++iblk) { //Simulation
			Reset(iblk);   //Reset block averages
			for(int istep=1; istep <= nstep; ++istep) {
				Move();
				Measure();
				Accumulate(); //Update block averages
			}
			Averages(iblk);   //Print results for current block
		}
	}
	rnd.SaveSeed();
	return 0;
}


void Input(void)	{

	system("rm *.txt");		
	
	ifstream ReadInput, ReadDelta;

	cout << "1D Variational Monte Carlo        " << endl;
	cout << "Particle in a potential V(x) = x^4 - 5/2 x^2             " << endl << endl;
	cout << "Trial wave function of the form Psi(x) = exp(-((x-mu)^2)/(2sigma^2)) + exp(-((x+mu)^2)/(2sigma^2))" << endl;
	cout << "Parameters to optimize: mu, sigma" << endl << endl;
	
	ReadInput.open("input.dat");
	ReadInput >> equilibration;

	
//Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
  
//Read input informations

	ReadInput >> mu;
	cout << "Starting with mu = " << mu << endl;

	ReadInput >> sigma;
	cout << "Starting with sigma = " << sigma << endl;
	if (equilibration!=1) {
		ReadDelta.open("delta.dat");
		ReadDelta >> delta;
		ReadDelta.close();
		cout << "The program perform Metropolis moves with uniform translations" << endl;
		cout << "Moves parameter = " << delta << endl;
	}
	else delta = 0.1;
	
	ReadInput >> x0;
	cout << "Starting position = " << x0 << endl;
	x = x0;
	
	ReadInput >> histo_xmin;
	ReadInput >> histo_xmax;
	ReadInput >> nbins;
	
	bin_size = (histo_xmax-histo_xmin)/(double)(nbins);

	ReadInput >> nblk;
	ReadInput >> nstep;

	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;
	ReadInput.close();
	
	if (equilibration == 1) cout << "Equilibration Phase" << endl;

//Prepare arrays for measurements
	ie = 0; //Energy
	n_props = 1;
	
	ip = 1; //|Psi|^2
	n_props = n_props + nbins; //Number of observables
	
	for (int i=0; i<n_props; i++) walker[i] = 0.0;

}

void Move(void) {
	double p, prob_old, prob_new;
	double xold, xnew;

	//Old
	xold = x;

	prob_old = Psi_Squared(xold);

	//New
	xnew = x + delta*(rnd.Rannyu() - 0.5);

	prob_new = Psi_Squared(xnew);

	//Metropolis test
	p = prob_new/prob_old;
	if(p >= rnd.Rannyu()) {
		//Update
		x = xnew;
    
		accepted = accepted + 1.0;
	}
	attempted = attempted + 1.0;
	
}

double Psi_Squared(double y) {

	double psi = exp(-((y-mu)*(y-mu))/(2*sigma*sigma)) + exp(-((y+mu)*(y+mu))/(2*sigma*sigma));
	return psi*psi;
	
}


void Reset(int iblk) { //Reset block averages

	if(iblk == 1) {
		for(int i=0; i<n_props; ++i) {
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i) {
		blk_av[i] = 0;
	}
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}

void Accumulate(void) { //Update block averages

	for(int i=0; i<n_props; ++i) {
		blk_av[i] = blk_av[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) { //Print results for current block

	ofstream Psi, Ene;
	const int wd=20;
    
	cout << "Block number " << iblk << endl;
	cout << "Acceptance rate " << accepted/attempted << endl << endl;
	
	stima_ene = blk_av[ie]/blk_norm; //Potential energy
	glob_av[ie] += stima_ene;
	glob_av2[ie] += stima_ene*stima_ene;
	err_ene=Error(glob_av[ie],glob_av2[ie],iblk);
    

    Ene.open("energy.txt",ios::app);
	
	
	//Energy
	Ene << iblk << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_ene << endl;
	
	if (iblk == nblk) { //Psi Squared
		Psi.open("psi_squared.txt");
		Normalize_Histo();
		for (int i=ip; i<n_props; i++) {
			Psi << histo_xmin + bin_size*(2*(i-ip)+1)/2.0 << setw(wd) << walker[i] << endl;
		}
		Psi.close();
	}

	Ene.close();
	cout << "----------------------------" << endl << endl;
}

void Measure() {
	double v, k;
	
	v = pow(x, 4) - 2.5*x*x;
	k = (-0.5/(sigma*sigma))*((-1+(x-mu)*(x-mu)/(sigma*sigma))*exp(-((x-mu)*(x-mu))/(2*sigma*sigma)) + (-1+(x+mu)*(x+mu)/(sigma*sigma))*exp(-((x+mu)*(x+mu))/(2*sigma*sigma)));


	walker[ie] = (k+v*sqrt(Psi_Squared(x)))/sqrt(Psi_Squared(x));
	Histo(x);
}

void Equilibrate(void) { //equilibrates delta to obtain Metropolis acceptance 50% and prints the obtained value in delta.dat
	
	double acceptance;
	ofstream Del;
	do {
		cout << "Equilibrating..." << endl;
		delta += 0.1;
		accepted = 0;
		attempted = 0;
		for (int i=0; i<nblk; i++) Move();
		acceptance = accepted/attempted;
		cout << "acceptance rate: " << acceptance << endl;}
	while (acceptance < 0.47 or acceptance > 0.53);
	cout << "Equilibrated! Correct value of delta printed in file delta.dat" << endl;
	Del.open("delta.dat");
	Del << delta << endl;
	Del.close();
}

double Error(double sum, double sum2, int iblk) {
	return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Histo(double y) {

	int pos = (int)((y-histo_xmin)/bin_size);
	walker[ip+pos]++;

}

void Normalize_Histo(void) {

	double N = 0;
	for (int i=ip; i<n_props; i++) N += bin_size*walker[i];
	for (int i=ip; i<n_props; i++) walker[i] /= N;
	
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
