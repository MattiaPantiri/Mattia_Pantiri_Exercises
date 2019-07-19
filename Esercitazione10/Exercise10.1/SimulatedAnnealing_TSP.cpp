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
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include "SimulatedAnnealing_TSP.h"

using namespace std;

int main () {
	
	Input();
	
	Initialize("circumference");
	
	//int N = 1;
	cout << "	Simulating..." << endl;
	output.open("IterationBP_Circle.txt");
	for (int k=0; k<n_temp; k++) {
		/*for (int i=0; i<4; i++) {
			N_accept[i] = 0;
			N_tot[i] = 0;
		}*/
		for (int j=0; j<k*100; j++) {
			Metropolis();
			//output << N << "	" << TSP.GetLength() << endl;
			//N++;
		}
		output << k+1 << "	" << TSP.GetLength() << endl;
		//for (int i=0; i<4; i++) acc[i] = (double)N_accept[i]/(double)N_tot[i];
		Temp=Temp*Scale_Factor;
		beta = 1./Temp;
	}
	output.close();
	cout << "	Completed!" << endl;
	cout << "	Best Path Length:" << TSP.GetLength() << endl << endl;
	Print("BestPathCircle.txt", "circumference");
	
	
	Initialize("square");
	
	//N = 1;
	cout << "	Simulating..." << endl;
	output.open("IterationBP_Square.txt");
	for (int k=0; k<n_temp; k++) {
		/*for (int i = 0; i<4; i++) {
			N_accept[i] = 0;
			N_tot[i] = 0;
		}*/
		for (int j=0; j<k*100; j++) {
			Metropolis();
			//output << N << "	" << TSP.GetLength() << endl;
			//N++;
		}
		output << k+1 << "	" << TSP.GetLength() << endl;
		//for (int i=0; i<4; i++) acc[i] = (double)N_accept[i]/(double)N_tot[i];
		Temp=Temp*Scale_Factor;
		beta = 1./Temp;
	}
	output.close();
	cout << "	Completed!" << endl;
	cout << "	Best Path Length:" << TSP.GetLength() << endl << endl;
	Print("BestPathSquare.txt", "square");
	
	
return 0;
}


// functions for class "Individual"

Individual::Individual() {
}

Individual::Individual(int size) {
	vector<int> appo(size); 
	_path = appo;
}

void Individual::SetLength() {

	_length = 0.;
	for (int i=0; i<n_cities; i++) {
		_length += distances[GetCity(Pbc(i, n_cities))][GetCity(Pbc(i+1, n_cities))];
	}
	
}

void Individual::Pair_Permutation() {

	int i = (int)rnd.Rannyu(0., n_cities);
	int j = (int)rnd.Rannyu(0., n_cities);
	while (j==i) j = (int)rnd.Rannyu(0., n_cities);
	
	swap(_path[i], _path[j]);	

}

void Individual::Partial_Shift() {
	
	int n = (int)rnd.Rannyu(1., n_cities/2);
	int m = (int)rnd.Rannyu(1., n_cities/2);
	int j = (int)rnd.Rannyu(0., n_cities);
	
	vector<int> appo(n+m);

	for (int i=0; i<n+m; i++) appo[Pbc(i+n, n+m)] = _path[Pbc(j+i, n_cities)];

	for (int i=0; i<n+m; i++) _path[Pbc(j+i, n_cities)] = appo[i];

}

void Individual::Permutation() {

	int n = (int)rnd.Rannyu(0., n_cities/2);
	int j = (int)rnd.Rannyu(0., n_cities/2);
	
	for (int i=j; i<j+n; i++) swap(_path[i], _path[Pbc(i+n, n_cities)]);
	
}

void Individual::Inversion() {

	int n = (int)rnd.Rannyu(0., n_cities);
	int j = (int)rnd.Rannyu(0., n_cities);
	
	vector<int> appo(n);
	
	for (int i=0; i<n; i++) appo[i] = _path[Pbc(j+n-i-1, n_cities)];
	
	for (int i=j; i<j+n; i++) _path[Pbc(i, n_cities)] = appo[i-j];
	
}


//other funcions

void Input() {

	ifstream ReadInput;

	//read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
	
	//read input information
	ReadInput.open("input.dat");
  
	ReadInput >> n_cities;
	ReadInput >> R;
	ReadInput >> L;
	ReadInput >> Temp_in;
	ReadInput >> n_temp;
	ReadInput >> Scale_Factor;
	
	ReadInput.close();
	
	cout << endl << "Travelling Salesman Problem (TSP) - Parallel Simulated Annealing" << endl << endl;
	cout << "--> Number of cities = " << n_cities << endl;
	cout << "--> Initial temperature = " << Temp_in << endl;
	cout << "--> The program visits " << n_temp << " temperatures" << endl;
	cout << "--> Temperature scales with a factor = " << Scale_Factor << endl << endl;
	
	x_circ.resize(n_cities);
	y_circ.resize(n_cities);
	x_square.resize(n_cities);
	y_square.resize(n_cities);
	
	for (int i=0; i<n_cities; i++) {
		theta = rnd.Rannyu(0., 2*M_PI);
		x_circ[i] = R*cos(theta);
		y_circ[i] = R*sin(theta);
	}
	
	for (int i=0; i<n_cities; i++) {
		x_square[i] = rnd.Rannyu(-L, L);
		y_square[i] = rnd.Rannyu(-L, L);
	}
	
}

void Initialize(string mode) {
	
	distances.resize(n_cities);
	for (int i=0; i<n_cities; i++) distances[i].resize(n_cities);
	
	
	if (mode == "circumference") cout << "1) The cities are set randomly on a circumference of radius " << R << endl << endl;
	if (mode == "square") cout << "2) The cities are set randomly inside a square of side " << 2*L << endl << endl;
	cout << "	Computing distances between cities..." << endl;
	
	
	if (mode == "circumference") {
		for (int i=0; i<n_cities; i++) {
			for (int j=i; j<n_cities; j++) {
				distances[i][j] = sqrt((x_circ[i] - x_circ[j])*(x_circ[i] - x_circ[j]) + (y_circ[i] - y_circ[j])*(y_circ[i] - y_circ[j]));
				distances[j][i] = distances[i][j];
			}
		}
	}
	
	if (mode == "square") {
		for (int i=0; i<n_cities; i++) {
			for (int j=i; j<n_cities; j++) {
				distances[i][j] = sqrt((x_square[i] - x_square[j])*(x_square[i] - x_square[j]) + (y_square[i] - y_square[j])*(y_square[i] - y_square[j]));
				distances[j][i] = distances[i][j];
			}
		}
	}
		
	cout << "	Completed!" << endl << endl;
	TSP.SetPathSize(n_cities);
	for (int j=0; j<n_cities; j++) TSP.SetCity(j, j);
	TSP.Shuffle();
	TSP.SetLength();
	
	Temp = Temp_in;
	beta = 1./Temp;
	
}

void Metropolis() {

	r = rnd.Rannyu();
	New_TSP = TSP;
	
	New_TSP.Pair_Permutation();
	New_TSP.SetLength();
	if (New_TSP.GetLength() > TSP.GetLength()) prob = exp(-beta*(New_TSP.GetLength() - TSP.GetLength()));
	else prob = 1.;
	if (r<prob) {
		TSP = New_TSP;
		//N_accept[0]++;
	}
	//N_tot[0]++;
		
	New_TSP = TSP;
	New_TSP.Permutation();
	New_TSP.SetLength();
	if (New_TSP.GetLength() > TSP.GetLength()) prob = exp(-beta*(New_TSP.GetLength() - TSP.GetLength()));
	else prob = 1.;
	if (r<prob) {
		TSP = New_TSP;
		//N_accept[1]++;
	}
	//N_tot[1]++;
	
	New_TSP = TSP;
	New_TSP.Partial_Shift();
	New_TSP.SetLength();
	if (New_TSP.GetLength() > TSP.GetLength()) prob = exp(-beta*(New_TSP.GetLength() - TSP.GetLength()));
	else prob = 1.;
	if (r<prob) {
		TSP = New_TSP;
		//N_accept[2]++;
	}
	//N_tot[2]++;
		
	New_TSP = TSP;
	New_TSP.Inversion();
	New_TSP.SetLength();
	if (New_TSP.GetLength() > TSP.GetLength()) prob = exp(-beta*(New_TSP.GetLength() - TSP.GetLength()));
	else prob = 1.;
	if (r<prob) {
		TSP = New_TSP;
		//N_accept[3]++;
	}
	//N_tot[3]++;
	
}

void Print (string filename, string mode) {

	output.open(filename);
	
	if (mode == "circumference" ) {
		for (int i=0; i<n_cities; i++) output << x_circ[TSP.GetCity(i)] << "	" << y_circ[TSP.GetCity(i)] << endl;
		output << x_circ[TSP.GetCity(0)] << "	" << y_circ[TSP.GetCity(0)] << endl;
	}
	
	if (mode == "square" ) {
		for (int i=0; i<n_cities; i++) output << x_square[TSP.GetCity(i)] << "	" << y_square[TSP.GetCity(i)] << endl;
		output << x_square[TSP.GetCity(0)] << "	" << y_square[TSP.GetCity(0)] << endl;
	}
	

	output.close();

}


int Pbc (int n, int size) {

	return n - size*(int)((double)n/(double)size);
	
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
