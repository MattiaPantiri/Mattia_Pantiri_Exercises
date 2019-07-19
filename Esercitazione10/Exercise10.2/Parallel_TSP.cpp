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
#include "Parallel_TSP.h"
#include "mpi.h"

using namespace std;

int main (int argc, char** argv) {

	MPI::Init(argc, argv);
	Size = MPI::COMM_WORLD.Get_size();
	Rank = MPI::COMM_WORLD.Get_rank();
	
	Input(Rank);
	
	Initialize("circumference", Rank);
	//int N = 1;
	if (Rank == 0) cout << "	Simulating..." << endl;
	/*string name = "Iteration_BP_Circle";
	string extension = ".txt";
	output.open(name + to_string(Rank) + extension);*/
	for (int k=1; k<=n_temp; k++) {
		/*for (int i=0; i<4; i++) {
			N_accept[i] = 0;
			N_tot[i] = 0;
		}*/
		for (int j=0; j<k*100; j++) {
			Metropolis();
			//output << N << "	" << TSP.GetLength() << endl;
			//N++;
		}
		lengths[k-1] = TSP.GetLength();
		//for (int i=0; i<4; i++) acc[i] = (double)N_accept[i]/(double)N_tot[i];
		Temp=Temp*Scale_Factor;
		beta = 1./Temp;
	}
	//output.close();
	if (Rank == 0) cout << "	Completed!" << endl << endl;
	//name = "BestPathCircle";
	//Print(name + to_string(Rank) + extension, "circumference");
	
	double irecv[Size];
	for (int i=0; i<Size; i++) irecv[i] = 0.0;
	double x = TSP.GetLength();
	MPI_Gather(&x, 1, MPI_DOUBLE_PRECISION, irecv, 1, MPI_DOUBLE_PRECISION, 0, MPI::COMM_WORLD);
	best_rank = pos_min(irecv, Size);
	if (Rank == 0) {
		cout << endl << "	The 4 cores have found this best path lengths:" << endl;
		cout << "	--> ";
		for (int i=0; i<Size; i++){
			cout << irecv[i] << "	";
		}
		cout << endl;
		cout << "	--> Minimum: " << min_value(irecv, 4) << endl;
		cout << "	--> By core number: " << best_rank << endl << endl;
	}
	if (Rank == best_rank) Print("circumference");
	
	
	
	Initialize("square", Rank);
	//N = 1;
	if (Rank == 0) cout << "	Simulating..." << endl;
	/*name = "Iteration_BP_Square";
	extension = ".txt";
	output.open(name + to_string(Rank) + extension);*/
	for (int k=1; k<=n_temp; k++) {
		/*for (int i=0; i<4; i++) {
			N_accept[i] = 0;
			N_tot[i] = 0;
		}*/
		for (int j=0; j<k*100; j++) {
			Metropolis();
			//output << N << "	" << TSP.GetLength() << endl;
			//N++;
		}
		lengths[k-1] = TSP.GetLength();
		//for (int i=0; i<4; i++) acc[i] = (double)N_accept[i]/(double)N_tot[i];
		Temp=Temp*Scale_Factor;
		beta = 1./Temp;
	}
	//output.close();
	if (Rank == 0) cout << "	Completed!" << endl << endl;
	//name = "BestPathSquare";
	//Print(name + to_string(Rank) + extension, "square");

	for (int i=0; i<Size; i++) irecv[i] = 0.0;
	x = TSP.GetLength();
	MPI_Gather(&x, 1, MPI_DOUBLE_PRECISION, irecv, 1, MPI_DOUBLE_PRECISION, 0, MPI::COMM_WORLD);
	
	if (Rank == 0) {
		cout << endl << "	The 4 cores have found this best path lengths:" << endl;
		cout << "	--> ";
		for (int i=0; i<Size; i++){
			cout << irecv[i] << "	";
		}
		cout << endl;
		cout << "	--> Minimum: " << min_value(irecv, 4) << endl;
		cout << "	--> By core number: " << best_rank << endl << endl;
	}
	if (Rank == best_rank) Print("square");
	
	
	MPI::Finalize();

return 0;
}


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

	int i = (int)rnd[Rank].Rannyu(0., n_cities);
	int j = (int)rnd[Rank].Rannyu(0., n_cities);
	while (j==i) j = (int)rnd[Rank].Rannyu(0., n_cities);
	
	swap(_path[i], _path[j]);	

}

void Individual::Partial_Shift() {
	
	int n = (int)rnd[Rank].Rannyu(1., n_cities/2);
	int m = (int)rnd[Rank].Rannyu(1., n_cities/2);
	int j = (int)rnd[Rank].Rannyu(0., n_cities);
	
	vector<int> appo(n+m);

	for (int i=0; i<n+m; i++) appo[Pbc(i+n, n+m)] = _path[Pbc(j+i, n_cities)];

	for (int i=0; i<n+m; i++) _path[Pbc(j+i, n_cities)] = appo[i];

}

void Individual::Permutation() {

	int n = (int)rnd[Rank].Rannyu(0., n_cities/2);
	int j = (int)rnd[Rank].Rannyu(0., n_cities/2);
	
	for (int i=j; i<j+n; i++) swap(_path[i], _path[Pbc(i+n, n_cities)]);
	
}

void Individual::Inversion() {

	int n = (int)rnd[Rank].Rannyu(0., n_cities);
	int j = (int)rnd[Rank].Rannyu(0., n_cities);
	
	vector<int> appo(n);
	
	for (int i=0; i<n; i++) appo[i] = _path[Pbc(j+n-i-1, n_cities)];
	
	for (int i=j; i<j+n; i++) _path[Pbc(i, n_cities)] = appo[i-j];
	
}


//other funcions

void Input(int a) {

	ifstream ReadInput;

	//read seed for city setting
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	s.SetRandom(seed,p1,p2);
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
	if (a == 0) {
		cout << endl << "Travelling Salesman Problem (TSP) - Parallel Simulated Annealing" << endl << endl;
		cout << "--> Number of cities = " << n_cities << endl;
		cout << "--> Initial temperature = " << Temp_in << endl;
		cout << "--> The program visits " << n_temp << " temperatures" << endl;
		cout << "--> Temperature scales with a factor = " << Scale_Factor << endl << endl;
	}
	
	lengths.resize(n_temp);
	
	//set cities on a circumference and on a square (N.B. done with the same seed and Primes in order to give the same TSP problem to all cores)
	x_circ.resize(n_cities);
	y_circ.resize(n_cities);
	x_square.resize(n_cities);
	y_square.resize(n_cities);
	
	for (int i=0; i<n_cities; i++) {
		theta = s.Rannyu(0., 2*M_PI);
		x_circ[i] = R*cos(theta);
		y_circ[i] = R*sin(theta);
	}
	
	for (int i=0; i<n_cities; i++) {
		x_square[i] = s.Rannyu(-L, L);
		y_square[i] = s.Rannyu(-L, L);
	}
	
	
//Read seed for random numbers in the cores simulations
	int p[8];
	Primes.open("Primes");
	for (int i=0; i<8; i++) Primes >> p[i];
	Primes.close();

	for (int i=0; i<4; i++) rnd[i].SetRandom(seed,p[2*i],p[2*i+1]);
	
}

void Initialize(string mode, int a) {
	
	distances.resize(n_cities);
	for (int i=0; i<n_cities; i++) distances[i].resize(n_cities);
	
	if (a == 0) {
		if (mode == "circumference") cout << "1) The cities are set randomly on a circumference of radius " << R << endl << endl;
		if (mode == "square") cout << "2) The cities are set randomly inside a square of side " << 2*L << endl << endl;
		cout << "	Computing distances between cities..." << endl;
	}
	
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
		
	if (a == 0) cout << "	Completed!" << endl;
	TSP.SetPathSize(n_cities);
	for (int j=0; j<n_cities; j++) TSP.SetCity(j, j);
	TSP.Shuffle();
	TSP.SetLength();
	
	Temp = Temp_in;
	beta = 1./Temp;
	
}

void Metropolis() {

	r = rnd[Rank].Rannyu();
	New_TSP = TSP;
	
	New_TSP.Pair_Permutation();
	New_TSP.SetLength();
	if (New_TSP.GetLength() > TSP.GetLength()) prob = exp(-beta*(New_TSP.GetLength() - TSP.GetLength()));
	else prob = 1.;
	if (r<prob) {
		TSP = New_TSP;
		N_accept[0]++;
	}
	N_tot[0]++;
		
	New_TSP = TSP;
	New_TSP.Permutation();
	New_TSP.SetLength();
	if (New_TSP.GetLength() > TSP.GetLength()) prob = exp(-beta*(New_TSP.GetLength() - TSP.GetLength()));
	else prob = 1.;
	if (r<prob) {
		TSP = New_TSP;
		N_accept[1]++;
	}
	N_tot[1]++;
	
	New_TSP = TSP;
	New_TSP.Partial_Shift();
	New_TSP.SetLength();
	if (New_TSP.GetLength() > TSP.GetLength()) prob = exp(-beta*(New_TSP.GetLength() - TSP.GetLength()));
	else prob = 1.;
	if (r<prob) {
		TSP = New_TSP;
		N_accept[2]++;
	}
	N_tot[2]++;
		
	New_TSP = TSP;
	New_TSP.Inversion();
	New_TSP.SetLength();
	if (New_TSP.GetLength() > TSP.GetLength()) prob = exp(-beta*(New_TSP.GetLength() - TSP.GetLength()));
	else prob = 1.;
	if (r<prob) {
		TSP = New_TSP;
		N_accept[3]++;
	}
	N_tot[3]++;
	
}

void Print (string mode) {
	
	if (mode == "circumference" ) {
		output.open("Best_Path_Circle.txt");
		for (int i=0; i<n_cities; i++) output << x_circ[TSP.GetCity(i)] << "	" << y_circ[TSP.GetCity(i)] << endl;
		output << x_circ[TSP.GetCity(0)] << "	" << y_circ[TSP.GetCity(0)] << endl;
		output.close();
		output.open("Iteration_BP_Circle.txt");
		for (int i=0; i<n_temp; i++) output << i << "	" << lengths[i] << endl;
		output.close();
	}
	
	if (mode == "square" ) {
		output.open("Best_Path_Square.txt");
		for (int i=0; i<n_cities; i++) output << x_square[TSP.GetCity(i)] << "	" << y_square[TSP.GetCity(i)] << endl;
		output << x_square[TSP.GetCity(0)] << "	" << y_square[TSP.GetCity(0)] << endl;
		output.close();
		output.open("Iteration_BP_Square.txt");
		for (int i=0; i<n_temp; i++) output << i << "	" << lengths[i] << endl;
		output.close();
	}

}


int Pbc (int n, int size) {

	return n - size*(int)((double)n/(double)size);
	
}

double min_value (double* v, int size) {
	
	double minimum = v[0];
	for (int i=0; i<size; i++) {
		if (v[i] < minimum) minimum = v[i];
	}
	return minimum;
}

int pos_min (double*v, int size) {

	int pos = 0;
	double minimum = v[0];
	for (int i=0; i<size; i++) {
		if (v[i] < minimum) {
			pos = i;
			minimum = v[i];
		}
	}
	
	return pos;
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
