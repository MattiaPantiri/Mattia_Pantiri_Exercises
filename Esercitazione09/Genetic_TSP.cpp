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
#include "Genetic_TSP.h"

using namespace std;

int main () {
	
	Input();
	
	
	Initialize("circumference");
	if (Check() == false) {
		cout << endl << "You visited the same city twice!!! The program will shut down!" << endl << endl;
		return -1;
	}
	SortPopulation();
	
	output.open("GenerationBP_Circle.txt");
	cout << "	Simulating the evolution of the population..." << endl;
	for (int k=0; k<n_generations; k++) {
		GenerateNewPopulation();
		if (Check() == false) {
			cout << endl << "You visited the same city twice!!! The program will shut down!" << endl << endl;
			return -1;
		}	
		for (int i=0; i<n_elements; i++) Population[i].SetLength();
		SortPopulation();
		ave = Average();
		output << k << "	" << Population[0].GetLength() << "	" << ave << endl;
	}
	cout << "	Completed!" << endl << endl;
	output.close();
	
	PrintBestPath("BestPathCircle.txt");
	
	
	Initialize("square");
	if (Check() == false) {
		cout << endl << "You visited the same city twice!!! The program will shut down!" << endl << endl;
		return -1;
	}
	SortPopulation();
	
	output.open("GenerationBP_Square.txt");
	cout << "	Simulating the evolution of the population..." << endl;
	for (int k=0; k<n_generations; k++) { 
		GenerateNewPopulation();
		if (Check() == false) {
			cout << endl << "You visited the same city twice!!! The program will shut down!" << endl << endl;
			return -1;
		}	
		for (int i=0; i<n_elements; i++) Population[i].SetLength();
		SortPopulation();
		ave = Average();
		output << k << "	" << Population[0].GetLength() << "	" << ave << endl;
	}
	cout << "	Completed!" << endl << endl;
	output.close();
	
	PrintBestPath("BestPathSquare.txt");
	

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

void Individual::Shift() {

	int n = (int)rnd.Rannyu(1., n_cities);
	vector<int> appo = _path;
	for (int i=0; i<n_cities; i++) _path[Pbc(i+n, n_cities)] = appo[i];
	
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
	int j = (int)rnd.Rannyu(0., n_cities);
	int k = (int)rnd.Rannyu(0., n_cities);
	
	for (int i=j; i<j+n; i++) {
		swap(_path[Pbc(i, n_cities)], _path[Pbc(k, n_cities)]);
		k++;
	}
	
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
	
	cout << endl << "Travelling Salesman Problem (TSP) - Genetic Algorithm" << endl << endl;


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
	ReadInput.open("input.dat");
  
	ReadInput >> n_cities;
	cout << "--> Number of cities = " << n_cities << endl;
	
	ReadInput >> n_elements;  
	cout << "--> Number of population individuals = " << n_elements << endl;
	
	ReadInput >> n_generations;
	cout << "--> Number of generations = " << n_generations << endl << endl;
	
	ReadInput >> R;
	ReadInput >> L;

	ReadInput.close();
	
}

void Initialize(string mode) {
	
	x.resize(n_cities);
	y.resize(n_cities);
	
	if (mode == "circumference") {
	
		cout << "1) The cities are set randomly on a circumference of radius " << R << endl << endl;
	
		for (int i=0; i<n_cities; i++) {
			theta = rnd.Rannyu(0., 2*M_PI);
			x[i] = R*cos(theta);
			y[i] = R*sin(theta);
		}
	
	}
	
	if (mode == "square") {
	
		cout << "2) The cities are set randomly inside a square of side " << L << endl << endl;
	
		for (int i=0; i<n_cities; i++) {
			x[i] = rnd.Rannyu(-L/2.0, L/2.0);
			y[i] = rnd.Rannyu(-L/2.0, L/2.0);
		}
	
	}
	
	distances.resize(n_cities);
	for (int i=0; i<n_cities; i++) distances[i].resize(n_cities);
	
	cout << "	Computing distances between cities..." << endl;
	
	for (int i=0; i<n_cities; i++) {
		for (int j=i; j<n_cities; j++) {
			distances[i][j] = sqrt((x[i] - x[j])*(x[i] - x[j]) + (y[i] - y[j])*(y[i] - y[j]));
			distances[j][i] = distances[i][j];
		}
	}
	
	cout << "	Completed!" << endl << endl;	
	
	Population.resize(n_elements);
	
	cout << "	Loading starting population..." << endl;
	
	for (int i=0; i<n_elements; i++) {
		Population[i] = Individual(n_cities);
		for (int j=0; j<n_cities; j++) {
			Population[i].SetCity(j, j);
		}
		Population[i].Shuffle();
		Population[i].SetLength();
	}
	
	cout << "	Loaded!" << endl << endl;
	
}

bool Check() {

	for (int i=0; i<n_elements; i++) {
		for (int j=0; j<n_cities; j++) {
			for (int k=j+1; k<n_cities; k++) {
				if (Population[i].GetCity(j) == Population[i].GetCity(k)) return false;
			}
		}
	}
	
	return true;

}

void SortPopulation() {

	sort(Population.begin(), Population.end(), mySort);
	
}


Individual Select() {
	
	int j = (int) (n_elements*pow(rnd.Rannyu(), 2));
	return Population[j];

}

void Generate() {

	Father.SetPathSize(n_cities);
	Mother.SetPathSize(n_cities);
	Son.SetPathSize(n_cities);
	Daughter.SetPathSize(n_cities);
	
	Father = Select();
	Mother = Select();
	
	int counter;
	int j = (int)rnd.Rannyu(0., n_cities-1);
	
	for (int i=0; i<j; i++) {
		Son.SetCity(Father.GetCity(i), i);
		Daughter.SetCity(Mother.GetCity(i), i);
	}
		
	for (int i=j; i<n_cities; i++) {
		for (int k=0; k<n_cities; k++) {
			counter = 0;
			for (int h=0; h<i; h++) {
				if (Daughter.GetCity(h) != Father.GetCity(k)) counter++;
			}
			if (counter == i) {
				Daughter.SetCity(Father.GetCity(k), i);
				break;				
			}
		}
	}
	
	for (int i=j; i<n_cities; i++) {
		for (int k=0; k<n_cities; k++) {
			counter = 0;
			for (int h=0; h<i; h++) {
				if (Son.GetCity(h) != Mother.GetCity(k)) counter++;
			}
			if (counter == i) {
				Son.SetCity(Mother.GetCity(k), i);
				break;
			}
		}
	}
			
}

void GenerateNewPopulation() {
	
	double r1, r2;
	New_Population.resize(n_elements);		

	for (int i=0; i<n_elements; i+=2) {
		r1 = rnd.Rannyu();
		r2 = rnd.Rannyu();
		Generate();
		if(r1<0.8 and r1>=0.75) Son.Pair_Permutation();
		if(r1<0.85 and r1>=0.8) Son.Shift();
		if(r1<0.9 and r1>=0.85) Son.Partial_Shift();
		if(r1<0.95 and r1>=0.9) Son.Permutation();
		if(r1>=0.95) Son.Inversion();
		if(r2<0.8 and r2>=0.75) Daughter.Pair_Permutation();
		if(r2<0.85 and r2>=0.8) Daughter.Shift();
		if(r2<0.9 and r2>=0.85) Daughter.Partial_Shift();
		if(r2<0.95 and r2>=0.9) Daughter.Permutation();
		if(r2>=0.95) Daughter.Inversion();
		New_Population[i] = Son;
		New_Population[i+1] = Daughter;
	}
	
	Population = New_Population;

}
	
void PrintBestPath (string filename) {

	output.open(filename);
	
	for (int i=0; i<n_cities; i++) {
		output << x[Population[0].GetCity(i)] << "	" << y[Population[0].GetCity(i)] << endl;
	}
	
	output << x[Population[0].GetCity(0)] << "	" << y[Population[0].GetCity(0)] << endl;

	output.close();

}


double Average() {
	
	double sum = 0.;
	int half = n_elements/2;
	for (int i=0; i<half; i++) sum += Population[i].GetLength();
	
	return sum/(double)half;
	
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
