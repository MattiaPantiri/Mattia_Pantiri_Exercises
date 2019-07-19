/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <vector>

using namespace std;

class RW { //class Random Walk

public:
	RW () {m_x = 0.; m_y = 0.; m_z = 0.;}
	
	void incrX () {m_x++;}
	void incrY () {m_y++;}
	void incrZ () {m_z++;}
	void decrX () {m_x--;}
	void decrY () {m_y--;}
	void decrZ () {m_z--;}
	double Square_Modulus() {return m_x*m_x + m_y*m_y + m_z*m_z;}
	void Spherical(double r, double theta, double phi) {m_x += r*sin(theta)*cos(phi); m_y += r*sin(theta)*sin(phi); m_z += r*cos(theta);}

private:
	double m_x, m_y, m_z;
};

//variables, parameters
int n_steps, n_simulations, n_blocks, block_length;
double r, theta, phi;
vector<double> prog_sum, prog_blkaverage, prog_squared_blkaverage;
vector<double> blkaverage, squared_blkaverage;
vector<double> average, squared_average;
vector<double> error;
int h;

//random numbers
#include "random.h"
int seed[4];
Random rnd;

//output stream
ofstream output;

//functions
void Input();
void Initialize();
void Block_Average(int i);
void Print_FinalRW();
void Simulate_DiscreteRW(RW* randw);
void Simulate_ContinuosRW(RW* randw);


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
