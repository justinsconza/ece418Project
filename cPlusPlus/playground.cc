#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include "imagefft.h"
#include "helpers.h"

using namespace std;



int main (int argc, char* argv[]){

	ImageTemplate<double> input;	
	input.LoadPng(argv[1]);

	ImageTemplate<double> interpolated;
	ImageTemplate<double> decimated;

	double BW = 2.5;

	double M = atof(argv[2]);

	int N = atoi(argv[3]);

	double dUp = 10.0;

	int u = (int) dUp*round(dUp*pow(M, 1.0/(double)N))/dUp;
	int d = dUp;

	cout << "u: " << u << endl;
	cout << "d: " << d << "\n\n";

	int theGcd = gcd(u,d);
	u /= theGcd;
	d /= theGcd;

	cout << "u: " << u << endl;
	cout << "d: " << d << "\n\n";

	cout << "(" << M << ")**(1/" << N << "): " << pow(M, 1.0/(double)N) << endl;
	cout << "u/d: " << (double)u/(double)d << "\n\n";

	cout << "interpolating..." << endl;
	interpolate(&input, &interpolated, u);	

	cout << "decimating..." << endl;
	downsampleRBJ(&interpolated, &decimated, d, BW);

	decimated.SavePng("playgroundOutput.png");
	
	return 0;
}







