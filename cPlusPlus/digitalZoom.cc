#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include "imagefft.h"
#include "helpers.h"

using namespace std;

// argv inputs:
//
// 1         2              3                4            5                      6
// input.png upsampleFactor downsampleFactor numberFrames destinationCoordinateX destinationCoordinateY

int main (int argc, char* argv[]){

	// load the input image
	// ImageTemplate<double> input;
	// input.LoadPng (argv[1]);

	// interpolation factor
	int U = atoi(argv[2]);
	int D = atoi(argv[3]);

	double BW = 2.5;

	// make images for each frame of movie
	int N = atoi(argv[4]);
	ImageTemplate<double> movie[N];	

	double alpha = (double)U / (double)D;
	// double scale = pow(alpha, 1.0/(double)N);
	// printf("alpha: %f\tscale: %f\n", alpha, scale);
	
	ImageTemplate<double> interpolated;
	ImageTemplate<double> decimated;

	movie[0].LoadPng(argv[1]);
	int W = movie[0].Width();
	int H = movie[0].Height();

	

	// origin and destination coordinates
	double xi = (double)W/2.0;
	double yi = (double)H/2.0;
	int alphaXi = (int)round(alpha*xi);
	int alphaYi = (int)round(alpha*yi);

	double xf = (double)atoi(argv[5]);
	double yf = (double)atoi(argv[6]);

	double xDelta = (xf-xi)/(double)(N-1);
	double yDelta = (yf-yi)/(double)(N-1);
	printf("xDelta: %f\tyDelta: %f\n", xDelta, yDelta);

	// mark the center and the destination as a target
	int B = 5;
	target(&movie[0], (int)round(xf), (int)round(yf), B, 200.0);
	target(&movie[0], (int)round(xi), (int)round(yi), B, 0.0);

	// mark the deltas on original
	for (int n=1; n<N; n++){
		target(&movie[0], (int)round(xi+n*xDelta), (int)round(yi+n*yDelta), 3, 0.0);
	}

	movie[0].SavePng("0.png");

	cout << endl;
	
	// intermediate coordinates updated each itealphan
	double xn = xDelta;
	double yn = yDelta;
	int xnInt;
	int ynInt;

	for(int n=1; n<N; n++){

		interpolate(&movie[n-1], &interpolated, U);
		// printf("\n\n  upsize: %d x %d\n", interpolated.Width(), interpolated.Height());
		
		downsampleRBJ(&interpolated, &decimated, D, BW);
		// printf("updnsize: %d x %d\n", decimated.Width(), decimated.Height());

		// need position of next center in terms of the blown-up image
		xn = alphaXi + pow(alpha, n)*xDelta;
		xnInt = (int)round(xn);

		yn = alphaYi + pow(alpha, n)*yDelta;
		ynInt = (int)round(yn);

		target(&decimated, xnInt, ynInt, 5, 255.0);

		crop(&decimated, &movie[n], W, H, xnInt, ynInt);

		cout << "saving movie[" << n << "]" << endl;
		movie[n].SavePng(std::to_string(n)+".png");
		
	}

	std::string command = "python /Users/justinsconza/Documents/ECE418/project418/cPlusPlus/generateAnimation.py " + std::to_string(N);
	cout << endl << "calling: " << command << endl << endl;;	
	system(command.c_str());

	return 0;
}
















/*

	// interpolated image
	ImageTemplate<double> interpolated;
	ImageTemplate<double> decimated;
	
	interpolate(&input, &interpolated, U);
	downsampleRBJ(&interpolated, &decimated, D, BW);

*/





