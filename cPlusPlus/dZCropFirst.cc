#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include "imagefft.h"
#include "helpers.h"

using namespace std;

#define TARGET_FLAG 0

// argv inputs:
//
// 1         | 2                    | 3             | 4                        | 5
// input.png | scale factor (float) | number frames | destination coordinate x | destination coordinate y

int main (int argc, char* argv[]){

	// scale factor (float)
	double M = atof(argv[2]);

	// number of frames in output movie
	int N = atoi(argv[3]);

	// array of images for movie
	ImageTemplate<double> movie[N];	
	movie[0].LoadPng(argv[1]);
	int W = movie[0].Width();
	int H = movie[0].Height();

	// resolution of intermediate up and down factors
	double dUp = 10.0;

	// intermediate up and down factors
	int u = (int) dUp*round(dUp*pow(M, 1.0/(double)N))/dUp;
	int d = dUp;

	// reduce u and d factors to be relatively prime
	int theGcd = gcd(u,d);
	u /= theGcd;
	d /= theGcd;

	if(u==1){
		cout << "\n\nToo many frames. Choose a lower value for N\n\n";
		return 1;
	}	

	// image scale factor each iteration
	double alpha = (double)u / (double)d;
	int wOverAlpha = (int)round((double)W/alpha);
	int hOverAlpha = (int)round((double)H/alpha);

	// bandwidth for IIR anti-aliasing filter, lower --> sharper final resolution
	double BW;
	if(TARGET_FLAG)
		BW = 2.5;
	else
		BW = 1.5;

	// origin coordinates are the middle of the image
	double xi = (double)W/2.0;
	double yi = (double)H/2.0;

	// final coordinates
	double xf = (double)atoi(argv[4]);
	double yf = (double)atoi(argv[5]);

	if (xf < 0 || xf > W || yf < 0 || yf > H){
		cout << "\n\nFinal coordinates must be within range of original image.\n\n";
		return 1;
	}

	// step size relative to original image
	double xDelta = (xf-xi)/(double)(N-1);
	double yDelta = (yf-yi)/(double)(N-1);

	cout << endl;

	cout << "u: " << u << endl;
	cout << "d: " << d << "\n\n";

	cout << "actual alpha : " << pow(M, 1.0/(double)N) << endl;
	cout << "rounded alpha: " << alpha << "\n\n";	

	if(TARGET_FLAG) {

		// mark the center, each step and final destination as targets
		target(&movie[0], (int)round(xf), (int)round(yf), 5, 200.0);
		target(&movie[0], (int)round(xi), (int)round(yi), 5, 0.0);
		for (int n=1; n<N; n++){
			target(&movie[0], (int)round(xi+n*xDelta), (int)round(yi+n*yDelta), 3, 0.0);
		}

	}
	
	// intermediate coordinates updated each iteration
	double xn;
	double yn;
	int xnInt;
	int ynInt;

	ImageTemplate<double> cropped;
	ImageTemplate<double> interpolated;

	movie[0].SavePng("0.png");

	for(int n=1; n<N; n++){

		// get crop coordinates
		xn = xi + pow(alpha, n-1)*xDelta;		
		xnInt = (int)round(xn);

		yn = yi + pow(alpha, n-1)*yDelta;
		ynInt = (int)round(yn);

		// do the crop
		crop(&movie[n-1], &cropped, wOverAlpha, hOverAlpha, xnInt, ynInt);
		
		// upsample then downsample
		interpolate(&cropped, &interpolated, u);		
		downsampleRBJ(&interpolated, &movie[n], d, W, H, BW);

		if (movie[n].Height() != H)
			movie[n].Resize(W,H);

		// copy paste cols and rows along all four edges out to the edges
		zeroOrderHold(&movie[n], 5);

		if (TARGET_FLAG)
			target(&movie[n], xi, yi, 5, 255.0);

		// save the frame, start over
		cout << "saving movie[" << n << "]" << endl;
		movie[n].SavePng(std::to_string(n)+".png");
		
	}

	// call python script to generate animated gif from frames
	std::string command = "python /Users/justinsconza/Documents/ECE418/project418/cPlusPlus/generateAnimation.py " + std::to_string(N);
	cout << endl << "calling: " << command << endl << endl;;	
	system(command.c_str());
	

	return 0;
}







