#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include "imagefft.h"
#include "helpers.h"

using namespace std;



void crop(ImageTemplate<double>* input, ImageTemplate<double>* output, int W, int H, int xn, int yn) {

	output->Resize(W,H);

	target(input, xn, yn, 9, 255.0);

	for(int y=yn-H/2; y<yn+H/2; y++){

		for(int x=xn-W/2; x<xn+W/2; x++){

			output->Pixel(x-(xn-W/2),y-(yn-H/2)) = input->Pixel(x,y);
		}
	}

}

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

	double ratio = (double)U / (double)D;
	double scale = pow(ratio, 1.0/(double)N);
	// printf("ratio: %f\tscale: %f\n", ratio, scale);
	
	ImageTemplate<double> interpolated;
	ImageTemplate<double> decimated;

	movie[0].LoadPng(argv[1]);
	

	// origin and destination coordinates
	double xI = (double)movie[0].Width()/2.0;
	double yI = (double)movie[0].Height()/2.0;
	double xF = (double)atoi(argv[5]);
	double yF = (double)atoi(argv[6]);

	// mark the center and the destination as a target
	int B = 9;
	target(&movie[0], xF, yF, B, 255.0);
	target(&movie[0], movie[0].Width()/2, movie[0].Height()/2, B, 0.0);

	movie[0].SavePng("0.png");

	
	// intermediate coordinates updated each iteration
	double xn;
	double yn;

	for(int n=1; n<N; n++){

		interpolate(&movie[n-1], &interpolated, U);
		// printf("\n\n  upsize: %d x %d\n", interpolated.Width(), interpolated.Height());
		
		downsampleRBJ(&interpolated, &decimated, D, BW);
		// printf("updnsize: %d x %d\n", decimated.Width(), decimated.Height());

		xn = ratio * (xI + (1.0/(double)N) * ratio * (xF - xI)); 
 		yn = ratio * (yI + (1.0/(double)N) * ratio * (yF - yI));

		crop(&decimated, &movie[n], movie[0].Width(), movie[0].Height(), (int)round(xn), (int)round(yn));

		cout << "saving movie[" << n << "]" << endl;
		movie[n].SavePng(std::to_string(n)+".png");
		
	}

	return 0;
}
















/*

	// interpolated image
	ImageTemplate<double> interpolated;
	ImageTemplate<double> decimated;
	
	interpolate(&input, &interpolated, U);
	downsampleRBJ(&interpolated, &decimated, D, BW);

*/





