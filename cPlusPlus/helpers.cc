#include "helpers.h"

void target(ImageTemplate<double>* input, int x0, int y0, int B, float color){
	for (int y=y0-B; y<y0+B; y++){
		for (int x=x0-B; x<x0+B; x++){
			input->Pixel(x,y) = color;
		}
	}	
}

void scale(ImageTemplate<double>* input) {

	// min and max
	double m = input->Pixel(0,0);
	double M = input->Pixel(0,0);

	int x; int y;

	for (y=0; y<input->Height(); y++){
		for (x=0; x<input->Width(); x++){
			if (input->Pixel(x,y) < m)
				m = input->Pixel(x,y);
			if (input->Pixel(x,y) > M)
				M = input->Pixel(x,y);
		}
	}

	for (y=0; y<input->Height(); y++){
		for (x=0; x<input->Width(); x++){
			input->Pixel(x,y) = (255.0/(M-m)) * (input->Pixel(x,y) - m);
		}
	}


}

void crop(ImageTemplate<double>* input, ImageTemplate<double>* output, int W, int H, int xn, int yn) {

	output->Resize(W,H);

	for(int y=yn-H/2; y<yn+H/2; y++){
		for(int x=xn-W/2; x<xn+W/2; x++){

			output->Pixel(x-(xn-W/2),y-(yn-H/2)) = input->Pixel(x,y);
		}
	}

}

int gcd(int n, int m) {
	int gcd, remainder;
 
	while (n != 0) {
		remainder = m % n;
		m = n;
		n = remainder;
	}
	gcd = m;
	return gcd;
}

void interpolate(ImageTemplate<double>* input, ImageTemplate<double>* output, int U) {

	// input dimensions
	int X = input->Width();
	int Y = input->Height();

	// resize output
	output->Resize(U*X,U*Y);

	// two pixels used for interpolating between
	double a;
	double b;

	// interpolation fraction
	double f = 0.0;

	// iterating variables
	int x;
	int y;
	int z;

	// populate the sparse output image
	for (y=0; y<Y; y++){
		for (x=0; x<X; x++){
			output->Pixel(U*x, U*y) = input->Pixel(x,y);
		}
	}

	// interpolate rows by moving along horizontals
	for (y=0; y<U*Y-U; y+=U){

		for(x=0; x<U*X-U; x+=U){

			a = output->Pixel(x,y);
			b = output->Pixel(x+U,y);

			for (z=1; z<U; z++){
				f = (double)z / (double)U;
				output->Pixel(x+z,y) = (1.0-f)*a + f*b;
			}
		}
	}

	// interpolate columns by moving down verticals
	for (y=0; y<U*Y-U; y+=U) {

		for(x=0; x<U*X; x++){

			a = output->Pixel(x,y);
			b = output->Pixel(x,y+U);

			for (z=1; z<U; z++){
				f = (double)z / (double)U;
				output->Pixel(x,y+z) = (1.0-f)*a + f*b;
			}
		}
	}

	// skipping ZOH stuff cause cropping anyway

}

void downsampleRBJ(ImageTemplate<double>* input, ImageTemplate<double>* output, int D, float BW) {

	// input dimensions
	int X = input->Width();
	int Y = input->Height();	

	// iterating variables
	int x;
	int y;
	int z;

	// images needed for filtering steps
	ImageTemplate<double> filterRow;
	filterRow.Resize(X,Y);
	for(y=0; y<filterRow.Height(); y++){
		for(x=0; x<filterRow.Width(); x++){
			filterRow.Pixel(x,y) = 0.0;
		}
	}

	ImageTemplate<double> decimateCol;
	decimateCol.Resize(X/D,Y);
	ImageTemplate<double> filterCol;
	filterCol.Resize(X/D,Y);

	// define output size
	output->Resize(X/D,Y/D);

	double wc = M_PI / D;
	double cos_wc = cos(wc);
	double sin_wc = sin(wc);

	double alpha = sin_wc * sinh( (log(2.0)/2.0) * BW * (wc/sin_wc) );

	// can optimize since a[0] same as a[2]
	double b[3] = {(1.0 - cos_wc)/2.0, 1.0 - cos_wc, (1 - cos_wc)/2.0};
	double a[3] = {1.0 + alpha, -2.0*cos_wc, 1.0 - alpha};
	double c[5] = {b[0]/a[0],b[1]/a[0],b[2]/a[0],a[1]/a[0],a[2]/a[0]};

	/*
	printf("alpha = %f\n",alpha);

	for (int i=0; i<3; i++){
		printf("b[%d] = %f, a[%d] = %f\n", i,b[i],i,a[i]);
	}

	for (int i=0; i<5; i++){
		printf("c[%d] = %f\n",i,c[i]);
	}
	*/
	
	// filter each row
	for (y=0; y<Y; y++){
		for (x=2; x<X; x++){

			filterRow.Pixel(x,y) = c[0]*input->Pixel(x,y) +
								   c[1]*input->Pixel(x-1,y) +
								   c[2]*input->Pixel(x-2,y) -
								   c[3]*filterRow.Pixel(x-1,y) -
								   c[4]*filterRow.Pixel(x-2,y);
		}
	}
		
	// filterRow.SavePng("1filterRow.png");
	
	// decimate the rows that were just filtered
	for (y=0; y<Y; y++){
		for (x=0; x<Y/D; x++){
			decimateCol.Pixel(x,y) = filterRow.Pixel(x*D,y);
		}
	}

	// decimateCol.SavePng("2decimateCol.png");

	// filter each column
	for (x=0; x<X/D; x++) {
		for (y=2; y<Y; y++){
			filterCol.Pixel(x,y) = c[0]*decimateCol.Pixel(x,y) +
								   c[1]*decimateCol.Pixel(x,y-1) +
								   c[2]*decimateCol.Pixel(x,y-2) -
								   c[3]*filterCol.Pixel(x,y-1) -
								   c[4]*filterCol.Pixel(x,y-1);
		}
	}

	// filterCol.SavePng("3filterCol.png");
	
	// decimate the columns that were just filtered
	for (y=0; y<Y/D; y++){
		for (x=0; x<X/D; x++){
			output->Pixel(x,y) = filterCol.Pixel(x,y*D);
		}
	}
	
	// output->SavePng("4output.png");
	
}