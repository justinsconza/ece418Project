	// input center to zoom-in on coordinates (these should be passed in)
	// int xC = input->Width()/2;
	// int yC = input->Height()/2;

	// upper corner of input cropped coordinates
	// int x0 = xC - output->Width()/2;
	// int y0 = yC - output->Height()/2;

	// for(int y=y0; y<y0+H; y++){
	// 	for(int x=x0; x<x0+W; x++){
	// 		output->Pixel(x-x0,y-y0) = input->Pixel(x,y);
	// 	}
	// }

	// printf("xn:%d\tyn:%d\n",xn,yn);



// radius and angle to destination
	double R = sqrt( pow(yF-yI,2.0) + pow(xF-xI,2.0) );
	double theta;
	if ((int)xF == (int)xI) 
		theta = 0.0;
	else
		theta = abs( atan( (yF-yI) / (xF-xI) ) );

	// printf("R=%f\ttheta=%f\n\n",R,theta);

	// destination quadrant
	int quadrant;
	if ((int)xF > (int)xI & (int)yF > (int)yI)
		quadrant = 1;
	if ((int)xF < (int)xI & (int)yF > (int)yI) {
		quadrant = 2;
		theta = M_PI - theta;
	}
	if (xF < xI & yF < yI) {
		quadrant = 3;
		theta = theta - M_PI;
	}
	if (xF > xI & yF < yI) {
		quadrant = 4;
		theta *= -1.0;
	}