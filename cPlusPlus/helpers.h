#ifndef HELPERS_H
#define HELPERS_H

#include <stdlib.h>
#include <iostream>
#include <vector>
#include "imagefft.h"
#include "helpers.h"

void target(ImageTemplate<double>* input, int x0, int y0, int B, float color);
void scale(ImageTemplate<double>* input);
void crop(ImageTemplate<double>* input, ImageTemplate<double>* output, int W, int H, int xn, int yn);
int gcd(int n, int m); 
void interpolate(ImageTemplate<double>* input, ImageTemplate<double>* output, int U);
void downsampleRBJ(ImageTemplate<double>* input, ImageTemplate<double>* output, int D, float BW);

#endif 