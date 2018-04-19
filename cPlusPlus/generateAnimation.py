import sys

import numpy as np
from numpy import sin as sin
from numpy import cos as cos
from numpy import pi as PI
import scipy as sp

import time

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as img

import imageio

images = []
for i in range(int(sys.argv[1])):
    file_str = str(i) + '.png'
    images.append(imageio.imread(file_str))

output_file = 'outputAnimation.gif'
imageio.mimsave(output_file, images)

