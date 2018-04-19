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

from helpers418 import *

#==============================================================================
# ok here's the deal: decimate_rbj is WAY better than decimate_FD BUT
# you pay for it.  decimate_rbj takes significantly longer to compute
#==============================================================================

#==============================================================================
# matplotlib.pyplot.close("all")
#==============================================================================

image_path = '/Users/justinsconza/Documents/ECE418/project418/python/lena.png'
image = img.imread(image_path,format=int)
I = image.shape[0]
J = image.shape[1]
   

U = 5                       # upsample factor
D = 2                       # downsample factor
N = 3                       # number of frames in movie

x = round((U/D)**(1/N),1)   # scale ratio each iteration, rounded to nearest tenth

u = int(x*10)               # upsample factor each iteration
d = 10                      # downsample factor each iteration

movie = np.zeros([I,J,N])   # movie of image arrays
movie[:,:,0] = image        # first frame is original

temp = image                # first image is original
for n in range(1,N):        # loop through remaining frames
    
    temp = decimate_FD(interpolate(upsample(temp,u),u),d)
    i = temp.shape[0]//2
    j = temp.shape[1]//2
    movie[:,:,n] = temp[i-I//2:i+I//2,j-J//2:j+J//2]


images = []
for i in range(N):
    file_str = 'output' + str(i) + '.png'
    sp.misc.imsave(file_str,movie[:,:,i])
    images.append(imageio.imread(file_str))

output_file = 'output_animation.gif'
imageio.mimsave(output_file, images)




#==============================================================================
# M_up = 47
# M_dn = 10
# 
# up = interpolate(upsample(image,M_up),M_up)
# decimated_FD = decimate_FD(up,M_dn)
# decimated_TD = decimate_rbj(up,M_dn ,2.5) 
#==============================================================================

#==============================================================================
# Comparison of times
#==============================================================================
#==============================================================================
# up = interpolate(upsample(image,5),5)
# decimated_FD = decimate_FD(up,4)
# decimated_TD = decimate_rbj(up,4 ,2.5) 
#==============================================================================

#==============================================================================
# fig,ax = plt.subplots(2,1, figsize=(14,7))
# ax[0].imshow(filtered1, interpolation='none', cmap='gray')
#==============================================================================

#==============================================================================
# sp.misc.imsave('output_FD.png',decimated_FD)
#==============================================================================
#==============================================================================
# sp.misc.imsave('output_TD.png',decimated_TD)
#==============================================================================
