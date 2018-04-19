import numpy as np
from numpy import sin as sin
from numpy import cos as cos
from numpy import pi as PI
import scipy as sp

import time

#==============================================================================
# UPSAMPLE
#==============================================================================

def upsample(image_in,N):
    
    t0 = time.time()
    
    I = image_in.shape[0]
    J = image_in.shape[1]
    image_out = np.zeros([N*I, N*J])
    for i in range(I):
        for j in range(J):
            image_out[N*i, N*j] = image_in[i,j]
    
    t1 = time.time()            
    print('upsample time: ',t1-t0)
    
    return image_out

#==============================================================================
# INTERPOLATE
#==============================================================================

def interpolate(image_in,N):
    
    t0 = time.time()
    
    I = image_in.shape[0]
    J = image_in.shape[1]
    
    image_out = np.zeros([I,J])
    
    # need to figure out scope with variable in functions...
    for i in range(I):
        for j in range(J):
            image_out[i,j] = image_in[i,j]
    
    # interpolate the rows by moving along the horizontals
    # move down a column in steps of N
    for i in range(0,I,N):
        
        # move along a row in steps of N and avoid overshooting with b
        for j in range(0,J-N,N):
            a = image_in[i,j]
            b = image_in[i,j+N]
            
            # do the row-wise interpolation
            for k in range(1,N):
                f = k/N
                image_out[i,j+k] = (1-f)*a + f*b
                
    # interpolate the columns by moving down the verticals
    for i in range(0,I-N,N):
        
        # move over along the row, one column at a time
        for j in range(J):
            a = image_out[i,j]
            b = image_out[i+N,j]
            
            # do the column-wise interpolation
            for k in range(1,N):
                f = k/N
                image_out[i+k,j] = (1-f)*a + f*b

            
    t1 = time.time()            
    print('interpolate time: ',t1-t0)
    
    return image_out
     

#==============================================================================
# DOWNSAMPLE W/O ANTI-ALIASING
#==============================================================================

def downsample(image_in,N):
    
    t0 = time.time()
    I = image_in.shape[0]
    J = image_in.shape[1]
    image_out = np.zeros([I//N,J//N])
    for i in range(I//N):
        for j in range(J//N):
            image_out[i,j] = image_in[i*N,j*N]
            
    t1 = time.time()            
    print('downsample time: ',t1-t0)
    return image_out


#==============================================================================
# DOWNSAMPLE W/ ANTI-ALIASING, SEPARABLE TIME-DOMAIN FILTER
#==============================================================================

# decimate one dimension at a time
# filter rows --> downsample rows --> filter cols --> downsample cols
def decimate_rbj(image_in,N,band_width):
    
    t0 = time.time()
    
    I = image_in.shape[0]
    J = image_in.shape[1]

    image_filter_row = np.zeros([I,J])
    image_decimate_col = np.zeros([I,J//N])
    image_filter_col = np.zeros([I,J//N])
    image_out = np.zeros([I//N,J//N])
            
    wc = PI/N
    cos_wc = cos(wc)
    sin_wc = sin(wc)
    
    # resonance, BW version
    # BW at 0 makes super peaky cutoff, BW at 1 makes very broad rolloff
    # use a BW of 1
    BW = band_width
    alpha = sin_wc*np.sinh( (np.log(2)/2) * BW * (wc/sin_wc) )

    # coefficients from RBJ cookbook
    b = np.array([(1 - cos_wc)/2, 1 - cos_wc, (1 - cos_wc)/2])
    a = np.array([1 + alpha, -2*cos_wc, 1 - alpha])
    
    
#==============================================================================
#     d = np.zeros(I)
#     d[0] = 1
#     h = np.zeros(I)
#==============================================================================
    
    # move down each row filter it
    for i in range(I):
        # filter each row
        # starting from second output b/c need two past samples for each input
        for j in range(2,J):
            image_filter_row[i,j] = (b[0]/a[0])*image_in[i,j] + \
                                    (b[1]/a[0])*image_in[i,j-1] + \
                                    (b[2]/a[0])*image_in[i,j-2] - \
                                    (a[1]/a[0])*image_filter_row[i,j-1] - \
                                    (a[2]/a[0])*image_filter_row[i,j-2]
                                    
#==============================================================================
#             h[j] = (b[0]/a[0])*d[j] + \
#                    (b[1]/a[0])*d[j-1] + \
#                    (b[2]/a[0])*d[j-2] - \
#                    (a[1]/a[0])*h[j-1] - \
#                    (a[2]/a[0])*h[j-2]
#                                     
#     H = np.abs(np.fft.fft(h))[0:I//2]
#==============================================================================
                   
    # decimate the rows that were just filtered by moving along columns
    for i in range(I):
        for j in range(J//N):
            image_decimate_col[i,j] = image_filter_row[i,j*N]
    
    # filter along the columns
    # move across each column
    for j in range(J//N):
        # filter each column
        # starting from second output b/c need two past samples for each input
        for i in range(2,J):
            image_filter_col[i,j] = (b[0]/a[0])*image_decimate_col[i,j] + \
                                    (b[1]/a[0])*image_decimate_col[i-1,j] + \
                                    (b[2]/a[0])*image_decimate_col[i-2,j] - \
                                    (a[1]/a[0])*image_filter_col[i-1,j] - \
                                    (a[2]/a[0])*image_filter_col[i-2,j]
    
    # decimate the columns that were just filtered by moving along rows
    for i in range(I//N):
        for j in range(J//N):
            image_out[i,j] = image_filter_col[i*N,j]
        
    t1 = time.time()            
    print('decimate_rbj time: ',t1-t0)    
    return image_out#, H

#==============================================================================
# DECIMATE IN FREQUENCY DOMAIN W/ ANTI-ALIASING ALL IN ONE STEP
#==============================================================================

def decimate_FD(image_in,N):
    
    
    t0 = time.time()
    
    I = image_in.shape[0]
    J = image_in.shape[1]
    
    image_in_fft = np.fft.fft2(image_in)

    small_fft = np.zeros([I//N,J//N],dtype=complex)
    
    for i in range(I//N//2):
        for j in range(J//N//2):
            small_fft[i,j] = image_in_fft[i,j]
            small_fft[i,J//N-1-j] = image_in_fft[i,J-1-j]
            small_fft[I//N-1-i,J//N-1-j] = image_in_fft[I//N-1-i,J-1-j]
            small_fft[I//N-1-i,j] = image_in_fft[I-1-i,j]
   
    #plt.imshow(np.log(abs(small_fft)))

    small_image = np.fft.ifft2(small_fft)
#==============================================================================
#     small_image = abs(small_image) / N*N
#==============================================================================

    for i in range(I//N):
        for j in range(J//N):
            if np.real(small_image[i,j]) < 0:
                small_image[i,j] = 0 + 0j
            else:
                small_image[i,j] /= (N*N)
    
    small_image = abs(small_image)
    
    t1 = time.time()
    print('decimate_FD time: ', t1-t0)
    
    return small_image