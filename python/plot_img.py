"""
/*
QSTEM - image simulation for TEM/STEM/CBED
    Copyright (C) 2000-2010  Christoph Koch
    Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
"""


import struct
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from os.path import splitext
from pyE17 import utils as U
def binread2D(filename, printFlag=True):
    comment=""
    # open the file and define the file ID (fid):
    f=open(filename,"rb")
    #with open(filename,'rb') as f:
    if True:
        # there are 8 ints followed by 3 doubles
        header = struct.unpack("iiiiiiiiddd",f.read(56))

        print header

        headerSize = header[0]

        # read additional parameters from file, if any exist:
        paramSize = header[1]

        commentSize = header[2]

        Nx = header[3]
        Ny = header[4]

        complexFlag = header[5]
        dataSize = header[6]
        doubleFlag = (dataSize==8*(complexFlag+1))
        complexFlag = bool(complexFlag)

        version = header[7]
        
        thicknessOrDefocus=header[8]

        dx = header[9]
        dy = header[10]

        if (paramSize > 0):
            params = np.fromfile(file=f, dtype=np.float64, count=paramSize);
            if printFlag:
                print '%d Parameters:'%paramSize
                print params

        # read comments from file, if any exist:
        if (commentSize > 0):
            comment = struct.unpack("%ds"%commentSize,f.read(commentSize))[0]

    if printFlag:
        print('binread2D %s: %d x %d pixels'%(filename,Nx,Ny))

    if complexFlag:
        if doubleFlag:
            if printFlag:
                print '64-bit complex data, %.3fMB)\n'%(Nx*Ny*16/1048576)
            img = np.fromfile(file=f, dtype=np.complex128, count=Nx*Ny)
        else:
            if printFlag:
                fprintf('32-bit complex data, %.3fMB)\n',Nx*Ny*8/1048576);
            img = np.fromfile(file=f, dtype=np.complex64, count = Nx*Ny)
    else:
        if doubleFlag:
            if printFlag:
                fprintf('64-bit real data, %.3fMB)\n',Nx*Ny*8/1048576);
            img = np.fromfile(file=f, dtype=np.float64, count=Nx*Ny)
        else:
            if printFlag:
                print '32-bit real data, %.3fMB)\n'%(Nx*Ny*4/1048576)
            img = np.fromfile(file=f, dtype=np.float32, count=Nx*Ny)
    img=img.reshape(Ny,Nx)
    
    return img, comment, thicknessOrDefocus, dx, dy
    
def plot_img(filename):
    print 'saving file ' + filename
    img, comment, t, dx, dy = binread2D(filename, False)
    fig = plt.figure()
    fig.suptitle(comment, fontsize=14, fontweight='bold')
    ax=fig.add_subplot(111)
    extent = [0, img.shape[0]*dx, 0, img.shape[1]*dy]
    ax.imshow(U.imsave(img), extent=extent, interpolation="nearest")
    ax.set_title("Thickness = %.3fA"%t)
    ax.set_xlabel("Angstroms")
    ax.set_ylabel("Angstroms")
    plt.savefig(splitext(filename)[0]+".png",bbox_inches=0,)

def plot_diff(filename):
    print 'saving file ' + filename
    img, comment, t, dx, dy = binread2D(filename, False)
    fig = plt.figure()
    fig.suptitle(comment, fontsize=14, fontweight='bold')
    ax=fig.add_subplot(111)
    extent = [0, img.shape[0]*dx, 0, img.shape[1]*dy]
    i = np.fft.fftshift(np.fft.fft2(img))
    intens = np.real(i*np.conj(i))
    cax = ax.imshow(np.log10(intens), extent=extent, interpolation="nearest", cmap=plt.get_cmap('binary'))
    cbar = fig.colorbar(cax)
    ax.set_title("Thickness = %.3fA"%t)
    ax.set_xlabel("Angstroms")
    ax.set_ylabel("Angstroms")
    plt.savefig(splitext(filename)[0]+"_diff.png",bbox_inches=0,)
def plot_imgs(path):
    print path+'*.img'
    files = glob.glob(path+'/*.img')
    for f in files:
        if 'wave' in f:
            plot_img(f)
        #plot_diff(f)
#        elif 'diff' in f:
#            plot_diff(f)   
            
if __name__=="__main__":
    import sys
    import glob
    plot_imgs(sys.argv[1])
    
     
        
            
    
    