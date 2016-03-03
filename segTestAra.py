###################################################### license
#
################
import numpy as np
import scipy as sp
from scipy import ndimage as ndi
import matplotlib.pyplot as plt

#from skimage import filters

from skimage.filters import sobel
from skimage.segmentation import slic, join_segmentations
from skimage.morphology import watershed
from skimage.color import label2rgb
from skimage import data, img_as_float
from scipy import ndimage
from skimage import restoration as rest

from skimage import io as sio

### P4D ###
import multiprocessing
import numpy as np
import os
import sys
import time
import getopt

sys.path.append('P4D_SourceCode')
import P4D_help                                 # implementations of various auxiliary functions
reload(P4D_help)
#"""import P4D_segment                              # function for segmentation of leafs and plants and computation of positions and shape factors
#reload(P4D_segment)
#import P4D_analyze                              # function for tracking of leaves, leaf number assignment and export of leaf-specific time series data
#reload(P4D_analyze)
#import P4D_track                                # implementation of Crocker-Grier tracking algorithm by Thomas A. Caswell under GLP3
#reload(P4D_track)
#import P4D_plot                                 # function for automated plotting of plant and leaf features for one plant
#reload(P4D_plot)
#import P4D_multi                                # function for generation of customized summary plots including multiple plants and genotypes
#reload(P4D_multi)
#"""
import scipy as sp
import skimage

import P3D_segment 
reload(P3D_segment)


def main(argv):
   ffile = ''
   nfile = ''
   pfile = ''
   optlist, args = getopt.getopt(sys.argv[1:],':p:n:f:h')
   for opt, par in optlist:
      if opt == '-h':
         print 'python segtestAra.py -f <Focus> -n <Nums> -p <Pics>'
         sys.exit()
      elif opt == "-f" :
         ffile = par
      elif opt == "-n" :
         nfile = par
      elif opt == "-p" :
         pfile = par
      else :
         raise StandardError, 'error "%s"' % opt

   print 'Input file is "', ffile
   print 'Output file is "', nfile
   print 'Output file is "', pfile
   segmentLeaf(ffile, nfile, pfile)


def segmentLeaf(f,n,p):
   #folder = os.path.join(os.getenv("HOME"), 'repo/leaftrack');

  # Set parameters
   
   dirF=os.path.join(f)
   dirN=os.path.join(n)
   dirP=os.path.join(p)

   SegBGupper=25                                   # upper threshold for image background
   SegFGlower=80                                   # lower threshold for image foreground
   SegSigmaGauss=50.0                              # sigma of Gaussian filter for smoothing depth image
   SegSigmaCanny=1.0                               # sigma of Canny filter for detection of leaf edges             
   SegThresSlope=0.1                               # slope of radially increasing threshold to find watershed seeds from distance transformed image
   SegThresAbsci=10.0                              # abscissa of radially increasing threshold to find watershed seeds from distance 
   SegRadiusOriginHeight= 2.0                     # size of disk for determining height of origin in pixels
   SegRadiusLeafHeight=5.0                        # size of disk for determining height of leaf positions in pixels
   SegRadiusStemEraser=5.0                        # size of disk for removing leaf petioles

   for i in [dirN,dirP]:                                 # generation output folders if not present
       if(not os.path.isdir(i)): os.mkdir(i)

   for id1 in range(10):
      imid = 'test'+str(id1+1)
      imR = sio.imread(os.path.join(dirF, imid+'.png')) # Raw image
      imF = sio.imread(os.path.join(dirF, imid+'.png'), as_grey = True) # grey image
      imF = imF*255 #image at grey scale, value 0-255    
      imD = np.ones(imF.shape)  # No depth image
      SegBGupper = np.percentile(imF, 75)
      ly, lx = imF.shape
      CropRadius = min(ly, lx)/2 - 10
      inp= imid, imR, imF, imD, dirP, dirN, SegBGupper,SegFGlower, CropRadius, SegSigmaGauss,SegSigmaCanny,SegThresSlope,SegThresAbsci,SegRadiusOriginHeight,SegRadiusLeafHeight,SegRadiusStemEraser
      paramFile = P3D_segment.segment(inp)
    #para=np.load(paramFile)

if __name__ == "__main__":
   #import sys
   main(sys.argv[:])

