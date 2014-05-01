from math import pi as PI
from astropy.io import fits
import numpy
from pylab import *
import glob
import scipy.stats
import os

def HDUCaller(Filelist,name):

    Events=[]
    Elist=[]

    for filename in Filelist:                       # Append data to initially empty numpy arrays; two-step process to avoid data being overwritten
        hdulist = fits.open(filename)
        tbdata= hdulist[1].data                     
        evts=tbdata.field(name) 
        Events = append(Events, evts)               # Appends data to array
        if name == "L":
            Events = Events - 360.*(Events>180.)
        hdulist.close()
    return Events


def ROI(Filelist,name):
    Events = HDUCaller(Filelist,name)               # Calls HDUCaller for events.
    GalL=[]
    GalB=[]
    Llist=[]
    Blist=[]
    Elist=[]
    Zlist =[]
    Zenith=[]
        
    r=20                                            # From Weniger paper
    phi=60                                          # From Weniger paper
        
    for filename in Filelist:                       # Append data to initially empty numpy arrays; two-step process to avoid data being overwritten
        hdulist = fits.open(filename)
        tbdata= hdulist[1].data
##
##        gall= tbdata.field('L')                     # Redundant when name == "L" or name == "B", but necessary for equations to create ROI
##        GalL = append(GalL, gall)
##        GalL= GalL - 360.*(GalL>180.)               # Removes discontinuity at 0 boundary if L is greater than 180.
##
##        
##        galb= tbdata.field('B')   
##        GalB = append(GalB, galb)

        zen = tbdata.field("ZENITH_ANGLE")
        Zenith=append(Zenith,zen)

        hdulist.close()                             # Closes fits file

               
#    for (x,y,z,angle) in zip(GalL,GalB, Events, Zenith):          # Equations from Weniger paper to produce hourglass ROI
#           if x**2 + y**2 <= r**2 and \
#              x**2 <= y**2*math.tan(60*PI/180)**2+9 and \
#              angle >=100:
    for(event,angle) in zip(Events, Zenith):
        if angle >=100:
            Elist.append(event)                      # Appends Events that fall within ROI
    return Elist

######################################              # End of ROI 


allFiles=glob.glob('*.fits')

Energy1=HDUCaller(allFiles,"ENERGY")                # Energy all data
E1=ROI(allFiles,"ENERGY")                           # Energy in ROI
#bincount(E1, bin1)                                  # bincount for Energy in ROI
Z1=HDUCaller(allFiles,"ZENITH_ANGLE")               # Zenith Angle all data
L1=ROI(allFiles,"L")                                # L in ROI
B1=ROI(allFiles,"B")                                # B in ROI



f=figure()
bins = numpy.arange(0,185,5)
hist(Z1,bins) 
title('Zenith Angle')
xlabel("Angle (degrees)")
ylabel("Count")
savefig('Zenith Angle')

f=figure()
bins = numpy.arange(0,305000, 5000)
hist(E1,bins) 
title('Energy')
xlabel("Energy MeV")
ylabel("Count")
savefig('Energy')

show()
