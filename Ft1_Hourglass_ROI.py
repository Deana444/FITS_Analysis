# Example program to analysis LAT data from a FITS file, using Python and Matplotlib
# Compares HDU information against Weniger paper

from math import pi as PI
from astropy.io import fits
import numpy
from pylab import *
import glob
import scipy.stats

r=20 # From Weniger paper
phi=60 # From Weniger paper
bsize=5000.# histogram bin size in MeV

allFiles=glob.glob('*.fits') # list of fits file(s) being used for analysis
print "Files:\n"
for name in allFiles:
       print name

### Empty arrays for each .fits data to be appended to
Energy=[]
GalL=[]
GalB=[]
Llist=[]
Blist=[]
Elist=[]
Time=[]
Tlist=[]

for item in allFiles:        # Append data to initially empty numpy arrays; two-step process to avoid data being overwritten

       hdulist = fits.open(item)
       tbdata= hdulist[1].data #from .fits file(s)

       energy= tbdata.field('ENERGY') # takes energy data from HDU
       Energy = append(Energy, energy)# appends data to array

       time= tbdata.field('TIME')
       Time = append(Time, time)
       
       gall= tbdata.field('L')   # takes L data from HDU
       GalL = append(GalL, gall) # appends data to array
       
       galb= tbdata.field('B')   # takes B data from HDU
       GalB = append(GalB, galb) # appends data to array

       GalL= GalL - 360.*(GalL>180.) # removes discontinuity at 0 boundary if L is greater than 180.
       for (x,y,z,w) in zip(GalL,GalB, Energy,Time): # equations from Weniger paper to produce hourglass ROI
              if x**2 + y**2 <= r**2:
                     if x**2 <= y**2*math.tan(60*PI/180)**2+9:
                            Llist.append(x) # appends L, B, Energy values that fall within ROI
                            Blist.append(y)
                            Elist.append(z)
                            Tlist.append(w)
       hdulist.close() # closes fits file

f=figure()
scatter(Llist,Blist,0.1) # creates hourglass B vs L ROI
xlabel('l')
ylabel('b')
savefig('B vs L')

f=figure()
bin1 = numpy.arange(Time[0], Time[len(Time)-1],2000000)
hist(Time, bin1) # histogram in 60 second intervals
xlabel('MET')
savefig('Mission Elapsed Time')

f=figure()
scatter(Tlist,Elist,0.1)
xlabel("Time")
ylabel("Energy")
savefig('Energy vs Time')

f=figure()
bin1 = numpy.arange(20000,300000,bsize) # Histogram of energies in 20-300GeV ROI
hist(Elist, bin1)
xlabel('Energy MeV')
savefig('Hourglass Energy.png')

f=figure()
bin2 = numpy.arange(80000,200000,bsize) # Histogram of energies in hourglass 80-200GeV ROI
hist(Elist, bin2)
xlabel('Energy MeV')
savefig('80-200 Energy.png')

show() # shows the figures

count = histogram(Elist,bin2) # outputs the count of the bin in Energy(GeV) hourglass histogram; 80-200GeV intervals
(row,column) = count # separates the matrix, row refers to the count and column revers to the energy intervals
gev=[] # empty array for GeV values
for item in column:
    gev.append(float(item)/1000) # fills aray with values; converts from MeV to GeV
print "\n\n\nE bins [GeV]:\t#Counts\n"
for (count,bins) in zip(row, gev):# for loop that prints the counts in each bin
       print bins,"-",bins+(bsize/1000.),"\t",count #prints the bincount in 5GeV intervals

