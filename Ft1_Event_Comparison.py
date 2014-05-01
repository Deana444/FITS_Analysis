# Example program to analysis LAT data from a FITS file, using Python and Matplotlib
# Takes different .fits file(s) and compares them to each other.
# Checks to see which events are the same/different in the files

from math import pi as PI
from astropy.io import fits
import numpy
from pylab import *
import glob
import scipy.stats
import os

############################

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
######################################              # End of HDUCaller

def ROI(Filelist,name):
    Events = HDUCaller(Filelist,name)               # Calls HDUCaller for events.
    GalL=[]
    GalB=[]
    Llist=[]
    Blist=[]
    Elist=[]
        
    r=20                                            # From Weniger paper
    phi=60                                          # From Weniger paper
        
    for filename in Filelist:                       # Append data to initially empty numpy arrays; two-step process to avoid data being overwritten
        hdulist = fits.open(filename)
        tbdata= hdulist[1].data

        gall= tbdata.field('L')                     # Redundant when name == "L" or name == "B", but necessary for equations to create ROI                    
        GalL = append(GalL, gall)
        
        galb= tbdata.field('B')   
        GalB = append(GalB, galb) 

        hdulist.close()                             # Closes fits file

        GalL= GalL - 360.*(GalL>180.)               # Removes discontinuity at 0 boundary if L is greater than 180.
        
    for (x,y,z) in zip(GalL,GalB, Events):          # Equations from Weniger paper to produce hourglass ROI
           if x**2 + y**2 <= r**2 and \
              x**2 <= y**2*math.tan(60*PI/180)**2+9:
               Elist.append(z)                      # Appends Events that fall within ROI
               Llist.append(x)
               Blist.append(y)

    return Elist

######################################              # End of ROI 

def bincount(Elist,bins):
    count = histogram(Elist,bins)                   # Reports bin counts in Energy(GeV) ROI; 80-200GeV intervals
    (row,column) = count                            # Separates the matrix, row refers to the count and column revers to the energy intervals
    gev=[]                                          # Empty array for GeV values
    for item in column:
        gev.append(float(item)/1000)                # Converts from MeV to GeV
    print "\n\n\nE bins [GeV]:\t#Counts\n"
    for (count,bins) in zip(row, gev):
           print bins,"-",bins+(5.),"\t",count      # Prints the Bin Count in 5GeV intervals

    f=figure()
    scatter(gev[:-1],row)
    xlabel("Energy GeV")
    ylabel("Counts")
    show()

    

######################################              # End of bincount

def EvtMatch(Filelist1, Filelist2, events, runs):

    EID_1 = HDUCaller(Filelist1,events)
    EID_2 = HDUCaller(Filelist2,events)

    RID_1 = HDUCaller(Filelist1,runs)
    RID_2 = HDUCaller(Filelist2,runs)

    set1 = set(zip(EID_1,RID_1))
    set2 = set(zip(EID_2,RID_2))     
    overlap = list(set1.intersection(set2))

    return overlap

######################################              # End of EvtMatch

def UniqueEvts(Filelist,Filelist1, Filelist2, events, runs):
    set2 = set(EvtMatch(Filelist1, Filelist2, events, runs))
    
    EID = HDUCaller(Filelist,events)
    RID = HDUCaller(Filelist,runs)

    set1 = set(zip(EID,RID))
    unique = list(set1.difference(set2))
    
    return unique

######################################              # End of UniqueEvts

def FileList(Files):
    print "\nFiles:"
    for name in Files:
        print name
######################################              # End of FileList
    
binsize=5000
bin1=numpy.arange(80000,205000,binsize)


allFiles=glob.glob('*.fits')                        # gets files
FileList(allFiles)                                  # prints file names
Energy1=HDUCaller(allFiles,"ENERGY")                # Energy all data
E1=ROI(allFiles,"ENERGY")                           # Energy in ROI
bincount(E1, bin1)                                  # bincount for Energy in ROI
T1=ROI(allFiles,"TIME")                             # MET in ROI
MET1=HDUCaller(allFiles,"TIME")                     # MET all data
Z1=HDUCaller(allFiles,"ZENITH_ANGLE")               # Zenith Angle all data
CT1=HDUCaller(allFiles,"CONVERSION_TYPE")           # CT all data
L1=ROI(allFiles,"L")                                # L in ROI
B1=ROI(allFiles,"B")                                # B in ROI



#
# Asks user input for .fits folder comparison
#
##print "\n\nSecond Comparison\nEnter complete path to folder with .fits files:\n"
##test = False
##while test== False:
##    pathname=raw_input()
##    if pathname[len(pathname)-1] != "/":
##        pathname = pathname + "/"
##    if pathname[0]!="/":
##        pathname="/"+pathname
##    test = os.path.isdir(pathname)
##    if test == False:
##        print "That path does not exist. Reenter path:\n"
##    allFiles=glob.glob(pathname+"*.fits")    
##    if test == True and len(allFiles)==0:
##        print "No .fits files found.\nEnter different path: "
##        test == False

#newFiles=glob.glob("/Users/deanadelvecchio/Desktop/FERMI/Data/Galactic Center/Pass8 Astro Server/Complete Range/Feb6/AstroServer00004/filtered/"+"*.fits")
newFiles=glob.glob("/Users/deanadelvecchio/Desktop/FERMI/Data/Galactic Center/Pass8 Astro Server/43 Month Range/Feb14/Filtered/"+"*.fits")

FileList(newFiles)                                  # prints file names
Energy2=HDUCaller(newFiles,"ENERGY")                # Energy all data
E2=ROI(newFiles,"ENERGY")                           # Energy in ROI
bincount(E2,bin1)                                   # bincount for Energy in ROI
T2=ROI(newFiles,"TIME")                             # MET in ROI
MET2=HDUCaller(newFiles,"TIME")                     # MET all data
Z2=HDUCaller(newFiles,"ZENITH_ANGLE")               # Zenith Angle all data
CT2=HDUCaller(newFiles,"CONVERSION_TYPE")           # CT all data
L2=ROI(newFiles,"L")                                # L in ROI
B2=ROI(newFiles,"B")                                # B in ROI


Overlap=EvtMatch(allFiles,newFiles,"EVENT_ID","RUN_ID") # Looks for matching events
(events, runs) = zip(*Overlap)

Unique1 = UniqueEvts(allFiles,allFiles,newFiles,"EVENT_ID", "RUN_ID")
(evt_un1, run_un1) = zip(*Unique1)

Unique2 = UniqueEvts(newFiles,allFiles,newFiles,"EVENT_ID", "RUN_ID")
(evt_un2, run_un2) = zip(*Unique2)

print "Overlapping events in Pass7 and Pass8: ",len(events)
print "Unique events in Pass7: ",len(evt_un1)
print "Unique events in Pass8: ",len(evt_un2)
                                                    # Graphs Created

##f=figure()
##scatter(events,runs,0.1)
##title("Overlapping Events in Pass7 and Pass8")
##xlabel("Event_ID")
##ylabel("Run_ID")
##savefig('Overlap Events')
##
##f=figure()
##scatter(evt_un1,run_un1,0.1, label='Pass7')
##title("Unique Events")
##xlabel("Event_ID")
##ylabel("Run_ID")
##legend()
##savefig('Unique Events P7')
##
##f=figure()
##scatter(evt_un2,run_un2,0.1, label='Pass8')
##title("Unique Events")
##xlabel("Event_ID")
##ylabel("Run_ID")
##legend()
##savefig('Unique Events P8')

f=figure()
scatter(T1,E1,0.1, color='red', label='Pass7')
scatter(T2,E2,0.1,color='blue', label='Pass8')
title("Energy Vs Time in ROI")
xlabel("Time")
ylabel("Energy")
legend()
savefig('Energy vs Time Comparison ROI')

f=figure()
hist([E1,E2],bin1, color=['red','blue'], label=['Pass7','Pass8'], rwidth=1,histtype='bar')
title('Energy Comparison in ROI 80-200GeV')
xlabel("Energy (GeV)")
ylabel("Count")
legend()
savefig('Energy Comparison ROI')

f=figure()
bin2=numpy.arange(T1[0],T1[len(T1)-1],2000000)
hist([T1,T2],bin2, color=['red','blue'], label=['Pass7','Pass8'], rwidth=1,histtype='bar')
title('MET Comparison ROI')
xlabel("Time (s)")
ylabel("Count")
savefig('MET Comparison ROI')
legend()

f=figure()
hist([Z1,Z2],50, color=['red','blue'], label=['Pass7','Pass8'], rwidth=1,histtype='bar')
title('Zenith Angle Comparison')
xlabel("Angle (degrees)")
ylabel("Count")
legend()
savefig('Zenith Angle Comparison')

f=figure()
hist([CT1,CT2],2, color=['red','blue'], label=['Pass7','Pass8'], rwidth=1,histtype='bar')
xlabel("Conversion Type")
ylabel("Count")
xticks(np.arange(0,2,1),('front=0', 'back=1'))
title('Conversion Type Comparison')
legend()
savefig('Conversion Type Comparison')

f=figure()
bin5=numpy.arange(MET1[0], MET1[len(MET1)-1],2000000)
hist([MET1,MET2],bin5, color=['red','blue'], label=['Pass7','Pass8'], rwidth=1,histtype='bar')
title('MET Comparison')
xlabel("Time (s)")
ylabel("Count")
savefig('MET Comparison (all)')
legend()

f=figure()
scatter(L1,B1,0.1) # creates hourglass B vs L ROI
xlabel('l')
ylabel('b')
title('Pass7 B vs L ROI')
savefig('Pass7 B vs L ROI')

f=figure()
scatter(L2,B2,0.1) # creates hourglass B vs L ROI
xlabel('l')
ylabel('b')
title('Pass8 B vs L ROI')
savefig('Pass8 B vs L ROI')

f=figure()
hist(E1,bin1, color='red', label='Pass7')
title('Energy in ROI 80-200GeV')
xlabel("Energy (GeV)")
ylabel("Count")
legend()
savefig('Energy Pass7 ROI')

f=figure()
hist(E2,bin1, color='blue', label='Pass8')
title('Energy in ROI 80-200GeV')
xlabel("Energy (GeV)")
ylabel("Count")
legend()
savefig('Energy Pass8 ROI')

show()


