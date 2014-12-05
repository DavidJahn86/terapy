#!/usr/bin/env python
import time
import argparse
import sys
import pylab
import glob
import Terapy
import TeraData

parser = argparse.ArgumentParser(description='Calculate optical constants from THz-TD Data')
parser.add_argument('--calcLength',action='store_false',help='switch length calculation off')
parser.add_argument('--silent',action='store_true',help='switch save results off')
parser.add_argument('--noSVMAF',default=5,nargs='?',type=int,help='No of SVMAF iterations')
parser.add_argument('--outname','-o',nargs='?',type=str,help='prefix output filenames')
parser.add_argument('--isample','-is',nargs='*',help='list of sample filenames')
parser.add_argument('--ireference','-ir',nargs='*',help='list of reference filenames')
parser.add_argument('--mode','-m',type=str,default='INRIM',choices=['INRIM','Marburg','lucastestformat'],help='format of the datafiles')
parser.add_argument('--thickness','-t',type=float,help='sample thickness')
parser.add_argument('--savePlots','-s',action='store_false',help='turn off saving TD and FD plots')
parser.add_argument('--workpath','-w',type=str,default='',help='specify a base folder')
args = parser.parse_args()

starttime=time.time()
ireffiles=args.ireference
isamfiles=args.isample
mode=args.mode
thickness=args.thickness
basefolder=args.workpath        

reffiles=[]
samfiles=[]

for i in range(len(ireffiles)):
    tf=glob.glob(basefolder+ireffiles[i])   
    reffiles+=tf

for i in range(len(isamfiles)):
    tf=glob.glob(basefolder+isamfiles[i])
    samfiles+=tf

if len(reffiles)==0:
    print "no Reference File specified"
    sys.exit()
    
if len(samfiles)==0:
    print "no Sample File specified"
    sys.exit()
        
if mode=='lucastestformat':
    reftd=TeraData.THzTdData(reffiles)
    samtd=TeraData.THzTdData(samfiles)
    
if mode=='Marburg':
    reftd=TeraData.ImportMarburgData(reffiles)
    samtd=TeraData.ImportMarburgData(samfiles)

if mode=='INRIM':
    reftd=TeraData.ImportInrimData(reffiles)
    samtd=TeraData.ImportInrimData(samfiles)

#    #initialize the fd_data objects        
ref_fd=TeraData.FdData(reftd)
sam_fd=TeraData.FdData(samtd)

##    #initialize the mdata object (H,and so on)

mdata=Terapy.HMeas(ref_fd,sam_fd)

#mdata.manipulateFDData(-11e9,[200e9,2.2e12])

myana=Terapy.teralyz(mdata,thickness,20e-6,30)
myana.doCalculation(args.calcLength,args.noSVMAF,args.silent)


if args.outname==None:
    args.outname=myana.getFilenameSuggestion()

args.outname+='_'

if args.savePlots:
    pylab.ioff()
    reftd.doPlotWithunc()
    samtd.doPlotWithunc()
    pylab.legend(('Reference','Sample'))
    pylab.savefig(args.workpath+args.outname + 'Time-Domain.png')
    pylab.close()
    
    ref_fd.doPlot()
    sam_fd.doPlot()
    pylab.figure('FD-ABS-Plot')
    pylab.legend(('Reference','Sample'))
    pylab.savefig(args.outname + 'ABS-Frequency-Domain.png')
    pylab.close()
    pylab.figure('FD-PHASE-Plot')
    pylab.legend(('Reference','Sample'))
    pylab.savefig(args.workpath+args.outname + 'PHASE-Frequency-Domain.png')
    pylab.close()
    
    mdata.doPlot()
    pylab.savefig(args.workpath+args.outname + 'TransferFunction.png')
    pylab.close()
#
myana.plotRefractiveIndex(1,1,args.workpath+args.outname)
myana.saveResults(args.workpath+args.outname)
#
endtime=time.time()
print "Consumed Time: " + str(endtime-starttime)
pylab.show()
