#!/usr/bin/env python

import argparse
import pylab
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

args = parser.parse_args()

reffiles=args.ireference
samfiles=args.isample
mode=args.mode
thickness=args.thickness

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

mdata=Terapy.teradata(ref_fd,sam_fd)

mdata.manipulateFDData(-11e9,[300e9,2.2e12])

myana=Terapy.teralyz(mdata,thickness,20e-6,30)
myana.doCalculation(args.calcLength,args.noSVMAF,args.silent)
#
#myana.plotInits(mdata.H, thickness)

if args.outname==None:
    args.outname=myana.getFilenameSuggestion()

args.outname+='_'

if args.savePlots:
    pylab.ioff()
    reftd.doPlotWithunc()
    samtd.doPlotWithunc()
    pylab.legend(('Reference','Sample'))
    pylab.savefig(args.outname + 'Time-Domain.png')
    pylab.close()
    
    ref_fd.doPlot()
    sam_fd.doPlot()
    pylab.figure('FD-ABS-Plot')
    pylab.legend(('Reference','Sample'))
    pylab.savefig(args.outname + 'ABS-Frequency-Domain.png')
    pylab.close()
    pylab.figure('FD-PHASE-Plot')
    pylab.legend(('Reference','Sample'))
    pylab.savefig(args.outname + 'PHASE-Frequency-Domain.png')
    pylab.close()
#    
    mdata.doPlots()
    pylab.savefig(args.outname + 'TransferFunction.png')
    pylab.close()

myana.plotRefractiveIndex(1,1,args.outname)
myana.saveResults(args.outname)

pylab.show()