#!/usr/bin/env python

### ToDo
'''  A script to fast analyze relevant measurement parameters, without calculating n, for comparism of different spectrometer
The script should work for single input files, and also for a series of sample and reference measurements
runGetSignalInfos should print:
    maxDR of fourier and timedomain data
    max SNR of fourier and timedomain data
    Bandwidth, highest and lowest accesible frequency
    if sample and reference measurements are specified this should be done for both of them, and also for H
runGetSingalInfos should plot:
    TD Plot along with uncertainties
    FD Plot along with uncertainties
    SNR Plot in TD and frequency domain
    DR Plot in TD and frequency domain
    if sample and reference measurements are specified this should be done for both of them, and also for H Comparison option
Maybe it would be also favorable to specify a series of reference files from one spectrometer (i.e. Marburg ,or INRIM on 20.10.2014 and compare it with another series of reference files)
'''
### Done
'''
SNR Plot in FD&TD, DR Plot in FD&TD, TD with u
'''

import argparse
import sys
import matplotlib.pyplot as plt
import glob
import Terapy
import TeraData

parser = argparse.ArgumentParser(description='Calculate optical constants from THz-TD Data')
parser.add_argument('--silent',action='store_true',help='switch save results off')
parser.add_argument('--outname','-o',nargs='?',type=str,help='prefix output filenames')
parser.add_argument('--isample','-is',nargs='*',help='list of sample filenames')
parser.add_argument('--ireference','-ir',nargs='*',help='list of reference filenames')
parser.add_argument('--mode','-m',type=str,default='INRIM',choices=['INRIM','Marburg','lucastestformat'],help='format of the datafiles')
parser.add_argument('--thickness','-t',type=float,help='sample thickness')
parser.add_argument('--savePlots','-s',action='store_false',help='turn off saving TD and FD plots')
parser.add_argument('--workpath','-w',type=str,default='./',help='specify a base folder')
args = parser.parse_args()

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
H=Terapy.HMeas(ref_fd,sam_fd)

fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.set_xlabel('frequency (GHz)')
ax.set_ylabel('SNR')
ax.grid(True)
ax.semilogy(ref_fd.getfreqsGHz(),ref_fd.getSNR(),sam_fd.getfreqsGHz(),sam_fd.getSNR(),H.getfreqsGHz(),H.getSNR())
ax.legend(('ref', 'sam', 'transfer function'))
plt.title('SNR')

#fig2 = plt.figure()
ax2 = fig.add_subplot(2,1,2)
ax2.set_xlabel('time (ps)')
ax2.set_ylabel('SNR')
ax2.grid(True)
ax2.semilogy(reftd.getTimesPs(),reftd.getSNR(),samtd.getTimesPs(),samtd.getSNR())
ax2.legend(('ref', 'sam'))
#plt.title('SNR')


fig1 = plt.figure()
ax1 = fig1.add_subplot(2,1,1)
ax1.set_xlabel('frequency (GHz)')
ax1.set_ylabel('dynamic range')
ax1.grid(True)
ax1.semilogy(ref_fd.getfreqsGHz(),ref_fd.getDR(),sam_fd.getfreqsGHz(),sam_fd.getDR(),H.getfreqsGHz(),H.getDR())
ax1.legend(('ref', 'sam', 'transfer function'))
plt.title('dynamic range')

#fig3 = plt.figure()
ax3 = fig1.add_subplot(2,1,2)
ax3.set_xlabel('time (ps)')
ax3.set_ylabel('dynamic range')
ax3.grid(True)
ax3.semilogy(reftd.getTimesPs(),reftd.getDR(),samtd.getTimesPs(),samtd.getDR())
ax3.legend(('ref', 'sam'))
#plt.title('dynamic range')

fig2 = plt.figure()
ax4 = fig2.add_subplot(2,1,1)
ax4.set_xlabel('time (ps)')
ax4.set_ylabel('X channel (V)')
#ax4.grid(True)
no_std=2
ax4.plot(reftd.getTimesPs(),reftd.getEX(),\
    reftd.getTimesPs(),reftd.getEX() + no_std*reftd.getUncEX(),'g--',\
    reftd.getTimesPs(),reftd.getEX() - no_std*reftd.getUncEX(),'g--')
#ax4.legend(('ref'))
plt.title('reference signal with uncertainty')

#fig5 = plt.figure()
ax5 = fig2.add_subplot(2,1,2)
ax5.set_xlabel('time (ps)')
ax5.set_ylabel('X channel (V)')
#ax4.grid(True)
no_std=2
ax5.plot(samtd.getTimesPs(),samtd.getEX(),\
    samtd.getTimesPs(),samtd.getEX() + no_std*samtd.getUncEX(),'g--',\
    samtd.getTimesPs(),samtd.getEX() - no_std*samtd.getUncEX(),'g--')
#ax5.legend(('sam'))
plt.title('sample signal with uncertainty')

fig3 = plt.figure()
ax6 = fig3.add_subplot(2,1,1)
ax6.set_xlabel('frequency (GHz)')
ax6.set_ylabel('dynamic range')
ax6.grid(True)
ax6.semilogy(ref_fd.getfreqsGHz(),ref_fd.getFAbs(),\
    ref_fd.getfreqsGHz(), ref_fd.getFAbs() + ref_fd.getFAbsUnc(), 'g--',\
        ref_fd.getfreqsGHz(), ref_fd.getFAbs() - ref_fd.getFAbsUnc(), 'g--')
#ax6.legend(('ref'))
plt.title('ABS with U')

ax7 = fig3.add_subplot(2,1,2)
ax7.set_xlabel('frequency (GHz)')
ax7.set_ylabel('dynamic range')
ax7.grid(True)
ax7.semilogy(ref_fd.getfreqsGHz(),ref_fd.fdData[:,3],\
    ref_fd.getfreqsGHz(), ref_fd.getFPh() + ref_fd.getFPhUnc(), 'g--',\
        ref_fd.getfreqsGHz(), ref_fd.getFPh() - ref_fd.getFPhUnc(), 'g--')
#ax7.legend(('ref'))
plt.title('PHASE with U')

fd = file( basefolder + 'SignalInfo.log', 'w')
fd.write('max DR in FD - ref\t' + str(max(ref_fd.getDR())) + '\n'\
    'max DR in FD - sam\t' + str(max(sam_fd.getDR())) + '\n'\
    'max DR in TD - ref\t' + str(max(reftd.getDR())) + '\n'\
    'max DR in TD - sam\t' + str(max(samtd.getDR())) + '\n\n'\
    'max SNR in FD - ref\t' + str(max(ref_fd.getSNR())) + '\n'\
    'max SNR in FD - sam\t' + str(max(sam_fd.getSNR())) + '\n'\
    'max SNR in TD - ref\t' + str(max(reftd.getSNR())) + '\n'\
    'max SNR in TD - sam\t' + str(max(samtd.getSNR())) + '\n')
fd.close()

'''    maxDR of fourier and timedomain data
    max SNR of fourier and timedomain data
    Bandwidth, highest and lowest accesible frequency
'''

'''if args.outname==None:
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
#'''
plt.show()
