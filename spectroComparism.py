#!/usr/bin/env python

'''
SNR Plot in FD&TD, DR Plot in FD&TD, TD with u for reference files from two different spectrometers/runs
'''

import argparse
import sys
import matplotlib.pyplot as plt
import glob
import Terapy
import TeraData
from matplotlib2tikz import save as tikz_save

plt.rc('text',usetex=True)
parser = argparse.ArgumentParser(description='Do a spectrometer comparism')
parser.add_argument('--outname','-o',nargs='?',type=str,help='prefix output filenames')
parser.add_argument('--ireference1','-ir1',nargs='*',help='list of reference filenames, dataset 1')
parser.add_argument('--ireference2','-ir2',nargs='*',help='list of reference filenames, dataset 2')
parser.add_argument('--mode1','-m1',type=str,default='INRIM',choices=['INRIM','Marburg','lucastestformat'],help='format of the datafiles of dataset 1')
parser.add_argument('--mode2','-m2',type=str,default='INRIM',choices=['INRIM','Marburg','lucastestformat'],help='format of the datafiles of dataset 2')
parser.add_argument('--workpath','-w',type=str,default='',help='specify a base folder')
args = parser.parse_args()

ireffiles1=args.ireference1
ireffiles2=args.ireference2
mode1=args.mode1
mode2=args.mode2
basefolder=args.workpath        

reffiles1=[]
reffiles2=[]

for i in range(len(ireffiles1)):
    tf=glob.glob(basefolder+ireffiles1[i])   
    reffiles1+=tf

for i in range(len(ireffiles2)):
    tf=glob.glob(basefolder+ireffiles2[i])   
    reffiles2+=tf

if len(reffiles1)==0:
    print "no Reference File for spectrometer one specified"
    sys.exit()
    
if len(reffiles2)==0:
    print "no Reference File for spectrometer two specified"
    sys.exit()
        
if mode1=='lucastestformat':
    reftd1=TeraData.THzTdData(reffiles1)
elif mode1=='Marburg':
    reftd1=TeraData.ImportMarburgData(reffiles1)
elif mode1=='INRIM':
    reftd1=TeraData.ImportInrimData(reffiles1)


if mode2=='lucastestformat':
    reftd2=TeraData.THzTdData(reffiles2)
elif mode2=='Marburg':
    reftd2=TeraData.ImportMarburgData(reffiles2)
elif mode2=='INRIM':
    reftd2=TeraData.ImportInrimData(reffiles2)

#    #initialize the fd_data objects        
ref_fd1=TeraData.FdData(reftd1)
ref_fd2=TeraData.FdData(reftd2)

#shift to same time!
reftd1.tdData[:,0]-=reftd1.getPeakPosition()
reftd2.tdData[:,0]-=reftd2.getPeakPosition()

fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.set_xlabel('frequency (GHz)')
ax.set_ylabel('SNR')
ax.grid(True)
ax.semilogy(ref_fd1.getfreqsGHz(),ref_fd1.getSNR(),ref_fd2.getfreqsGHz(),ref_fd2.getSNR())
ax.legend((mode1, mode2))
plt.title('SNR')

#fig2 = plt.figure()
ax2 = fig.add_subplot(2,1,2)
ax2.set_xlabel('time (ps)')
ax2.set_ylabel('SNR')
ax2.grid(True)
ax2.semilogy(reftd1.getTimesPs(),reftd1.getSNR(),reftd2.getTimesPs(),reftd2.getSNR())
ax2.legend((mode1, mode2))
fig.savefig(basefolder+'SNR-Plot.png')
tikz_save(basefolder+'SNR-Plot.tikz',figureheight='\\figureheight',figurewidth='\\figurewidth')
#plt.title('SNR')


fig1 = plt.figure()
ax1 = fig1.add_subplot(2,1,1)
ax1.set_xlabel('frequency (GHz)')
ax1.set_ylabel('dynamic range')
ax1.grid(True)
ax1.semilogy(ref_fd1.getfreqsGHz(),ref_fd1.getDR(),ref_fd2.getfreqsGHz(),ref_fd2.getDR())
ax1.legend((mode1,mode2))
plt.title('dynamic range')

#fig3 = plt.figure()
ax3 = fig1.add_subplot(2,1,2)
ax3.set_xlabel('time (ps)')
ax3.set_ylabel('dynamic range')
ax3.grid(True)
ax3.semilogy(reftd1.getTimesPs(),reftd1.getDR(),reftd2.getTimesPs(),reftd2.getDR())
ax3.legend((mode1,mode2))
#plt.title('dynamic range')

fig1.savefig(basefolder+'DR-Plot.png')
tikz_save(basefolder+'DR-Plot.tikz',figureheight='\\figureheight',figurewidth='\\figurewidth')

 
fig2 = plt.figure()
ax4 = fig2.add_subplot(2,1,1)
ax4.set_xlabel('time (ps)')
ax4.set_ylabel('X channel (V)')
#ax4.grid(True)
no_std=2
ax4.plot(reftd1.getTimesPs(),reftd1.getEX(),\
    reftd1.getTimesPs(),reftd1.getEX() + no_std*reftd1.getUncEX(),'g--',\
    reftd1.getTimesPs(),reftd1.getEX() - no_std*reftd1.getUncEX(),'g--')
#ax4.legend(('ref'))
plt.title('Reference spectrometer ' + mode1+ ' with uncertainty')

#fig5 = plt.figure()
ax5 = fig2.add_subplot(2,1,2)
ax5.set_xlabel('time (ps)')
ax5.set_ylabel('X channel (V)')
#ax4.grid(True)
no_std=2
ax5.plot(reftd2.getTimesPs(),reftd2.getEX(),\
    reftd2.getTimesPs(),reftd2.getEX() + no_std*reftd2.getUncEX(),'g--',\
    reftd2.getTimesPs(),reftd2.getEX() - no_std*reftd2.getUncEX(),'g--')
#ax5.legend(('sam'))
plt.title('Reference spectrometer ' + mode2+ ' with uncertainty')
fig2.savefig(basefolder+'TDSignal-Plot.png')
tikz_save(basefolder+'TDSignal-Plot.tikz',figureheight='\\figureheight',figurewidth='\\figurewidth')

fig3 = plt.figure()
ax6 = fig3.add_subplot(2,1,1)
ax6.set_xlabel('frequency (GHz)')
ax6.set_ylabel('dynamic range')
ax6.grid(True)
ax6.semilogy(ref_fd1.getfreqsGHz(),ref_fd1.getFAbs(),\
    ref_fd1.getfreqsGHz(), ref_fd1.getFAbs() + ref_fd1.getFAbsUnc(), 'g--',\
        ref_fd1.getfreqsGHz(), ref_fd1.getFAbs() - ref_fd1.getFAbsUnc(), 'g--',
    ref_fd2.getfreqsGHz(),ref_fd2.getFAbs(),\
    ref_fd2.getfreqsGHz(), ref_fd2.getFAbs() + ref_fd2.getFAbsUnc(), 'g--',\
        ref_fd2.getfreqsGHz(), ref_fd2.getFAbs() - ref_fd2.getFAbsUnc(), 'g--')
#ax6.legend(('ref'))
plt.title('ABS with U')

ax7 = fig3.add_subplot(2,1,2)
ax7.set_xlabel('frequency (GHz)')
ax7.set_ylabel('dynamic range')
ax7.grid(True)
ax7.plot(ref_fd1.getfreqsGHz(),ref_fd1.getFPh(),\
    ref_fd1.getfreqsGHz(), ref_fd1.getFPh() + ref_fd1.getFPhUnc(), 'g--',\
        ref_fd1.getfreqsGHz(), ref_fd1.getFPh() - ref_fd1.getFPhUnc(), 'g--',
    ref_fd2.getfreqsGHz(),ref_fd2.getFPh(),\
    ref_fd2.getfreqsGHz(), ref_fd2.getFPh() + ref_fd2.getFPhUnc(), 'g--',\
        ref_fd2.getfreqsGHz(), ref_fd2.getFPh() - ref_fd2.getFPhUnc(), 'g--')

fig3.savefig(basefolder+'FDSignal.png')
tikz_save(basefolder+'FDSignal.tikz',figureheight='\\figureheight',figurewidth='\\figurewidth')

#ax7.legend(('ref'))
plt.title('PHASE with U')

fd = file( basefolder + 'SignalInfo.log', 'w')
fd.write('max DR in FD - ref1\t' + str(max(ref_fd1.getDR())) + '\n'\
    'max DR in FD - ref2\t' + str(max(ref_fd2.getDR())) + '\n'\
    'max DR in TD - ref1\t' + str(max(reftd1.getDR())) + '\n'\
    'max DR in TD - ref2\t' + str(max(reftd2.getDR())) + '\n\n'\
    'max SNR in FD - ref1\t' + str(max(ref_fd1.getSNR())) + '\n'\
    'max SNR in FD - ref2\t' + str(max(ref_fd2.getSNR())) + '\n'\
    'max SNR in TD - ref1\t' + str(max(reftd1.getSNR())) + '\n'\
    'max SNR in TD - ref2\t' + str(max(reftd2.getSNR())) + '\n')
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
