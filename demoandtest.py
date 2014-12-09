import Terapy
import TeraData
import glob
import pylab as py


#provide a path to the reference and sample files
#a list of strings or a glob expression can be used
samfiles=glob.glob('/home/jahndav/Dropbox/THz-Analysis/rehi/Sample_?.txt')
reffiles=glob.glob('/home/jahndav/Dropbox/THz-Analysis/rehi/Ref*.txt')

#use the appropriate importer for the Data,
#initialize the TimeDomain Data Object with a list of filenames
sam_td=TeraData.ImportMarburgData(samfiles)
ref_td=TeraData.ImportMarburgData(reffiles)

#the Frequency domain data is calculated by passing a TimeDomain data object to the FdData class
sam_fd=TeraData.FdData(sam_td)
ref_fd=TeraData.FdData(ref_td)

#the Transferfunction is a object of type HMeas, which is nothing than a special frequency domain 
#data type, such that HMeas inherits from FdData
H=Terapy.HMeas(ref_fd,sam_fd)

##################################################################################################
# until here was a demonstration how the objects can be initialized
##################################################################################################

#the time domain data can be manipulatd through the member functions of THzTdData, so for instance
#windowing:
sam_td.setTDData(sam_td.getWindowedData(1e-12))
ref_td.setTDData(ref_td.getWindowedData(1e-12))

#create new frequency domain data, with the windowed timedomain data as input
sam_fd_windowed=TeraData.FdData(sam_td)
ref_fd_windowed=TeraData.FdData(ref_td)

H_windowed=Terapy.HMeas(ref_fd_windowed,sam_fd_windowed)

#Zeropadding in time domain
sam_td.zeroPaddData(sam_td.getLength()*2)
ref_td.zeroPaddData(ref_td.getLength()*2)

#create new FdData array
sam_fd_zpandwind=TeraData.FdData(sam_td)
ref_fd_zpandwind=TeraData.FdData(ref_td)

H_zpandwind=Terapy.HMeas(ref_fd_zpandwind,sam_fd_zpandwind)

#do some plots
#py.figure(1)
H.doPlot()
H_windowed.doPlot()
H_zpandwind.doPlot()
py.figure('H-PHASE-Plot')
py.legend(('orig','windowed','w+zp'))

py.figure('H-ABS-Plot')
py.legend(('orig','windowed','w+zp'))

#py.plot(sam_fd.getfreqsGHz(),abs(sam_fd.getFAbs()-sam_fd_windowed.getFAbs())/sam_fd.getFAbs())
#py.plot(sam_fd.getfreqsGHz(),abs(sam_fd.getFPh()-sam_fd_windowed.getFPh())/sam_fd.getFPh())
#py.xlim((0,400))
#py.plot(sam_fd.getfreqsGHz(),sam_fd.getFAbs()-sam_fd_zpandwind.getFAbs())
