import Terapy
import TeraData
import glob
import pylab as py

path2='/home/jahndav/Dropbox/THz-Analysis/'    
#    samfiles=glob.glob(path2+'MarburgData/*_Lact1*')
samfiles=glob.glob('/home/jahndav/Dropbox/THz-Analysis/rehi/Sample_?.txt')
reffiles=glob.glob('/home/jahndav/Dropbox/THz-Analysis/rehi/Ref*.txt')

sam_td=TeraData.ImportMarburgData(samfiles)
ref_td=TeraData.ImportMarburgData(reffiles)

#the origginal fdData:
sam_fd=TeraData.FdData(sam_td)
ref_fd=TeraData.FdData(ref_td)
H=Terapy.HMeas(ref_fd,sam_fd)

#the windowed data
sam_td.setTDData(sam_td.getWindowedData(1e-12))
ref_td.setTDData(ref_td.getWindowedData(1e-12))

sam_fd_windowed=TeraData.FdData(sam_td)
ref_fd_windowed=TeraData.FdData(ref_td)

H_windowed=Terapy.HMeas(ref_fd_windowed,sam_fd_windowed)

#the zeropadded and windowed data
sam_td.zeroPaddData(sam_td.getLength()*2)
ref_td.zeroPaddData(ref_td.getLength()*2)

sam_fd_zpandwind=TeraData.FdData(sam_td)
ref_fd_zpandwind=TeraData.FdData(ref_td)

H_zpandwind=Terapy.HMeas(ref_fd_zpandwind,sam_fd_zpandwind)

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
