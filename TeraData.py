import pylab as py
from scipy.interpolate import interp1d

class THzTdData():

    def __init__(self,*kw,**kwargs):
        self._thzdata_raw=[]

        if not kwargs.has_key('existing'):
            #create instance by loading filenames
            filenames=kw[0]
            if len(kw)<2:
                params=[1,0,1,'.',1]
            else:
                params=kw[1]
            self.filename=filenames
            self.numberOfDataSets=len(self.filename)
            
     
            for fname in self.filename:
                self._thzdata_raw.append(self.importfile(fname,params))
            
            self.resetTDData()
        else:
            
            self.setTDData(kw[0])
            self._thzdata_raw=[self.tdData]            
            if len(kw)>1:
                self.filename=kw[1]
            else:
                self.filename=['None']
            
            self.numberOfDataSets=len(self.filename)                
            if len(kw)>2:
                self._thzdata_raw=kw[2]
    
    def calcTDData(self,tdDatas):
        #tdDatas should be an array of tdDataarrays
        meantdData=py.mean(tdDatas,axis=0)
        
        unc=self.calcunc(tdDatas)
        return py.column_stack((meantdData[:,:2],unc))       
    
    def getPreceedingNoise(self):
        precnoise=[]        
        for i in range(self.numberOfDataSets):
            #overthink this function!
            ratio=2
            ix_max=py.argmax(self._thzdata_raw[i][:,1])
            ix_min=py.argmin(self._thzdata_raw[i][:,1])
            earlier=min(ix_max,ix_min)
            #take only the first nts part            
            ix_start=int(earlier/ratio)
            precnoise=py.concatenate((precnoise, self._thzdata_raw[i][:ix_start,1]))
        return precnoise
    
    def calcunc(self,tdDatas):
        #tdDatas is a np array of tdData measurements
        if tdDatas.shape[0]==1:
            uncarray=py.zeros((len(tdDatas[0,:,0]),))
        else:
            uncarray=py.std(py.asarray(tdDatas),axis=0)      
        return uncarray
       
    def importfile(self,fname,params):
        # if even more sophisticated things are needed, just inherit THzTdData class
        #and override the importfile method
        
        if params[3]=='.':
            data=py.loadtxt(fname,usecols=(params[1],params[2]))
        elif params[3]==',':
            str2float=lambda val: float(val.replace(',','.'))
            data=py.loadtxt(fname,
            converters={params[1]:str2float,params[2]:str2float},
            usecols=(params[1],params[2]),skiprows=params[4])                
        data[:,0]*=params[0]
        
        if data[1,0]-data[0,0]<0:
            data=py.flipud(data)
     
        return data

    def _processRawData(self,tdDatas):
        #so far only offset removal, think of windowing etc...
        tempTDDatas=[]
        for tdData in tdDatas:
            t=self._removeLinearDrift(tdData)
            tempTDDatas.append(t)
        #before interpolating to a common time axis, this need to be a list of 
        #tdData arrays, since they might differ in length
        #first shift than interpolate 
        tempTDDatas,time_jitter=self._removeTimeShift(tempTDDatas)
        tempTDDatas=self._bringToCommonTimeAxis(tempTDDatas)
        return tempTDDatas
        
    def _bringToCommonTimeAxis(self,tdDatas):
        #What can happen: 
        #a) datalengthes are not equal due to missing datapoints 
        #   => no equal frequency bins
        #b) time positions might not be equal
            
        miss_points_max=10    
        #check for missing datapoints, allowing 
        #not more than miss_points_max points to miss 
        all_lengthes=[]
        for thisdata in tdDatas:
            all_lengthes.append(len(thisdata[:,0]))
        
        #do it always, just print a warning, if miss_points_max is exceeded
        if min(all_lengthes)!=max(all_lengthes):
            print "Datalength of suceeding measurements not consistent, try to fix"
            if max(all_lengthes)-min(all_lengthes)>miss_points_max:
                print "Warning: Data seems to be corrupted. \n" +\
                "The length of acquired data of repeated measurements differs by \n" + \
                    str(max(all_lengthes)-min(all_lengthes)) + ' datapoints'
        #interpolation does no harm, even if everything is consistent (no interpolation in this case)
        commonMIN=min(tdDatas[0][:,0])
        commonMAX=max(tdDatas[0][:,0])
        commonLENGTH=len(tdDatas[0][:,0])
        
        for i in range(1,self.numberOfDataSets):
            commonMIN=max(commonMIN,min(tdDatas[i][:,0]))
            commonMAX=min(commonMAX,max(tdDatas[i][:,0]))
            commonLENGTH=min(commonLENGTH,len(tdDatas[i][:,0]))
            
        for i in range(self.numberOfDataSets):
            tdDatas[i]=self.getInterData(tdDatas[i],commonLENGTH,commonMIN,commonMAX)
            tdDatas[i][:,0]-=commonMIN
        
        return py.asarray(tdDatas)
    
    def _removeTimeShift(self,tdDatas):
        #not sure if needed, maybe we want to correct the rawdata by shifting the maxima on top of each other
        #the indices of the maxima        
        peak_pos=[]
        for tdData in tdDatas:
            time_max_raw=tdData[py.argmax(tdData[:,1]),0]
            thisPeakData=self.getShorterData(tdData,time_max_raw-0.5e-12,time_max_raw+0.5e-12)
            thisPeakData=self.getInterData(thisPeakData,len(thisPeakData[:,0])*20,thisPeakData[0,0],thisPeakData[-1,0],'cubic')
            peak_pos.append(thisPeakData[py.argmax(thisPeakData[:,1]),0])
            
        for i in range(len(tdDatas)):
            tdDatas[i][:,0]-=peak_pos[i]
        return tdDatas,py.asarray(peak_pos).std()
  
    def _removeLinearDrift(self,tdData):
        #not yet implemented
        signal=tdData[:,1]
        tdData[:,1]=py.detrend(signal,key='linear')
        return tdData
    
    def getInterData(self,tdData,desiredLength,mint,maxt,tkind='linear'):
        #use cubic interpolation only, if you know, that the data is not too no
        intpdata=interp1d(tdData[:,0],tdData[:,1:],axis=0,kind=tkind)
        timeaxis=py.linspace(mint,maxt,desiredLength)
        longerData=intpdata(timeaxis)
        return py.column_stack((timeaxis,longerData))
#            return py.transpose(py.asarray([timeaxis,longerData]))
    
    def getLength(self):
        return len(self.tdData[:,0])
        
    def getfilename(self):
        return self.filename
    
    def getTimeWindowLength(self):
        peak=self.getPeakPosition()
        return abs(self.tdData[-1,0]-peak)
    
    def getPeakPosition(self):
        return self.tdData[py.argmax(self.tdData[:,1]),0]
    
    def zeroPaddData(self,desiredLength,paddmode='zero'):    
        if desiredLength>0:
            desiredLength=int(desiredLength)
            if paddmode=='gaussian':
                paddvec=py.normal(0,py.std(self.getPreceedingNoise())*0.05,desiredLength)
                                
            else:
                paddvec=py.ones((desiredLength,))*py.mean(self.tdData[-20:,1])
            
            uncpadd=py.ones((desiredLength,))*py.mean(self.tdData[-20:,2])
            longunc=py.append(self.tdData[:,2],uncpadd)
            longdata=py.append(self.tdData[:,1],paddvec)
            longtime=py.append(self.tdData[:,0],py.linspace(self.tdData[-1,0],self.tdData[-1,0]+desiredLength*self.dt,desiredLength))
            
            self.setTDData(py.column_stack((longtime,longdata,longunc)))
            
    def getFirstPuls(self,before,after):
        tcenter=self.getPeakPosition()
        #assume width of the pulse to be not more than 10ps!
        tmin=tcenter-before
        tmax=tcenter+after
        return self.getShorterData(self.tdData,tmin,tmax)
    
    def getShorterData(self,tdData,tmin,tmax):
        ix=py.all([tdData[:,0]>=tmin,tdData[:,0]<tmax],axis=0)
        return tdData[ix,:]
    
    
    def setTDData(self,tdData):
        #if tdData array is changed outside, use this function, to guarantee data integrity

        self.tdData=tdData
        self.dt=abs(py.mean(tdData[10:20,0]-tdData[9:19,0]))       
#        self.dt=(tdData[-1,0]-tdData[0,0])/len(tdData[:,0])
        self.num_points=len(tdData[:,0])
    
    def getWindowedData(self,windowlength):
        N=int(windowlength/self.dt)
        w=py.blackman(N*2)
        w=py.asarray(py.hstack((w[0:N],py.ones((self.getLength()-N*2),),w[N:])))
        windowedData=self.tdData
        windowedData[:,1]*=w        
        return windowedData
        
    def resetTDData(self):
        tempdata=self.calcTDData(self._processRawData(self._thzdata_raw))
        self.setTDData(tempdata)
    
    def doTdPlot(self,name='TD-Plot'):
        py.figure(name)
        py.plot(self.tdData[:,0]*1e12,self.tdData[:,1])      
        py.xlabel('Time in ps')
        py.ylabel('Amplitude, arb. units')
        
    def doPlotWithunc(self,no_std=2):
        self.doTdPlot('TD-UNC-Plot')
        py.plot(self.tdData[:,0]*1e12,self.tdData[:,1]+no_std*self.tdData[:,2],'g--')
        py.plot(self.tdData[:,0]*1e12,self.tdData[:,1]-no_std*self.tdData[:,2],'g--')

        
class ImportMarburgData(THzTdData):

    def __init__(self,filename):
        Params=self._filePreferences()        
        THzTdData.__init__(self,filename,Params)
        
    def _filePreferences(self):
        time_factor=1
        colon_time=0
        colon_data=1
        decimal_sep=','
        skiprows=0
        return [time_factor,colon_time,colon_data,decimal_sep,skiprows]

class ImportInrimData(THzTdData):
  
    def __init__(self,filename):
        Params=self._filePreferences()        
        THzTdData.__init__(self,filename,Params)
  
    def _filePreferences(self):        
        time_factor=1  #0
        colon_time=2    #1
        colon_data=3    #2
        decimal_sep='.' #3
        skiprows=0      #4
        return [time_factor,colon_time,colon_data,decimal_sep,skiprows]

class FdData():
    #hard frequency axis
    FMIN=0
    FMAX=5e12    
    
    def __init__(self,tdData,fbins=-1,fbnds=[FMIN,FMAX]):
        #FdData is always attached to a TdData Object
        
        #the original _tdData object should be private and not accessed from outside
        #if you want to know it use FdData.getassTDData()
        self._tdData=tdData
        
        #self._fdData_raw holds the raw _fdData should also not be altered from outside
        #fdData have always the following structure #col1=freq,col2=complexfft, col3=abs(col2), col4=ph(col2)    
        self.setFDData(self.calculatefdData(self._tdData))       
        self.resetfdData(fbins,fbnds)


    def getInterpolatedFDData(self,newfbins,strmode='linear'):
        oldfreqs=self.getfreqs()
        newfreqs=py.arange(min(oldfreqs),max(oldfreqs),newfbins)
        interfdData=interp1d(oldfreqs,self.fdData[:,1:],strmode,axis=0)
        
        return py.column_stack((newfreqs,interfdData(newfreqs)))
     
      
    def setTDData(self,tdData):
        #Neccessary? couldn't we just create a new object with init?
        fbins_old=self.getfbins()
        minf_old=min(self.getfreqs())
        maxf_old=max(self.getfreqs())
        self._tdData=tdData        
        self.setFDData(self.calculatefdData(self._tdData))
        
        self.resetfdData(fbins_old,[minf_old,maxf_old])

    def setFDData(self,fdData):
        self.fdData=fdData
        
    def getassTDData(self):
        return self._tdData
    
    def resetfdData(self,fbins=-1,fbnds=[FMIN,FMAX]):

        if fbins>0 and self.getfbins()>fbins:
            self.zeroPadd(fbins)
            
        self.setFDData(self.cropData(self.fdData,fbnds[0],fbnds[1]))

    def zeroPadd(self,fbins):
        spac=1/self._tdData.dt/fbins
        actlen=self._tdData.getLength()
        nozeros=py.ceil(spac-actlen)
        self._tdData.zeroPaddData(nozeros)
        #leave the old bnds in place
        bnds=[min(self.getfreqs()),max(self.getfreqs())]
        zpd=self.calculatefdData(self._tdData)
        self.setFDData(self.cropData(zpd,bnds[0],bnds[1]))
        
    def cropData(self,data,startfreq=FMIN,endfreq=FMAX):
        ix=py.all([data[:,0]>=startfreq,data[:,0]<=endfreq],axis=0)
        return data[ix,:]
    
    def getfbins(self):
        return abs(self.fdData[1,0]-self.fdData[0,0])
    
    def getfreqs(self):
        return self.fdData[:,0].real
        
    def getmaxfreq(self):
        #take care of constant offset
        i=1     
        fmax=30e9
        while fmax<100e9:
            fmax=self.fdData[py.argmax(self.fdData[i:,1]),0]
            i+=1       
        
        return fmax
        
    def getUnwrappedPhase(self):        
        return self.fdData[:,3]
    
    def calculatefdData(self,tdData):
        #this could also be avoided, since it doesn't change in most of the cases
        unc=self.calculateSTDunc()
        #no need to copy it before (access member variable is ugly but shorter)
    
        fd=py.fft(tdData.tdData[:,1])
        fdabs=abs(fd)
        
        fdph=abs(py.unwrap(py.angle(fd)))
        dfreq=py.fftfreq(tdData.num_points,tdData.dt)                
        fdph=self.removePhaseOffset(py.column_stack((dfreq,fdph)))
        t=py.column_stack((dfreq,fd,fdabs,fdph))
        t=self.cropData(t,0,unc[-1,0])
        intpunc=interp1d(unc[:,0],unc[:,1:],axis=0)        
        unc=intpunc(t[:,0])
        return py.column_stack((t[:,:4],unc))
    
    def removePhaseOffset(self,ph):
        #cut data to reasonable range:
        lower=200e9
        upper=1e12        
        ph_c=self.cropData(ph,lower,upper)
        p=py.polyfit(ph_c[:,0],ph_c[:,1],1)
        return ph[:,1]-p[1]

    def calculateSTDunc(self):
        #return asarray _tdData
        #first calculate white noise contribution
        #should be independent of underlying data change!
        precNoise=py.fft(self._tdData.getPreceedingNoise())
        
        #make sure, that no interpolated or zero padded data is used for calculating the 
        #freqeuncy uncertainty, not nice style!
        commonTdData=self._tdData._bringToCommonTimeAxis(self._tdData._thzdata_raw)
        noise_real=py.std(precNoise.real)*py.ones(commonTdData[0][:,0].shape)
        noise_imag=py.std(precNoise.imag)*py.ones(commonTdData[0][:,0].shape)
        
        
        #second calculate variation between measurements
        if self._tdData.numberOfDataSets<=1:
            repeat_noise_imag=0
            repeat_noise_real=0
        else:
            a=py.fft(commonTdData[:,:,1],axis=1)
            repeat_noise_real=py.std(a.real,axis=0)
            repeat_noise_imag=py.std(a.imag,axis=0)
       
        dfreq=py.fftfreq(len(commonTdData[0][:,0]),commonTdData[0][5,0]-commonTdData[0][4,0])
        t=py.column_stack((dfreq,py.sqrt(noise_real**2+repeat_noise_real**2),py.sqrt(noise_imag**2+repeat_noise_imag**2)))
        return self.cropData(t)
    
    def getSNR(self):
        noiselevel=py.mean(abs(py.fft(self._tdData.getPreceedingNoise())))
        maxfreq=self.getmaxfreq()
        signalmax=self.cropData(self.fdData,maxfreq-4*self.getfbins(),maxfreq+4*self.getfbins())
        return 20*py.log10(py.mean(signalmax[:,2])/noiselevel)        
        
    def getBandwidth(self):
        #this function should return the lowest trustable and highest trustable frequency, 
        # along with the resulting bandwidth
        SNR=self.getSNR()
        absdata=-20*py.log10(self.fdData[:,2]/max(self.fdData[:,2]))
        ix=SNR-10>absdata #dangerous, due to sidelobes, there might be some high freq component!
        tfr=self.fdData[ix,0]
        return min(tfr),max(tfr)

    def doPlot(self):

        py.figure('FD-ABS-Plot')
        py.plot(self.fdData[:,0].real/1e9,20*py.log10(abs(self.fdData[:,2])))
        py.xlabel('Frequency in GHz')
        py.ylabel('Amplitude, arb. units')
        
        py.figure('FD-PHASE-Plot')
        py.plot(self.fdData[:,0].real/1e9,self.fdData[:,3].real)                

if __name__=='__main__':
    import glob
    
    samfiles=glob.glob('/home/jahndav/Dropbox/THz-Analysis/rehi/Sample*')
 
    myTDData=THzTdData(samfiles)
#    myTDData.doTdPlot()
    
#    myTDData.doPlotWithunc()

    
#    myFDData=FdData(myTDData)
#    myFDData.doPlot()
#    print myFDData.getSNR()    
#    myFDData.doPlot()    
#    
#    myFDData.setFDData(myFDData.getInterpolatedFDData(5e9))
#    myFDData.doPlot()
    
    
    
    #myFDData.setFDData(myFDData.cropData(myFDData.fdData,150e9,4e12))
#    myFDData.doPlot()    
