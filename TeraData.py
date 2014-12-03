import pylab as py
from scipy.interpolate import interp1d
import scipy.signal as signal
import sys
from uncertainties import unumpy

class THzTdData():

    def __init__(self,*kw,**kwargs):
        self._thzdata_raw=[]

        if not 'existing' in kwargs:
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
    
    def getPreceedingNoise(self,timePreceedingSignal=-1):
        precnoise=[]
        
        #check for reasonable input
        #peak position in raw data is always at origin!
        nearestdistancetopeak=2.5e-12
        if min(self._thzdata_raw[0][:,0])+timePreceedingSignal>-nearestdistancetopeak:
            timePreceedingSignal=-min(self._thzdata_raw[0][:,0])-nearestdistancetopeak
            print("Time Interval of preceeding noise to long, reset done")
            #dont use user input if higher than peak
       #really neccessary for each dataset, or would it suffice to take the first?
        for i in range(self.numberOfDataSets):       
            starttime=min(self._thzdata_raw[i][:,0])
            if timePreceedingSignal==-1:
              #determine length automatically            

                ratio=2
                ix_max=py.argmax(self._thzdata_raw[i][:,1])
                ix_min=py.argmin(self._thzdata_raw[i][:,1])
                earlier=min(ix_max,ix_min)
                #take only the first nts part            
                endtime=self._thzdata_raw[i][int(earlier/ratio),0]
            else:
                    
                endtime=starttime+timePreceedingSignal
            
            noise=self.getShorterData(self._thzdata_raw[i],starttime,endtime)
            noise=py.detrend(noise[:,1])
#            noise=noise[:,1]
            precnoise=py.concatenate((precnoise,noise))
        return precnoise
    
    def calcunc(self,tdDatas):
        #tdDatas is a np array of tdData measurements
        if tdDatas.shape[0]==1:
            uncarray=py.zeros((len(tdDatas[0,:,0]),))
        else:
            uncarray=py.std(py.asarray(tdDatas[:,:,1]),axis=0)    
        return uncarray
       
    def importfile(self,fname,params):
        # if even more sophisticated things are needed, just inherit THzTdData class
        #and override the importfile method
        try:
            if params[3]=='.':
                data=py.loadtxt(fname,usecols=(params[1],params[2]))
            elif params[3]==',':
                str2float=lambda val: float(val.replace(',','.'))
                data=py.loadtxt(fname,
                converters={params[1]:str2float,params[2]:str2float},
                usecols=(params[1],params[2]),skiprows=params[4])                
        except IOError:
            print "File " + fname + " could not be loaded"
            sys.exit()

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
            print("Datalength of suceeding measurements not consistent, try to fix")
            if max(all_lengthes)-min(all_lengthes)>miss_points_max:
                print("Warning: Data seems to be corrupted. \n" +\
                "The length of acquired data of repeated measurements differs by \n" + \
                    str(max(all_lengthes)-min(all_lengthes)) + ' datapoints')
        
        #interpolation does no harm, even if everything is consistent (no interpolation in this case)
        commonMIN=max([thistdData[:,0].min() for thistdData in tdDatas])
        commonMAX=min([thistdData[:,0].max() for thistdData in tdDatas])
        commonLENGTH=min([thistdData[:,0].shape[0] for thistdData in tdDatas])
            
        for i in range(self.numberOfDataSets):
            tdDatas[i]=self.getInterData(tdDatas[i],commonLENGTH,commonMIN,commonMAX)
        
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
        
        peak_pos=py.asarray(peak_pos)
        mp=py.mean(peak_pos)
        for i in range(len(tdDatas)):
            tdDatas[i][:,0]-=(peak_pos[i]-mp)
        return tdDatas,py.std(peak_pos)
  
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
    
    def getDR(self):
        Emax=abs(self.tdData[:,1])
        noise=py.sqrt(py.mean(self.getPreceedingNoise()**2))
        return Emax/noise
    
    def getSNR(self):
        return abs(self.tdData[:,1])/self.tdData[:,2]
    
    def getTimeWindowLength(self):
        peak=self.getPeakPosition()
        return abs(self.tdData[-1,0]-peak)
    
    def getPeakPosition(self):
        return self.tdData[py.argmax(self.tdData[:,1]),0]
    
    def zeroPaddData(self,desiredLength,paddmode='zero',where='end'):    
        if desiredLength<0:
            return 0
            
        desiredLength=int(desiredLength)
        
        if paddmode=='gaussian':
            paddvec=py.normal(0,py.std(self.getPreceedingNoise())*0.05,desiredLength)
                                
        else:
            paddvec=py.ones((desiredLength,))*py.mean(self.tdData[-20:,1])

        if where=='end':
            uncpadd=py.ones((desiredLength,))*py.mean(self.tdData[-20:,2])
            longunc=py.append(self.tdData[:,2],uncpadd)
            longdata=py.append(self.tdData[:,1],paddvec)
            longtime=py.append(self.tdData[:,0],py.linspace(self.tdData[-1,0],self.tdData[-1,0]+desiredLength*self.dt,desiredLength))
        else:
            uncpadd=py.ones((desiredLength,))*py.mean(self.tdData[:20,2])
            longunc=py.append(uncpadd,self.tdData[:,2])
            longdata=py.append(paddvec,self.tdData[:,1])
            timepadd=py.linspace(self.tdData[0,0]-(desiredLength+1)*self.dt,self.tdData[0,0],desiredLength)
            longtime=py.append(timepadd,self.tdData[:,0])
        self.setTDData(py.column_stack((longtime,longdata,longunc)))
            
    def getFirstPuls(self,after):
        tcenter=self.getPeakPosition()
        #assume width of the pulse to be not more than 10ps!, further crop always from
        #the beginning of the time data, in order to avoid phase problems
        
        tmax=tcenter+after
        return self.getShorterData(self.tdData,self.tdData[0,0],tmax)
    
    def getShorterData(self,tdData,tmin,tmax):
        ix=py.all([tdData[:,0]>=tmin,tdData[:,0]<tmax],axis=0)
        return tdData[ix,:]
    
    
    def setTDData(self,tdData):
        #if tdData array is changed outside, use this function, to guarantee data integrity

        self.tdData=tdData
        self.dt=abs(py.mean(tdData[10:20,0]-tdData[9:19,0]))       
#        self.dt=(tdData[-1,0]-tdData[0,0])/len(tdData[:,0])
        self.num_points=len(tdData[:,0])
    
    def getWindowedData(self,windowlength_time):
        N=int(windowlength_time/self.dt)
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
        time_factor=1    #0
        colon_time=2   #1
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
        self.maxDR=max(self.getDR())

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
            
        self.setFDData(self.getcroppedData(self.fdData,fbnds[0],fbnds[1]))

    def zeroPadd(self,fbins):
        spac=1/self._tdData.dt/fbins
        actlen=self._tdData.getLength()
        nozeros=py.ceil(spac-actlen)
        self._tdData.zeroPaddData(nozeros)
        #leave the old bnds in place
        bnds=[min(self.getfreqs()),max(self.getfreqs())]
        zpd=self.calculatefdData(self._tdData)
        self.setFDData(self.getcroppedData(zpd,bnds[0],bnds[1]))
        
    def getcroppedData(self,data,startfreq=FMIN,endfreq=FMAX):
        ix=py.all([data[:,0]>=startfreq,data[:,0]<=endfreq],axis=0)
        return data[ix,:]
    
    def getfbins(self):
        return abs(self.fdData[1,0]-self.fdData[0,0])
    
    def getfreqs(self):
        return self.fdData[:,0].real

    def getfreqsGHz(self):
        return self.getfreqs()*1e-9         
    
    def getmaxfreq(self):
        #take care of constant offset
        cutted=self.getcroppedData(self.fdData,150e9,FdData.FMAX)
        fmax=cutted[py.argmax(cutted[:,2]),0]
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
        t=self.getcroppedData(t,0,unc[-1,0])
        intpunc=interp1d(unc[:,0],unc[:,1:],axis=0)        
        unc=intpunc(t[:,0])
        return py.column_stack((t[:,:4],unc))
    
    def removePhaseOffset(self,ph):
        #cut data to reasonable range:
        lower=200e9
        upper=1e12        
        ph_c=self.getcroppedData(ph,lower,upper)
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
        dfreq=py.fftfreq(len(commonTdData[0][:,0]),commonTdData[0][5,0]-commonTdData[0][4,0])
            
        noise_real=py.std(precNoise.real)*py.ones(commonTdData[0][:,0].shape)
        noise_imag=py.std(precNoise.imag)*py.ones(commonTdData[0][:,0].shape)
        noise_abs=py.std(abs(precNoise))*py.ones(commonTdData[0][:,0].shape)
        noise_ph=py.std(py.angle(precNoise))*py.ones(commonTdData[0][:,0].shape)
        
        
        #second calculate variation between measurements
        if self._tdData.numberOfDataSets<=1:
            repeat_noise_imag=0
            repeat_noise_real=0
            repeat_noise_abs=0
            repeat_noise_ph=0
            
        else:
            a=py.fft(commonTdData[:,:,1],axis=1)
            b=abs(a)
            c=[]
            for i in range(self._tdData.numberOfDataSets):
                c.append(self.removePhaseOffset(py.column_stack((dfreq,py.unwrap(py.angle(a[i,:]))))))
            c=py.asarray(c)
            
            repeat_noise_abs=py.std(b,axis=0)
            repeat_noise_ph=py.std(c,axis=0)
            repeat_noise_real=py.std(a.real,axis=0)
            repeat_noise_imag=py.std(a.imag,axis=0)
       
        t=py.column_stack((dfreq,py.sqrt(noise_real**2+repeat_noise_real**2),py.sqrt(noise_imag**2+repeat_noise_imag**2)))
        t=py.column_stack((t,py.sqrt(noise_abs**2+repeat_noise_abs**2),py.sqrt(noise_ph**2+repeat_noise_ph**2)))        
        return self.getcroppedData(t)    
    
    def getDR(self):
        
        noiselevel=py.sqrt(py.mean(abs(py.fft(self._tdData.getPreceedingNoise()))**2))
        window_size=5
        window=py.ones(int(window_size))/float(window_size)
        hlog=py.convolve(20*py.log10(self.fdData[:,2].real), window, 'valid')
        one=py.ones((2,))
        hlog=py.concatenate((hlog[0]*one,hlog,hlog[-1]*one))
        return hlog-20*py.log10(noiselevel)         
    
    def getSNR(self):
        return self.fdData[:,2]/self.fdData[:,6]
        
    def getBandwidth(self,dbDistancetoNoise=15):
        #this function should return the lowest trustable and highest trustable frequency, 
        # along with the resulting bandwidth
        
        absdata=-20*py.log10(self.fdData[:,2]/max(self.fdData[:,2]))
        ix=self.maxDR-dbDistancetoNoise>absdata #dangerous, due to sidelobes, there might be some high freq component!
        tfr=self.fdData[ix,0]

        return min(tfr),max(tfr)
        
    def findAbsorptionLines(self):
        #no interpolation here, just return the measured data peak index
        tresh=6 #db treshhold to detect a line
        #filter first
    
        #bring to logarithmic scale
        movav=self.getmovingAveragedData()        
        hlog=20*py.log10(movav[:,2])
        #smooth data
          #remove etalon minima
        ixlines_prob=signal.argrelmin(hlog)[0]
        ixlines=[]
        ixmax=signal.argrelmax(hlog)[0]        
        #check depth of minima compared two the two adjacent maxima
        
        if ixmax[0]>ixlines_prob[0]:
            #true if starts with minimum, treat this case separatedly
            if hlog[ixmax[0]]-hlog[ixlines_prob[0]]>tresh:
                ixlines.append(ixlines_prob[0])
                #and remove first minimum from ixlines_prob
            ixlines_prob=py.delete(ixlines_prob,0)
                
        if ixmax[-1]<ixlines_prob[-1]: #treat this case separatedly
            if hlog[ixmax[-1]]-hlog[ixlines_prob[-1]]>tresh:
                ixlines.append(ixlines_prob[0])
            ixlines_prob=py.delete(ixlines_prob,-1)

        #now the remaining elements of ixlines, ixmax should be in the order of max min max ... min max         
        leftdist=hlog[ixmax[:-1]]-hlog[ixlines_prob]
        rightdist=hlog[ixmax[1:]]-hlog[ixlines_prob]
        
        for i in range(len(ixlines_prob)):
            #check if distance is higher than treshhold
            if (leftdist[i]>tresh or rightdist[i]>tresh):
                ixlines.append(ixlines_prob[i])
        
        #remove inappropriate lines (by distance so far, at least 5 datapoints!
        ixlines=py.asarray(ixlines)
        ixlines=py.concatenate(([ixlines[0]],ixlines[py.diff(ixlines)>5]))
        return ixlines
    
    def findAbsorptionPeaks_TESTING(self):
#        movav=self.getmovingAveragedData()        
        hlog=-20*py.log10(self.fdData[:,2])
        etalon=myFDData.getEtalonSpacing()
        Ns=int(etalon/myFDData.getfbins())
        Ns=py.arange(max(1,Ns-10),Ns+10,1)        
        peaks=signal.find_peaks_cwt((hlog),Ns)
        return peaks
        
    def getEtalonSpacing(self):
        #how to find a stable range! ? 
#        etalonData=self.getFilteredData(20e9,3)
        bw=self.getBandwidth()
        rdata=self.getcroppedData(self.fdData,max(bw[0],250e9),min(bw[1],3e12))  
        #get etalon frequencies:
        #need to interpolate data!
        oldfreqs=rdata[:,0].real       
        intpdata=interp1d(oldfreqs,rdata[:,2],'cubic')
        
        fnew=py.arange(min(oldfreqs),max(oldfreqs),0.1e9)
        absnew=intpdata(fnew)
        ixmaxima=signal.argrelmax(absnew)[0]
        ixminima=signal.argrelmin(absnew)[0]
        
        fmaxima=py.mean(py.diff(fnew[ixmaxima]))
        fminima=py.mean(py.diff(fnew[ixminima]))
        
        df=(fmaxima+fminima)*0.5 #the etalon frequencies
        print(str(df/1e9) + " GHz estimated etalon frequency")
        return df

    def getFilteredData(self,windowlen=100e9,maxorder=3):
        #windowlength should be smaller the etalon frequency, in order to emphasize them
        #for the length algorithm
        #it should be of the order of typical absorption line widthes for their analysis
        #remove etalon minima
        
        #what happens if the savgol filter is applied to real and imaginary part of the fft?!        
        #maybe nice to check?!
        fbins=self.getfbins()
        #100e9 should be replaced by the Etalon frequency 
        N_min=int(windowlen/fbins)
        order=min(N_min-1,maxorder)
        
        absdata=signal.savgol_filter(self.fdData[:,2],N_min-N_min%2+1,order)
        phdata=signal.savgol_filter(self.fdData[:,3],N_min-N_min%2+1,order)
       
        return py.column_stack((self.fdData[:,0],self.fdData[:,1],absdata,phdata,self.fdData[:,4:]))
    
    def getmovingAveragedData(self,window_size_GHz=-1):
        #so far unelegant way of convolving the columns one by one
        if window_size_GHz<0.5e9:
            window_size=int(self.getEtalonSpacing()/self.getfbins())
        else:
            window_size=int(window_size_GHz/self.getfbins())
              
        window_size+=window_size%2+1
        window=py.ones(int(window_size))/float(window_size)

        dataabs=py.convolve(self.fdData[:,2], window, 'valid')
        dataph=py.convolve(self.fdData[:,2], window, 'valid')
        one=py.ones((window_size-1)/2,)
        dataabs=py.concatenate((dataabs[0]*one,dataabs,dataabs[-1]*one))
        dataph=py.concatenate((dataph[0]*one,dataph,dataph[-1]*one))
        return py.column_stack((self.fdData[:,0],self.fdData[:,1],dataabs,dataph,self.fdData[:,4:]))
    
    def doPlot(self):

        py.figure('FD-ABS-Plot')
        py.plot(self.getfreqsGHz(),20*py.log10(abs(self.fdData[:,2])))
        py.xlabel('Frequency in GHz')
        py.ylabel('Amplitude, arb. units')
        
        py.figure('FD-PHASE-Plot')
        py.plot(self.getfreqsGHz(),self.fdData[:,3].real)   
             
    def doPlotWithUnc(self):
        f=1

        uabs=unumpy.uarray(self.fdData[:,2].real,self.fdData[:,6].real)
        uabs=20*unumpy.log10(uabs)
        u_a=unumpy.nominal_values(uabs)
        u_s=unumpy.std_devs(uabs)
        py.figure('FD-ABS-UNC-Plot')
        py.plot(self.getfreqsGHz(),u_a)
        py.plot(self.getfreqsGHz(),u_a+u_s,'--',self.getfreqsGHz(),u_a-u_s,'--')

        py.figure('FD-PH-UNC-Plot')
        py.plot(self.getfreqsGHz(),self.fdData[:,3])
        py.plot(self.getfreqsGHz(),self.fdData[:,3]+f*self.fdData[:,7],'--',self.getfreqsGHz(),self.fdData[:,3]-f*self.fdData[:,7],'--')

if __name__=='__main__':
    import glob
    path2='/home/jahndav/Dropbox/THz-Analysis/'    
#    samfiles=glob.glob(path2+'MarburgData/*_Lact1*')
    samfiles=glob.glob('/home/jahndav/Dropbox/THz-Analysis/rehi/Sam*')
 
    myTDData=ImportMarburgData(samfiles)

    myFDData=FdData(myTDData)
    myFDData.doPlotWithUnc()
#    py.plot(myFDData.getfreqsGHz(),myFDData.fdData[:,2])
#    peaks=myFDData.findAbsorptionLines()
#    peaksS=myFDData.findAbsorptionPeaks_TESTING()
#    py.plot(myFDData.getfreqsGHz(),myFDData.fdData[:,2])
#    py.plot(myFDData.getfreqsGHz()[peaks],myFDData.fdData[peaks,2],'x')
#    py.plot(myFDData.getfreqsGHz()[peaksS],myFDData.fdData[peaksS,2],'*')
    
#    myFDData.setFDData(myFDData.getcroppedData(myFDData.fdData,200e9,2200e9))
#    etalon=myFDData.getEtalonSpacing()
#    Ns=int(etalon/myFDData.getfbins())
#    Ns=py.arange(max(1,Ns-10),Ns+10,1)
#    hlog=20*py.log10(myFDData.fdData[:,2])
#    peaks=signal.find_peaks_cwt((hlog**2),Ns)
#    f=myFDData.getfreqsGHz()
#        
#    py.plot(f,hlog)
#    py.plot(f[peaks],hlog[peaks],'*')