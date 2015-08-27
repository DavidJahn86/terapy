import pylab as py
from scipy.interpolate import interp1d
import scipy.signal as signal
import sys
from uncertainties import unumpy

class THzTdData():

    def __init__(self,*kw,**kwargs):
        '''
        This is a general TimeDomain-Class,
        Instances can be initialized by
        1) create TDData Object from measurement files, without passing 'existing'
        2) create TDData Object from existing TDData Object by passing 'existing'
        '''
        #_thzdata_raw is an array of measurments,dim=3, for each measurement
        # an array t,X,Y is stored
        self._thzdata_raw=[]
        #two different constructors: 
        if not 'existing' in kwargs:
            #create instance by loading from the filenames provided as first argument
            filenames=kw[0]
            if len(kw)<2:
                #if no parameter set is set, use the following dictionary
                params={'time_factor':1,
                        'time_col':0,
                        'X_col':1,
                        'Y_col':2,
                        'dec_sep':'.',
                        'skiprows':0}
            else:
                #else use the passed dictionary with the fileformat
                params=kw[1]
  
            #filenames  
            self.filename=filenames
            #number of measurements
            self.numberOfDataSets=len(self.filename)
            
            #import files
            for fname in self.filename:
                self._thzdata_raw.append(self.importfile(fname,params))
            
            #calculate TDData array along with uncertainties and so on
            self.resetTDData()
        else:
            #the second constructor (with keyword), provides at least a 
            #TDData array
            self.setTDData(kw[0])
            self._thzdata_raw=[self.tdData]
            
            #set filename (only to prevent conflicts)
            if len(kw)>1:
                self.filename=kw[1]
            else:
                self.filename=['None']
            
            self.numberOfDataSets=len(self.filename) 
            #set rawdata, if also passed               
            if len(kw)>2:
                self._thzdata_raw=kw[2]
    
    def calcTDData(self,tdDatas):
        #tdDatas is a a 3d array of measurements, along with their uncertainties
        #meantdData is the weighted sum of the different measurements
        #meantdData,sumofweights=py.average(tdDatas[:,:,1:3],axis=0,weights=1.0/tdDatas[:,:,3:]**2,returned=True)
        meantdData=py.average(tdDatas[:,:,1:3],axis=0)
        #use error propagation formula
        noise=py.sqrt(py.mean(self.getAllPrecNoise()[0]**2))
        if tdDatas.shape[0]==1:
            rep = py.zeros((len(tdDatas[0,:,0]),2))
        else:
            rep = py.std(tdDatas[:,:,1:3],axis=0, ddof=1)/py.sqrt(self.numberOfDataSets)
        unc = py.sqrt(rep**2+noise**2)
        #unc=py.sqrt(1.0/sumofweights)
        #time axis are all equal
        return py.column_stack((tdDatas[0][:,0],meantdData,unc))       
    
    def calcunc(self,tdDatas):
        #not used anymore, older version, should we remove it???
         #tdDatas is a np array of tdData measurements
        if tdDatas.shape[0]==1:
            repeatability=py.zeros((len(tdDatas[0,:,0]),2))
        else:
            repeatability=py.std(py.asarray(tdDatas[:,:,1:3]),axis=0, ddof = 1)/py.sqrt(self.numberOfDataSets)
        #this line is wrong
        elnNoise=tdDatas[0,:,3:]
        uncarray = py.sqrt(repeatability**2 + elnNoise**2)
        
        return uncarray
    
    def doTdPlot(self,name='TD-Plot'):
        #plots the X and Y channel of the processed Data
        py.figure(name)
        py.plot(self.getTimesPs(),self.getEX())      
        py.plot(self.getTimesPs(),self.getEY())      
        py.xlabel('Time in ps')
        py.ylabel('Amplitude, arb. units')
        
    def doPlotWithunc(self,no_std=2):
        #plots the X and the Y channel of the processed measurements, along
        #with its uncertainties
        self.doTdPlot('TD-UNC-Plot')
        py.plot(self.getTimesPs(),self.getEX()+no_std*self.getUncEX(),'g--')
        py.plot(self.getTimesPs(),self.getEX()-no_std*self.getUncEX(),'g--')
        py.plot(self.getTimesPs(),self.getEY()+no_std*self.getUncEY(),'g--')
        py.plot(self.getTimesPs(),self.getEY()-no_std*self.getUncEY(),'g--')    
    
    def getAllPrecNoise(self,timePreceedingSignal=-1):
        #returns the concatenated Preceeding Noise
        precNoise=py.array([])
        for tdData in self._thzdata_raw:
            tN=self.getPreceedingNoise(tdData,timePreceedingSignal)
            if precNoise.shape[0]==0:
                precNoise=tN
            else:
                precNoise=py.vstack((precNoise,tN))
        return precNoise
    
    def getDR(self):
        #returns the dynamic range
        Emax=abs(self.getEX())
        noise=py.sqrt(py.mean(self.getAllPrecNoise()[0]**2))
        return Emax/noise
    
    def getelnNoise(self,tdData):
        #returns the uncertainty due to electronic noise
        
        #signal preceeding the pulse (X and Y channel)  
        precNoise=self.getPreceedingNoise(tdData)
        #is this normalization really correct?!
        elnNoise = py.std(precNoise, ddof = 1,axis=0)/py.sqrt(precNoise.shape[0])
        return elnNoise
       
    def getEX(self):
        #returns the mean X Channel
        return self.tdData[:,1]
        
    def getEY(self):
        #returns the mean Y Channel        
        return self.tdData[:,2]

    def getfilename(self):
        #returns the list of filenames, from which the data was loaded
        return self.filename
        
    def getFirstPuls(self,after=5e-12):
        #returns only the first pulse         
        
        #finds the main peak
        tcenter=self.getPeakPosition()
        #without argument, assume width of the pulse to be not more than 10ps!
        #further crop always from the beginning of the time data, 
        #in order to avoid phase problems
        
        tmax=tcenter+after
        return self.getShorterData(self.tdData,self.getTimes()[0],tmax)
        
    def getInterData(self,tdData,desiredLength,mint,maxt,tkind='linear'):
        #use cubic interpolation only, if you know, that the data is not too no
        intpdata=interp1d(tdData[:,0],tdData[:,1:],axis=0,kind=tkind)
        timeaxis=py.linspace(mint,maxt,desiredLength)
        longerData=intpdata(timeaxis)
        return py.column_stack((timeaxis,longerData))

    def getLength(self):
        #returns the length of the tdData array
        return len(self.getTimes())

    def getPeakPosition(self):
        #gives the time, at which the signal is maximal
        return self.tdData[py.argmax(self.getEX()),0]
    
    def getPeakWidth(self):
        #use as definition: 10 % over and under maxima around peak position
        position=self.getPeakPosition()
        E=self.getEX()
        Ebigger=abs(E)>0.1*abs(py.amax(E))
        times=self.tdData[Ebigger,0]
        times=times[times>position-3e-12]
        times=times[times<position+3e-12]        
        #make maximal pulse width smaller than 10 ps
        
        return py.amax(times)-py.amin(times)
       
       
    def getPreceedingNoise(self,tdData,timePreceedingSignal=-1):
        #retruns the preceeding noise in X and Y channel of tdData
        #it uses timePreceedingSignal to evaluate the noise, if not given,
        #it autodetects a "good" time
        
        #we should crop the data at least 2.5 ps before the main peak
        nearestdistancetopeak=2.5e-12
        if min(tdData[:,0])+timePreceedingSignal>-nearestdistancetopeak:
            timePreceedingSignal=-min(tdData[:,0])-nearestdistancetopeak
            print("Time Interval of preceeding noise to long, reset done")
            #dont use user input if higher than peak
        
        #get the first time
        starttime=min(tdData[:,0])
        if timePreceedingSignal<0:
            #determine length automatically            
            ratio=2
            ix_max=py.argmax(tdData[:,1])
            ix_min=py.argmin(tdData[:,1])
            earlier=min(ix_max,ix_min)
                #take only the first nts part            
            endtime=tdData[int(earlier/ratio),0]
        else:
            endtime=starttime+timePreceedingSignal

        #cut the data from starttime to endtime            
        noise=self.getShorterData(tdData,starttime,endtime)
        #detrend the noise
        noise=signal.detrend(noise[:,1:3],axis=0)
        return noise

    def getSNR(self):
        #returns the SNR
        return abs(self.getEX())/self.getUncEX()

    def getShorterData(self,tdData,tmin,tmax):
        #cuts tdData from time tmin to time tmax
        ix=py.all([tdData[:,0]>=tmin,tdData[:,0]<tmax],axis=0)
        return tdData[ix,:]
    
    def getTimes(self):
        #returns the time colon
        return self.tdData[:,0]

    def getTimesPs(self):
        #returns the time colon in ps
        return self.getTimes()*1e12

    def getTimeWindowLength(self):
        #returns thetime from the signal peak to the end of the measurement
        peak=self.getPeakPosition()
        return abs(self.tdData[-1,0]-peak)
    
    def getUncEX(self):
        #returns the uncertainty of the X channel
        return self.tdData[:,3]
    
    def getUncEY(self):
        #returns the uncertainty of the Y channel
        return self.tdData[:,4]
    
    def getWindowedData(self,windowlength_time=1e-12):
        #returns the blackmanwindowed tdData
        N=int(windowlength_time/self.dt)
        w=py.blackman(N*2)
        w=py.hstack((w[0:N],py.ones((self.getLength()-N*2),),w[N:]))
        windowedData=self.tdData
        #this could be written more elegantly?!
        windowedData[:,1]*=w
        windowedData[:,2]*=w
        windowedData[:,3]*=w
        windowedData[:,4]*=w
        return windowedData
    
    def getPulseWidth(self):
        pass    

    def importfile(self,fname,params):
        # if even more sophisticated things are needed, just inherit THzTdData class
        #and override the importfile method
        #try to load the file with name fname
        #it should be possible to write this shorter
        print(fname)
        try:
            #if no Y_col is specified            
            if 'Y_col' in params:
                #import it right away
                if params['dec_sep']=='.':
                    data=py.loadtxt(fname,
                                usecols=(params['time_col'],
                                         params['X_col'],
                                         params['Y_col']),
                                skiprows=params['skiprows'])
                                
                elif params['dec_sep']==',':
                    #if the decimal separator is , do a replacement
                    str2float=lambda val: float(val.replace(',','.'))
                    data=py.loadtxt(fname,
                                converters={params['time_col']:str2float,
                                            params['X_col']:str2float,
                                            params['Y_col']:str2float},
                                usecols=(params['time_col'],
                                         params['X_col'],
                                         params['Y_col']),
                                skiprows=params['skiprows'])                
            else:
                #import it right away
                if params['dec_sep']=='.':
                    data=py.loadtxt(fname,
                                usecols=(params['time_col'],
                                         params['X_col']),
                                skiprows=params['skiprows'])
                                
                elif params['dec_sep']==',':
                    #if the decimal separator is , do a replacement
                    str2float=lambda val: float(val.replace(',','.'))
                    data=py.loadtxt(fname,
                                converters={params['time_col']:str2float,
                                            params['X_col']:str2float},
                                usecols=(params['time_col'],
                                         params['X_col']),
                                skiprows=params['skiprows'])
                dummy_Y=py.zeros((data.shape[0],1))
                data=py.column_stack((data,dummy_Y))
        except IOError:
            print("File " + fname + " could not be loaded")
            sys.exit()

        #scale the timaaxis
        data[:,0]*=params['time_factor']
        
        #if the measurement was taken in negative time direction, flip the data
        if data[1,0]-data[0,0]<0:
            data=py.flipud(data)
     
        return data
    
    def resetTDData(self):
        #recalculate the mean tdData
        #do some data preprocessing of the rawdata
        processedData=self._processRawData(self._thzdata_raw)
        #calculate the mean
        meanData=self.calcTDData(processedData)
        #set the tdData to the calculated meanData
        self.setTDData(meanData)
 
    def setTDData(self,tdData):
        #set the tdData, use this function! 
        self.tdData=tdData
        #we lost some steps,this is why the mean occurs, maybe snd line is better?
        self.dt=abs(py.mean(tdData[10:20,0]-tdData[9:19,0]))       
#        self.dt=(tdData[-1,0]-tdData[0,0])/len(tdData[:,0])
        self.num_points=len(tdData[:,0])

    def zeroPaddData(self,desiredLength,paddmode='zero',where='end'):    
        #zero padds the time domain data, it is possible to padd at the beginning,
        #or at the end, and further gaussian or real zero padding is possible        
        #might not work for gaussian mode!

        desiredLength=int(desiredLength)
        #escape the function        
        if desiredLength<0:
            return 0

        #calculate the paddvectors        
        if paddmode=='gaussian':
            paddvec=py.normal(0,py.std(self.getPreceedingNoise())*0.05,desiredLength)
        else:
            paddvec=py.ones((desiredLength,self.tdData.shape[1]-1))
            paddvec*=py.mean(self.tdData[-20:,1:])
            
        timevec=self.getTimes()
        if where=='end':
            #timeaxis:
            newtimes=py.linspace(timevec[-1],timevec[-1]+desiredLength*self.dt,desiredLength)
            paddvec=py.column_stack((newtimes,paddvec))
            longvec=py.row_stack((self.tdData,paddvec))
        else:
            newtimes=py.linspace(timevec[0]-(desiredLength+1)*self.dt,timevec[0],desiredLength)
            paddvec=py.column_stack((newtimes,paddvec))
            longvec=py.row_stack((paddvec,self.tdData))
            
        self.setTDData(longvec)
        
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
        
        #interpolate the data
        for i in range(self.numberOfDataSets):
            tdDatas[i]=self.getInterData(tdDatas[i],commonLENGTH,commonMIN,commonMAX)
        
        return py.asarray(tdDatas)
    
    def _determineLockinPhase(self,rawtdData):
        #determine the phase difference from X and Y channel (at maximum signal)
        
        ix_max=py.argmax(rawtdData[:,1])
        #take 9 datapoints to evaluate theta
        no=4
        XCs=rawtdData[max(0,ix_max-no):min(rawtdData.shape[0],ix_max+no),1]
        YCs=rawtdData[max(0,ix_max-no):min(rawtdData.shape[0],ix_max+no),2]
        return py.arctan(py.mean(YCs/XCs))
    
    def _processRawData(self,tdDatas):
        tempTDDatas=[]
        for tdData in tdDatas:
            #this rotates all signal to X, adds X and Y uncertainty to each
            #tdData
            t=self._rotateToXChannel(tdData)
            #this removes Linear Drifts in X-Channel
            t=self._removeLinearDrift(t)
            
            tempTDDatas.append(t)

        #before interpolating to a common time axis, this need to be a list of 
        #tdData arrays, since they might differ in length
        
        #first shift maxima on top, than interpolate, doesn't affect unc array
        tempTDDatas,time_jitter=self._removeTimeShift(tempTDDatas)
        #also uncertainty is interpolated
        tempTDDatas=self._bringToCommonTimeAxis(tempTDDatas)
        return tempTDDatas

    def _removeLinearDrift(self,tdData):
        #do this for x and y channel?
        #overthink use of detrend here!
        tdData[:,1:3]=signal.detrend(tdData[:,1:3],axis=0)
        #take care in unsymmetric pulses: (hard 20 ?!)
#        tdData[:,1:3]=tdData[:,1:3]-py.mean(tdData[:20,1:3])
        return tdData
    
    def _removeTimeShift(self,tdDatas):
        #not sure if needed, maybe we want to correct the rawdata by shifting the maxima on top of each other
        #the indices of the maxima        
        #takes at the moment only the X Channel data and corrects it (safer!)
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

    def _rotateToXChannel(self,tdData):
        #this function should remove all signal from Y-Channel
        
        #Calculate lock-in phase:
        unc_raw=self.getelnNoise(tdData)
        XC=unumpy.uarray(tdData[:,1],unc_raw[0])
        YC=unumpy.uarray(tdData[:,2],unc_raw[1])
        #go to pulse:
        phase=self._determineLockinPhase(tdData)
        
        #rotate to XChannel
        XC_new=XC*py.cos(phase)+YC*py.sin(phase)
        YC_new=-XC*py.sin(phase)+YC*py.cos(phase)
        tdData[:,1]=unumpy.nominal_values(XC_new)
        tdData[:,2]=unumpy.nominal_values(YC_new)
        
        unc_new=py.column_stack((unumpy.std_devs(XC_new),unumpy.std_devs(YC_new)))
        tdData=py.column_stack((tdData,unc_new))
        return tdData

        
class ImportMarburgData(THzTdData):
    #only an example how a importer could look like,
    #in case of inrim and Marburg data, not really needed, (just define params)
    def __init__(self,filename):
       
        params={'time_factor':1,
                'time_col':0,
                'X_col':1,
                'Y_col':2,
                'dec_sep':',',
                'skiprows':0}    
        THzTdData.__init__(self,filename,params)
        
class ImportInrimData(THzTdData):
  
    def __init__(self,filename):
        params={'time_factor':1,
                'time_col':2,
                'X_col':3,
                'Y_col':5,
                'dec_sep':'.',
                'skiprows':0}    
        THzTdData.__init__(self,filename,params)
  

class FdData():
    '''A general fourier data class
    for initialization pass a time domain data object
    '''
    #Class variables, that restrict the accesible frequency data from 0 to 5 THz
    FMIN=0
    FMAX=10e12    
    
    def __init__(self,tdData,fbins=-1,fbnds=[FMIN,FMAX]):
        #FdData is always attached to a TdData Object
        
        #the original _tdData object should be private and not accessed from outside
        #if you want to know it use FdData.getassTDData()
        self._tdData=tdData
        
        #calculate the fft and store it to the fdData array
        self.fdData=self._calculatefdData(self._tdData) 
        
        #crop it to user bounds fbnds (min freq, max freq) and maybe also to a 
        #desired frequency step width fbins
        self.resetfdData(fbins,fbnds)
        self.maxDR=max(self.getDR())

    def _calculatefdData(self,tdData):
        #no need to copy it before (access member variable is ugly but shorter)

        #calculate the fft of the X channel
        fd=py.fft(tdData.getEX())
        #calculate absolute and phase values
        fdabs=abs(fd)
        fdph=abs(py.unwrap(py.angle(fd)))
        #calculate frequency axis
        dfreq=py.fftfreq(tdData.num_points,tdData.dt)
        #extrapolate phase to 0 and 0 frequency
        fdph=self.removePhaseOffset(dfreq,fdph)
        #set up the fdData array
        t=py.column_stack((dfreq,fd.real,fd.imag,fdabs,fdph))
        
        #this is so not nice!
        self.fdData=t
      
        #this seems to take very long! so maybe we skip it, until needed!
        unc=self.calculateFDunc()
        
        t=self.getcroppedData(t,0,unc[-1,0])
        
        #interpolate uncertainty        
        intpunc=interp1d(unc[:,0],unc[:,1:],axis=0)        
        unc=intpunc(t[:,0])
        return py.column_stack((t,unc))

    def calculateFDunc_old(self):
        #Calculates the uncertainty of the FFT according to:
        #   - J. M. Fornies-Marquina, J. Letosa, M. Garcia-Garcia, J. M. Artacho, "Error Propagation for the transformation of time domain into frequency domain", IEEE Trans. Magn, Vol. 33, No. 2, March 1997, pp. 1456-1459
        #return asarray _tdData
        #Assumes tha the amplitude of each time sample is statistically independent from the amplitude of the other time
        #samples

        # Calculates uncertainty of the real and imaginary part of the FFT and ther covariance
        unc_E_real = []
        unc_E_imag = []
        cov = []
        
       
        for f in self.getfreqs():
            unc_E_real.append(py.sum((py.cos(2*py.pi*f*self._tdData.getTimes())*self._tdData.getUncEX())**2))
            unc_E_imag.append(py.sum((py.sin(2*py.pi*f*self._tdData.getTimes())*self._tdData.getUncEX())**2))
            cov.append(-0.5*sum(py.sin(4*py.pi*f*self._tdData.getTimes())*self._tdData.getUncEX()**2))
      
        unc_E_real = py.sqrt(py.asarray(unc_E_real))
        unc_E_imag = py.sqrt(py.asarray(unc_E_imag))
        cov = py.asarray(cov)
        # Calculates the uncertainty of the modulus and phase of the FFT
        unc_E_abs = py.sqrt((self.getFReal()**2*unc_E_real**2+self.getFImag()**2*unc_E_imag**2+2*self.getFReal()*self.getFImag()*cov)/self.getFAbs()**2)
        unc_E_ph = py.sqrt((self.getFImag()**2*unc_E_real**2+self.getFReal()**2*unc_E_imag**2-2*self.getFReal()*self.getFImag()*cov)/self.getFAbs()**4)
        
        t=py.column_stack((self.getfreqs(),unc_E_real,unc_E_imag,unc_E_abs,unc_E_ph))
        return self.getcroppedData(t)  
    
    def calculateFDunc(self):
        unc_ex_sq=py.asmatrix(self._tdData.getUncEX()**2)
        times=py.asmatrix(self._tdData.getTimes())
        freqs=py.asmatrix(self.getfreqs())

        unc_E_real=py.asmatrix(py.asarray(py.cos(2*py.pi*freqs.transpose()*times))**2)*unc_ex_sq.transpose()
        unc_E_imag=py.asmatrix(py.asarray(py.sin(2*py.pi*freqs.transpose()*times))**2)*unc_ex_sq.transpose()
        cov=(-0.5*py.sin(4*py.pi*freqs.transpose()*times)*unc_ex_sq.transpose())

        
        unc_E_real = py.sqrt(py.asarray(unc_E_real))[:,0]
        unc_E_imag = py.sqrt(py.asarray(unc_E_imag))[:,0]
        cov = py.asarray(cov)[:,0]

        # Calculates the uncertainty of the modulus and phase of the FFT
        unc_E_abs = py.sqrt((self.getFReal()**2*unc_E_real**2+self.getFImag()**2*unc_E_imag**2+2*self.getFReal()*self.getFImag()*cov)/self.getFAbs()**2)
        unc_E_ph = py.sqrt((self.getFImag()**2*unc_E_real**2+self.getFReal()**2*unc_E_imag**2-2*self.getFReal()*self.getFImag()*cov)/self.getFAbs()**4)

        t=py.column_stack((self.getfreqs(),unc_E_real,unc_E_imag,unc_E_abs,unc_E_ph))
        return self.getcroppedData(t)  
    

    def doPlot(self):
        #plot the absolute and phase of the fdData array
        py.figure('FD-ABS-Plot')
        py.plot(self.getfreqsGHz(),20*py.log10(self.getFAbs()))
        py.xlabel('Frequency in GHz')
        py.ylabel('Amplitude, arb. units')
        
        py.figure('FD-PHASE-Plot')
        py.plot(self.getfreqsGHz(),self.getFPh())   
             
    def doPlotWithUnc(self):
        #plot absolute and phase along with uncertainties
        f=1

        uabs=unumpy.uarray(self.getFAbs(),self.getFAbsUnc())
        #scale the uncertainty for the 20*log10 plot
        uabs=20*unumpy.log10(uabs)
        u_a=unumpy.nominal_values(uabs)
        u_s=unumpy.std_devs(uabs)
        
        py.figure('FD-ABS-UNC-Plot')
        py.plot(self.getfreqsGHz(),u_a)
        py.plot(self.getfreqsGHz(),u_a+u_s,'--',self.getfreqsGHz(),u_a-u_s,'--')

        py.figure('FD-PH-UNC-Plot')
        py.plot(self.getfreqsGHz(),self.getFPh())
        py.plot(self.getfreqsGHz(),self.getFPh()+f*self.getFPhUnc(),'--',self.getfreqsGHz(),self.getFPh()-self.getFPhUnc(),'--')

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
        #this method might lead to better results than the findAbsorptionLines method?
        #movav=self.getmovingAveragedData()        
        hlog=-20*py.log10(self.getFAbs())
        etalon=self.getEtalonSpacing()
        Ns=int(etalon/self.getfbins())
        Ns=py.arange(max(1,Ns-10),Ns+10,1)        
        peaks=signal.find_peaks_cwt((hlog),Ns)
        return peaks
        
    def fitAbsorptionLines(self):
        #this method should fit the absorption lines and return fwhm, depth, unfinished
        peaks=self.findAbsorptionLines()
        #naive idea: use always 1/n of the space between two peaks for cutting:
        peakfreqs=self.getfreqs()[peaks]
        #first and last one are special
        #first one:
        if len(peaks)==0:
            #nothing to do 
            return 0
        
        #if len(peaks)==1:
            #there is only one
            
        #if len(peaks)==2:
            
        #if len(peaks)>2:
        #(peakfreqs[1]-peakfreqs[0])
    
    def getassTDData(self):
        #return the underlying tdData object
        return self._tdData

    def getBandwidth(self,dbDistancetoNoise=5):
        #this function should return the lowest trustable and highest trustable
        #frequency, along with the resulting bandwidth
        
        absdata=20*py.log10(self.getFAbs()/max(self.getFAbs()))
        maxSNR=20*py.log10(max(self.getSNR()))
        ix=maxSNR-dbDistancetoNoise>-absdata #dangerous, due to sidelobes, there might be some high freq component!
        tfr=self.getfreqs()[ix]
        dtfr=py.diff(tfr)
        frequencygaps=py.where(dtfr>100e9)[0]
        if len(frequencygaps)==0:
            itemix=-1
        else:
            itemix=frequencygaps[0]
        
        return [min(tfr),tfr[itemix]]

    def getcroppedData(self,data,startfreq=FMIN,endfreq=FMAX):
        #this function returns the array data cropped from startfreq to endfreq
        ix=py.all([data[:,0]>=startfreq,data[:,0]<=endfreq],axis=0)
        return data[ix,:]

    def getDR(self):
        #this function should return the dynamic range
        #this should be the noiselevel of the fft
        noiselevel=py.sqrt(py.mean(abs(py.fft(self._tdData.getAllPrecNoise()[0]))**2))
        #apply a moving average filter on log
        window_size=5
        window=py.ones(int(window_size))/float(window_size)
        hlog=py.convolve(20*py.log10(self.getFAbs()), window, 'valid')
        one=py.ones((2,))
        hlog=py.concatenate((hlog[0]*one,hlog,hlog[-1]*one))
        return hlog-20*py.log10(noiselevel)         

    def getEtalonSpacing(self):
        #this should return the frequency of the Etalon
    
    
        #how to find a stable range!
        bw=self.getBandwidth()
        rdata=self.getcroppedData(self.fdData,max(bw[0],250e9),min(bw[1],2.1e12))  
        
        #need to interpolate data!
        oldfreqs=rdata[:,0]
        intpdata=interp1d(oldfreqs,rdata[:,3],'cubic')
        
        fnew=py.arange(min(oldfreqs),max(oldfreqs),0.1e9)
        absnew=intpdata(fnew)
        #find minimia and maxima        
        ixmaxima=signal.argrelmax(absnew)[0]
        ixminima=signal.argrelmin(absnew)[0]
        
        fmaxima=py.mean(py.diff(fnew[ixmaxima]))
        fminima=py.mean(py.diff(fnew[ixminima]))
        #calculate etalon frequencies
        df=(fmaxima+fminima)*0.5 #the etalon frequencies
        print(str(df/1e9) + " GHz estimated etalon frequency")
        return df

    def getFReal(self):
        #the real part of the fft(tdData)
        return self.fdData[:,1]
    def getFImag(self):
        #the imag part of the fft(tdData)
        return self.fdData[:,2]
    def getFAbs(self):
        #the absolute value of fft(tdData)
        return self.fdData[:,3]
    def getFPh(self):
        #the phase of the fourierdomaindata
        return self.fdData[:,4]
    def getFRealUnc(self):
        #the uncertainty of the real part of the fft
        return self.fdData[:,5]
    def getFImagUnc(self):
        #the uncertainty of the imag part of the fft        
        return self.fdData[:,6]
    def getFAbsUnc(self):
        #the uncertainty of the absolute value
        #maybe the calculation can also be done here, since it is not needed often
        return self.fdData[:,7]
    def getFPhUnc(self):
        #the uncertainty of the phase value
        return self.fdData[:,8]

    def getfbins(self):
        #the frequency spacing
        return abs(self.fdData[1,0]-self.fdData[0,0])

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
        
        absdata=signal.savgol_filter(self.getFAbs(),N_min-N_min%2+1,order)
        phdata=signal.savgol_filter(self.getFAbs(),N_min-N_min%2+1,order)
       
        return py.column_stack((self.fdData[:,:3],absdata,phdata,self.fdData[:,5:]))
    
    def getfreqs(self):
        #return the frequency axis
        return self.fdData[:,0]

    def getfreqsGHz(self):
        #return the frequency axis in GHz
        return self.getfreqs()*1e-9         

    def getInterpolatedFDData(self,newfbins,strmode='linear'):
        #not only for this function, also for others maybe useful: 
        #check if it is possible to give as an argument a slice, i.e.
        # '1:' for interpolating all data, '3' for only the absdata and so on 
        oldfreqs=self.getfreqs()
        newfreqs=py.arange(min(oldfreqs),max(oldfreqs),newfbins)
        
        interfdData=interp1d(oldfreqs,self.fdData[:,1:],strmode,axis=0)
        
        return py.column_stack((newfreqs,interfdData(newfreqs)))

    def getLength(self):
        #the length of the fdData array
        return self.fdData.shape[0]
    
    def getmaxfreq(self):
        #take care of constant offset what is the best lower range?
        cutted=self.getcroppedData(self.fdData,150e9,FdData.FMAX)
        fmax=cutted[py.argmax(cutted[:,3]),0]
        return fmax

    def getmovingAveragedData(self,window_size_GHz=-1):
        #so far unelegant way of convolving the columns one by one
        #even not nice, improvement possible?
        if window_size_GHz<0.5e9:
            window_size=int(self.getEtalonSpacing()/self.getfbins())
        else:
            window_size=int(window_size_GHz/self.getfbins())
              
        window_size+=window_size%2+1
        window=py.ones(int(window_size))/float(window_size)
        
        dataabs=py.convolve(self.getFAbs(), window, 'valid')
        dataph=py.convolve(self.getFPh(), window, 'valid')
        one=py.ones((window_size-1)/2,)
        dataabs=py.concatenate((dataabs[0]*one,dataabs,dataabs[-1]*one))
        dataph=py.concatenate((dataph[0]*one,dataph,dataph[-1]*one))
        return py.column_stack((self.fdData[:,:3],dataabs,dataph,self.fdData[:,5:]))
    
    def getSNR(self):
        #returns the signal to noise ratio
        return self.getFAbs()/self.getFAbsUnc()

    def removePhaseOffset(self,freqs,ph,startfreq=200e9,endfreq=1e12):
        #cut phase to reasonable range:
        ph_c=self.getcroppedData(py.column_stack((freqs,ph)),startfreq,endfreq)
        #determine the slope and the offset      
        p=py.polyfit(ph_c[:,0],ph_c[:,1],1)
        #return full phase-offset(croppedPhase)
        return ph-p[1]
     
    def resetfdData(self,fbins=-1,fbnds=[FMIN,FMAX]):
        #crop and zeropadd the fdData, starting from _tdData again!
        if fbins>0 and self.getfbins()>fbins:
            self.zeroPadd(fbins)
            
        self.setFDData(self.getcroppedData(self.fdData,fbnds[0],fbnds[1]))
      
    def setTDData(self,tdData):
        #Neccessary? couldn't we just create a new object with init?
        fbins_old=self.getfbins()
        minf_old=min(self.getfreqs())
        maxf_old=max(self.getfreqs())
        self._tdData=tdData        
        self.setFDData(self._calculatefdData(self._tdData))
        
        self.resetfdData(fbins_old,[minf_old,maxf_old])

    def setFDData(self,fdData):
        #sets the fdData array, obsolete?
        self.fdData=fdData

    def setPhase(self,newPh):
        #sets the Phase, needed for example for removing offsets
        if len(newPh)!=self.getLength():
            print( 'Setting phase not possible, wrong length')
        else:
            self.fdData[:,4]=newPh

    def zeroPadd(self,fbins):
        #zero padd the underlying tdData such that the fbins afterwards are fbins
        spac=1/self._tdData.dt/fbins
        actlen=self._tdData.getLength()
        nozeros=py.ceil(spac-actlen)
        self._tdData.zeroPaddData(nozeros)
        #leave the old bnds in place
        bnds=[min(self.getfreqs()),max(self.getfreqs())]
        zpd=self._calculatefdData(self._tdData)
        self.setFDData(self.getcroppedData(zpd,bnds[0],bnds[1]))
        
if __name__ == '__main__':
    myTDData=ImportInrimData(['2015-02-19T140227_fw0_.dat'])
    myFDData=FdData(myTDData)