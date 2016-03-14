import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.signal as signal
from scipy.interpolate import UnivariateSpline

from uncertainties import unumpy
import datetime

class TimeDomainData():
    '''
    A Class for representing the measured THz-Time Domain Data. 
    Often used operations on this data are implemented.

    Params
    ----------------
    timeaxis: numpy array
        The time axis. 
    efield: numpy array
        The measured Efield data. (should have same length as timeaxis)
    uncertainty: numpy array
        The uncertainty on the efield, Default=None
    datasetname:
        Name of the Dataset, Default=''
    propagateUncertainty:
        If set to true, all calculations will be performed using the
        uncertainty package, this will propagate the uncertainty but it will also 
        make operation with the time domain data slow, Default=False
    
    Example
    ------------------
    >>>t=np.arange(0,100e-12,30e-15)
    >>>P=(t-10e-12)*np.exp(-(t-10e-12)**2/0.15e-24)
    >>>td=TimeDomainData(t,P)
    >>>td.plotme()
    '''
    
    def fromFile(filename,decimalSeparator='.',dataColumns=[0,1],
                skipRows=0,timefactor=1,propagateUncertainty=False,flipData=False):
        '''
        Loads the file filename with the given format and returns an object of TimeDomainData Type
        If a YChannel is given the signal is rotated to only one Channel
        
        Parameters
        ---------------
        filename: string
            The name of the file to read.
        decimalSeparator: string
            The decimalö separator used, typically either . or , (optional)
        datacolumns: array of ints
            The column number of the time data, XChannel and optionally YChannel
        skipRows: int
            The number of rows to skip
        timefactor: float
            The time axis should be in units of (s), else input here the correct factor.
        propagateUncertainty: bool
            Switch uncertainty propagation on or off
        flipData: bool
            If the pulse arrives from late times (reversed delay line) the data can be flipped

        Returns
        -------------
        Time Domain Data object
        '''
        try:
            str2float=lambda val: float(val.decode("utf-8").replace(decimalSeparator,'.'))
            #Import also YChannel, else neglect it
            if len(dataColumns)>2 and dataColumns[2]>0:
                #import it right away
                data=np.loadtxt(filename,
                                converters={dataColumns[0]:str2float,
                                            dataColumns[1]:str2float,
                                            dataColumns[2]:str2float},
                                usecols=(dataColumns[0],dataColumns[1],dataColumns[2]),
                                skiprows=skipRows)
                timeaxis=data[:,0]
                efield=TimeDomainData._rotateToXChannel(data[:,1],data[:,2])
            else:
                #import it right away
                data=np.loadtxt(filename,
                                converters={dataColumns[0]:str2float,
                                            dataColumns[1]:str2float},
                                usecols=(dataColumns[0],dataColumns[1]),
                                skiprows=skipRows)
                timeaxis=data[:,0]
                efield=data[:,1]
        except IOError:
            print("File " + filename + " could not be loaded")
            return 0
            
        #scale the timaaxis
        timeaxis*=timefactor
        
        #if the measurement was taken in negative time direction, flip the data
        #INRIM Mode
        if timeaxis[1]-timeaxis[0]<0:
            timeaxis=timeaxis[::-1]
            efield=efield[::-1]
        
        #if the measurement was taken towards smaller times, Marburg case
        #if np.argmax(efield)>len(efield)/2:
        if flipData:
            #assume now that the data must get flipped, this must become optional
            efield=efield[::-1]
            
        if propagateUncertainty:
            sigma_BG=TimeDomainData.estimateBGNoise(timeaxis,efield)
        else:
            sigma_BG=None
        return TimeDomainData(timeaxis,efield,sigma_BG,filename,propagateUncertainty)
    
    def _rotateToXChannel(XChannel,YChannel):
        '''Removes all signal from Y-Channel, determines the phase difference from X and Y channel (at maximum signal)
            (discussion: Should for different times different thetas be allowed?)
        
            Parameters
            --------------
            XChannel: array
                The Xchannel raw data.
            YChannel: array
                The YChannel raw data.
            
            Returns
            ----------------
            XChannelNew: array
                The rotated data, YChannel should be empty now
        '''
        
        ix_max=np.argmax(XChannel)
        #take 9 datapoints to evaluate theta
        no=4
        XCs=XChannel[max(0,ix_max-no):min(len(XChannel),ix_max+no)]
        YCs=YChannel[max(0,ix_max-no):min(len(YChannel),ix_max+no)]
        phase=np.arctan(np.mean(YCs/XCs))

        #rotate to XChannel
        XC_new=XChannel*np.cos(phase)+YChannel*np.sin(phase)
        #YC_new=-XC*.sin(phase)+YC*py.cos(phase)
        return XC_new 
    
    def fromFrequencyDomainData(fdata):
        '''
        Uses a FrequencyDomainData object to calculate the TimeDomainData
        via inverse Fourier-Transform
        
        Parameter
        ----------
        fdata: FrequencyDomainData
            The FrequencyDomainData object from which the timeDomain signal should be calculated
        
        Returns
        --------------
        The corresponding TimeDomainData Object.
        '''
        #if no negative frequency components are present (standard case, generate the correct data length)
        timetrace=np.fft.irfft(fdata.getSpectrumRef())
        timeaxis=np.fft.rfftfreq(len(timetrace)*2-1,fdata.getfbins()/2)
        #so far uncertainty gets lost :-()        
        
        return TimeDomainData(timeaxis,timetrace)
        
    def estimateBGNoise(timeaxis,efield,timePreceedingSignal=-1):
        '''
        Estimate of the Background Noise

        Parameter
        -----------
        timeaxis: array
            The timeaxis of the timeDomain signal
        efield: array
            The time domain trace efield
        timePreceedingSignal: float
            The time until which the signal is considered to be noise
            (optional), if not given: half of the signal before the pulse is used

        Returns
        ----------
        The standard deviation of the Background Noise
        '''
        #we should crop the data at least 2.5 ps before the main peak
        nearestdistancetopeak=2.5e-12
        if min(timeaxis)+timePreceedingSignal>-nearestdistancetopeak:
            timePreceedingSignal=-min(timeaxis)-nearestdistancetopeak
            print("Time Interval of preceeding noise to long, reset done")
            #dont use user input if higher than peak

        #get the first time
        starttime=min(timeaxis)
        if timePreceedingSignal<0:
            #determine length automatically            
            ratio=2
            ix_max=np.argmax(efield)
            ix_min=np.argmin(efield)
            earlier=min(ix_max,ix_min)
                #take only the first nts part
            
            endtime=timeaxis[int(earlier/ratio)]
        else:
            endtime=starttime+timePreceedingSignal

        #cut the data from starttime to endtime            
        ix=np.all([timeaxis>=starttime,timeaxis<endtime],axis=0)
        noise=efield[ix]
        #detrend the noise
        noise=signal.detrend(noise)
        std_BGNoise=np.std(noise)
        return std_BGNoise
    
    def averageTimeDomainDatas(timeDomainDatas,propagateUncertainty=False):
        '''
        Averages the list of TimeDomainData Objects. If the time axis of the 
        elements in timeDomainDatas is not equal the timeDomainDatas will be interpolated 
        to a common time axis.
        
        Parameters
        ------------
        timeDomainDatas: list
            A list of timeDomainData objects
        propagateUncertainty: bool
            Switch to turn on the uncertainty propagation, default False
        
        Returns
        -------------
        TimeDomainData Object
        '''
        #in case of different data length the time axis will be interpolated
        timeDomainDatas=TimeDomainData._bringToCommonTimeAxis(timeDomainDatas)
        
        efields=[tdd.getEfield() for tdd in timeDomainDatas]
        av=np.average(efields,axis=0)
        std=np.std(efields,axis=0,ddof=1)/np.sqrt(len(timeDomainDatas))
        name='average'
        for tdd in timeDomainDatas:
            name+='_' + tdd.getDataSetName()
        
        return TimeDomainData(timeDomainDatas[0].getTimeAxis(),av,std,name,propagateUncertainty) 
            
    def _bringToCommonTimeAxis(timeDomainDatas,force=True):
        '''
        Interpolates the list of timeDomainDatas to a common time axis. The new time axis is from the smallest 
        common to the largest common time value. 
        
        Parameter
        -------------
        timeDomainDatas: list
            A list of timeDomainData objects.
            
        Returns
        ------------
            A list of timeDomainData objects that have the same time axis.
        '''
        dt=timeDomainDatas[0].getTimeStep()
       
        timeaxisarray=[tdd.getTimeAxisRef() for tdd in timeDomainDatas]
        samelength=all(tdd.getSamplingPoints()==timeDomainDatas[0].getSamplingPoints() for tdd in timeDomainDatas)
        if samelength and np.all((abs(timeaxisarray-timeaxisarray[0]))<dt/1000):
            #no interpolation neeeded
            return timeDomainDatas
        else:
            print("Time axis manipulation neccessary")
            miss_points_max=10
            #check for missing datapoints, allowing 
            #not more than miss_points_max points to miss 
            all_lengthes=[]
            for tdd in timeDomainDatas:
                all_lengthes.append(tdd.getSamplingPoints())
            
            #do it always, just print a warning, if miss_points_max is exceeded
            if min(all_lengthes)!=max(all_lengthes):
                print("Datalength of suceeding measurements not consistent, try to fix")
                if max(all_lengthes)-min(all_lengthes)>miss_points_max:
                    print("Warning: Data seems to be corrupted. \n" +\
                    "The length of acquired data of repeated measurements differs by \n" + \
                        str(max(all_lengthes)-min(all_lengthes)) + ' datapoints')
                    return 0
            #interpolation does no harm, even if everything is consistent (no interpolation in this case)
            commonMin=max([tdd.getTimeAxisRef().min() for tdd in timeDomainDatas])
            commonMax=min([tdd.getTimeAxisRef().max() for tdd in timeDomainDatas])
            commonLength=min(all_lengthes)
            
            #interpolate the data
            commonAxisDatas=[]
            for tdd in timeDomainDatas:
                commonAxisDatas.append(tdd.getInterpolatedTimeDomainData(commonLength,commonMin,commonMax))
            return commonAxisDatas
        
    def _removeTimeShift(timeDomainDatas):
        '''Shifts the maxima of several pulses on top of each other to correct for time jitters
        
        Parameter
        ------------
        timeDomainDatas: list
            A list of timeDomainData objects.
            
        Returns
        --------------
        A list of time-shifted TimeDomainDatas
        '''
        if len(timeDomainDatas)>1:
            peak_pos=[]
            for tdd in timeDomainDatas:
                peak_pos.append(tdd.getPeakPositionInterpolated())
            
            peak_pos=np.asarray(peak_pos)
            print('Peak Position standard deviation: {:2.2f} fs'.format(np.std(peak_pos*1e15)))
            mp=np.mean(peak_pos)
            
            shiftedDatas=[]
            for tdd,peak in zip(timeDomainDatas,peak_pos):
                newTimeAxis=tdd.getTimeAxis()-(peak-mp)
                shiftedDatas.append(TimeDomainData(newTimeAxis,tdd.getEfield(),tdd.getUncertainty(),tdd.getDataSetName(),tdd.uncertaintyEnabled))
            
            return shiftedDatas
        else:
            return timeDomainDatas
            
    def _preProcessData(timeDomainDatas,average=True,removeTimeShift=True,removeDrift=True):
        '''Preprocessing of the list of raw timeDomainDatas

        Parameters
        ---------------        
        timeDomainDatas: list
            A list of timeDomainData Objects
        average: bool
            If True the timeDomainData will be averaged
        removeTimeShift: bool 
            If true the peak positions of the timeDomainDatas will be shifted
            on top of each other
        removeDrift: bool
            If true a constant and linear drift will be substracted from the Efield.
        
        Returns
        ---------------
        In all cases a list of preprocessed TimeDomainData objects
        '''
        
        tempTDDatas=timeDomainDatas
        
        if removeDrift:
            tempTDDatas=[]
            for tdd in timeDomainDatas:
                #this removes Linear Drifts in X-Channel
                tempTDDatas.append(tdd._removeLinearDrift())
                
        if removeTimeShift and len(tempTDDatas)>1:
            #first shift maxima on top, than interpolate, doesn't affect unc array
            tempTDDatas=TimeDomainData._removeTimeShift(tempTDDatas)
        
        if average and len(tempTDDatas)>1:
            tempTDDatas=TimeDomainData.averageTimeDomainDatas(tempTDDatas)
        
        return tempTDDatas
            
    def importMultipleFiles(fns,importparams={},preprocessing={}):
        '''Imports multiple files and does preprocessing
        
        Parameters
        ----------------
        fns: list
            The list of filenames to be imported
        importparams: dictionary
            (optional) Options to be passed to the fromFile function. Read the listing of fromFile for Details.
        preprocessing: dictionary
            (optional) Options to be passed to the _preProcessData function.
        
        Returns
        ----------------
        A single TimeDomainData Object or if averaging is disabled a list of TimeDomainData Objects.
        
        Example
        ---------------
        >>>importparams={  'timefactor':1,
                        'dataColumns':[0,1,2],
                        'decimalSeparator':',',
                        'skipRows':0,
                        'propagateUncertainty':False,
                        'flipData':False}
        >>>preprocessing={'average':True,'removeTimeShift':True,'removeDrift':True}
        '''
        datas=[]
        if isinstance(fns,list):
            for fn in fns:
                datas.append(TimeDomainData.fromFile(fn,**importparams))
        elif isinstance(fns,str):
            datas.append(TimeDomainData.fromFile(fns,**importparams))
        else:
            print('Either a list of filenames or a single filename should be entered')
            return 0
           
        av_data=TimeDomainData._preProcessData(datas,**preprocessing)
        if isinstance(av_data,list):
            #if no averaging than a list should be returned
            #in all other cases only the TimeDomainData object!            
            if len(av_data)>1:
                return av_data
            else:
                return av_data[0]
        else:
            return av_data
        
    def __init__(self,timeaxis,efield,uncertainty=None,datasetname='',propagateUncertainty=False):
        
        self.timeaxis=np.copy(timeaxis)
        if uncertainty is None:
            self.uncertainty=0#TimeDomainData.estimateBGNoise(timeaxis,efield)
        else:
            self.uncertainty=uncertainty
        
        if propagateUncertainty:              
            self.efield=unumpy.uarray(efield,self.uncertainty)
        else:
            self.efield=efield

        self.uncertaintyEnabled=propagateUncertainty
        self.datasetname=datasetname

    def enableUncertaintyPropagation(self):
        #if not isinstance(self.efield,unumpy.uarray):
        if self.uncertaintyEnabled==False:        
            #should be something already
            #if self.uncertainty==None:
            #    self.uncertainty=self.getSigmaBG()
            
            self.efield=unumpy.uarray(self.efield,self.uncertainty)
            self.uncertaintyEnabled=True
            
    def disableUncertaintyPropagation(self):
        #if isinstance(self.efield,unumpy.uarray):
        if self.uncertaintyEnabled==True:        
            self.uncertainty=unumpy.std_devs(self.efield)
            self.efield=unumpy.nominal_values(self.efield)
            self.uncertaintyEnabled=False

    def getDataSetName(self):
        '''Returns the Name of the Dataset'''
        return self.datasetname
        
    def setDataSetName(self,newname):
        '''Sets the name of the Dataset'''
        self.datasetname=newname

    def getDynamicRange(self):
        
        #returns the dynamic range
        Emax=max(self.getEfield())
        noise=self.getSigmaBG()
        return Emax/noise
        
    def getEfield(self):
        '''Returns the EField'''
        return unumpy.nominal_values(self.efield)

    def getFirstPuls(self,after=5e-12):
        '''returns only the first pulse'''
        
        #finds the main peak
        tcenter=self.getPeakPosition()
        #without argument, assume width of the pulse to be not more than 10ps!
        #further crop always from the beginning of the time data, 
        #in order to avoid phase problems
        tmax=tcenter+after
        tmin=min(self.getTimeAxisRef())
        
        return self.getTimeSlice(tmin,tmax)

    def getInterpolatedTimeDomainData(self,desiredLength,mint,maxt,tkind='linear'):
        #use cubic interpolation only, if you know, that the data is not too no
        intpdata=interp1d(self.getTimeAxisRef(),np.asarray([self.getEfield(),self.getUncertainty()]),axis=1,kind=tkind)
        timeaxis=np.linspace(mint,maxt,desiredLength,endpoint=True)
        longerData=intpdata(timeaxis)
        return TimeDomainData(timeaxis,longerData[0,:],longerData[1,:],self.getDataSetName(),self.uncertaintyEnabled)

    def getPeakPosition(self):
        '''gives the time, at which the signal is maximal'''
        return self.getTimeAxisRef()[np.argmax(self.getEfield())]

    def getPeakPositionInterpolated(self,newresolution=0.1e-15):
        '''Interpolates cubic'''
        peakPosition=self.getPeakPosition()
        NoSteps=5
        
        mi=peakPosition-NoSteps*self.getTimeStep()
        ma=peakPosition+NoSteps*self.getTimeStep()
        reduced=self.getTimeSlice(mi,ma)
        
        interpolater=interp1d(reduced.getTimeAxisRef(),reduced.getEfield(),kind='cubic')
        x_new=np.arange(min(reduced.getTimeAxisRef())+newresolution,max(reduced.getTimeAxisRef())-newresolution,newresolution)
                
        y_new=interpolater(x_new)
        
        return x_new[np.argmax(y_new)]

    def getPeakPositionExtrapolated(self,newresolution=0.1e-15):
        '''Extrapolates the data using splines in order to resolve the peak position better'''
        peakPosition=self.getPeakPosition()
        NoSteps=5
        
        mi=peakPosition-NoSteps*self.getTimeStep()
        ma=peakPosition+NoSteps*self.getTimeStep()
        reduced=self.getTimeSlice(mi,ma)
        extrapolator=UnivariateSpline(reduced.getTimeAxisRef(),reduced.getEfield(),k=5)
        x_new=np.arange(mi,ma,newresolution)
        y_new=extrapolator(x_new)
        return x_new[np.argmax(y_new)]

    def getPeakWidth(self):
        efield=abs(self.getEfield())
        mefield=np.amax(efield)
        #take care, maybe peak resolution not enough... interpolation?
        m=self.getTimeAxis()[efield>mefield/2]
        return m[-1]-m[0]

    def getSamplingPoints(self):    
        return self.getTimeAxisRef().shape[0]
  
    def getSigmaBG(self):
        return TimeDomainData.estimateBGNoise(self.timeaxis,unumpy.nominal_values(self.efield))

    def getSigmaRepeatability(self):
        pass

    def getSNR(self):
        #returns the SNR
        return abs(self.getEfield())/self.getUncertainty()

    def getTimeAxis(self):
        return np.copy(self.timeaxis)
        
    def getTimeAxisRef(self):
        return self.timeaxis
    
    def getTimeStep(self):
        return self.timeaxis[1]-self.timeaxis[0]

    def getTimeSlice(self,tmin,tmax):
        newtaxis=self.getTimeAxis()
        ix=np.all([newtaxis>=tmin,newtaxis<=tmax],axis=0)
        return TimeDomainData(newtaxis[ix],self.getEfield()[ix],self.getUncertainty()[ix],self.getDataSetName(),self.uncertaintyEnabled)

    def getTimeWindowLength(self):
        #returns thetime from the signal peak to the end of the measurement
        peak=self.getPeakPosition()
        return abs(self.getTimeAxisRef()[-1]-peak)
        
    def getUncertainty(self):
        return unumpy.std_devs(self.efield)
    
    def getWindowedData(self,windowlength_time=5e-12,windowtype='blackman'):
        '''Applies a window function to the TimeDomainData. 
        
        Parameter
        ---------------
        windowlength_time: float, default:-1
            The time (s) of the rising edge of the window. If negative, no flat top is used. Default: 5e-12
        windowtype: String or tuple
            Read scipy.signal.get_window for detailed description. All windows listed there can be used.
            Default: blackman
            
        Returns
        -------------        
        TimeDomainData Object with windowfunction applied. 
        '''
        if windowlength_time>0:
            #check that N is not too large! needs a fix here
            N=int(windowlength_time/self.getTimeStep())
            if 2*N>self.getSamplingPoints():
                N=int(self.getSamplingPoints()/2)
                print("Window too large")
            w=signal.get_window(windowtype,N*2)
            w=np.hstack((w[:N],np.ones((self.getSamplingPoints()-N*2),),w[N:]))
        else:
            w=signal.get_window(windowtype,self.getSamplingPoints())
        
        windowedData=self.getEfield()
        windowedData*=w
            
        return TimeDomainData(self.getTimeAxisRef(),unumpy.nominal_values(windowedData),unumpy.std_devs(windowedData),self.getDataSetName(),self.uncertaintyEnabled)
        
    def _removeLinearDrift(self):
        '''
        Does a scipy.signal.detrend on the timedomain data. Removes constant and linear drifts.
        Returns
        ---------------
        TimeDomainData Object without drift on Efield
        '''
        newfield=signal.detrend(unumpy.nominal_values(self.efield))
            
        return TimeDomainData(self.getTimeAxisRef(),newfield,self.getUncertainty(),self.getDataSetName(),self.uncertaintyEnabled)
  
    def plotme(self,ax=None,*plotarguments,**plotDict):
        '''
        Plots the TimeDomain trace in the axes ax. If no axes are given use current axes.

        Parameter
        ---------------
        ax: AxesSubplot
            The axes to plot into. Set to plt.gca() or None if not needed. 
        plotarguments: 
            Additional arguments to be passed to plot
        plotDict:
            Additional arguments to be passed to plot
        Returns
        --------
        The used axes object.
        '''
        if ax is None:
            ax=plt.gca()
        
        ax.plot(self.getTimeAxisRef()*1e12,self.getEfield(),*plotarguments,**plotDict)
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('Amplitude (arb. units)')
        
        return ax        
        
    def zeroPaddToTargetFrequencyResolution(self,fbins=5e9,paddmode='zero',where='end'):
        '''
        As many zeros are added as need to achieve a frequency resolution of fbins. If
        this is needed only for FourierTransform than this can be done also in the constructor of
        FrequencyDomainData.
        
        Parameter
        --------------
        fbins: float
            Frequency resolution (Hz) after zeropadding and Fourier Transform.
        paddmode: string
            Can be 'zero' (default) or 'gaussian'. Padds zeros or gaussian noise.
        where: string
            Can be 'end' (default) or 'start'. Padds the zeros before or after the pulse.

        Returns
        -----------
        TimeDomainData Object with new timeaxis         
        '''
        spac=1/self.getTimeStep()/fbins
        actlen=self.getSamplingPoints()
        nozeros=np.ceil(spac-actlen)
        return self.zeroPaddData(nozeros,paddmode,where)
        
    def zeroPaddData(self,desiredLength,paddmode='zero',where='end'):    
        '''
        Adds desiredLength number of zeros to the Efield Data. If this is needed only
        for FourierTransform than this can be done also in the constructor of FrequencyDomainData.
        
        Parameter
        --------------
        desiredLength: int
            Number of zeros to be padded.
        paddmode: string
            Can be 'zero' (default) or 'gaussian'. Padds zeros or gaussian noise.
        where: string
            Can be 'end' (default) or 'start'. Padds the zeros before or after the pulse.

        Returns
        -----------
        TimeDomainData Object with new timeaxis         
        ''' 
        ##seems not to work for before padding (makes not often sense but yet where is the problem ? needs a fix!)
        desiredLength=int(desiredLength)
        #escape the function        
        if desiredLength<0:
            return self

        #calculate the paddvectors        
        if paddmode=='gaussian':
            paddvec=np.random.normal(0,self.getSigmaBG(),desiredLength)
        else:
            paddvec=np.zeros((desiredLength,))
            
        timevec=self.getTimeAxis()
        
        
        if where=='end':
            #timeaxis:
            newtimes=np.linspace(timevec[-1],timevec[-1]+desiredLength*self.getTimeStep(),desiredLength)
            timevec=np.concatenate((timevec,newtimes))
            if self.uncertaintyEnabled:
                newefield=np.concatenate((self.getEfield(),unumpy.uarray(paddvec,self.getSigmaBG())))
            else:
                newefield=np.concatenate((self.getEfield(),paddvec))
        else:
            newtimes=np.linspace(timevec[0]-(desiredLength+1)*self.getTimeStep(),timevec[0],desiredLength)
            timevec=np.concatenate((newtimes,paddvec))
            if self.uncertaintyEnabled:
                newefield=np.concatenate((unumpy.uarray(paddvec,self.getSigmaBG()),self.getEfield()))
            else:
                newefield=np.concatenate(paddvec,self.getEfield())
        return TimeDomainData(timevec,unumpy.nominal_values(newefield),unumpy.std_devs(newefield),self.getDataSetName(),self.uncertaintyEnabled)

def importMarburgData(filenames):
    params={'timefactor':1,
            'dataColumns':[0,1,2],
            'decimalSeparator':',',
            'skipRows':0,
            'propagateUncertainty':False,
            'flipData':False}
    return TimeDomainData.importMultipleFiles(filenames,params)

def importINRIMData(filenames):
    params={'time_factor':1,
            'time_col':2,
            'X_col':3,
            'Y_col':5,
            'dec_sep':'.',
            'skiprows':0}    
    return TimeDomainData.importMultipleFiles(filenames,params)

def getTimeAndDate(filename):
    '''Assume that filename is of marburg format'''
    fn=filename.split('/')[-1]
    year=int(fn[:4])
    month=int(fn[4:6])
    day=int(fn[6:8])
    hour=int(fn[9:11])
    minute=int(fn[11:13])
    second=int(fn[13:15])
    
    return datetime.datetime(year,month,day,hour,minute,second)
    
    
def outputInformation(tdData,fdData,filename,header=False):
    myfile=open(filename,'a')
    if header==True:
        headerstr='Name; Time Step; Peak Width; Peak Position; Max Efield; '
        headerstr+='Peak Position; TD SNR; Frequency bins; Bandwidth; FD SNR'
        myfile.write(headerstr+'\n')
    writestr=''
    writestr+=tdData.getDataSetName()+'; '
    writestr+=str(tdData.getTimeStep())+'; '
    writestr+=str(tdData.getPeakWidth())+'; '
    writestr+=str(tdData.getPeakPosition())+'; '
    writestr+=str(np.amax(abs(tdData.getEfield())))+'; '
    writestr+=str(tdData.getPeakPosition())+'; '
    writestr+=str(np.amax(tdData.getSNR())) +'; '
    writestr+=str(fdData.getfbins()) + '; '
    writestr+=str(fdData.getBandwidth()) + '; '
    writestr+=str(np.amax(fdData.getSNR())) + '; '
    myfile.write(writestr + '\n')
    myfile.close()
        
        
    

class FrequencyDomainData():
    '''
    A simple class for THz-TDS, it calculates the fft of the measured time domain trace and provides
    useful functions on the fft
    * a uncertainty calculation is also carried out
    '''
    
    FMIN=0      #minimal kept frequency
    FMAX=-1     #maximal kept frequency
    
    def fromTimeDomainData(tdd,freqresolution=0):
        '''creates a FrequencyDomainData object from a timedomaindata
        '''
        if freqresolution<=0:
            N=tdd.getSamplingPoints()
        else:
            N=int(1/tdd.getTimeStep()/freqresolution)
        frequencies=np.fft.rfftfreq(N,tdd.getTimeStep())
        spectrum=np.fft.rfft(tdd.getEfield(),N)
        phase=np.unwrap(np.angle(spectrum)) #
        return FrequencyDomainData(frequencies,spectrum,phase)
    
    def divideTwoSpectra(fdNumerator,fdDenominator):
        '''Calculate the Transferfunction fdNumerator/fdDenominator
            In our context, most of the time fdNumerator=Sample Spectrum and fdDenominator=ReferenceSpectrum
        '''
        if not np.all((abs(fdNumerator.getFrequenciesRef()-fdDenominator.getFrequenciesRef())<fdNumerator.getfbins()/1000):
            print("Frequency axis of the two inputs are not equal, try to fix")            
            return 0
        
        spectrum=fdNumerator.getSpectrum()/fdDenominator.getSpectrum()
        phase=fdNumerator.getPhases()-fdDenominator.getPhases()
        return FrequencyDomainData(fdNumerator.getFrequencies(),spectrum,phase)
    
    def __init__(self,frequencies,spectrum,phase=None):

        self.frequencies=np.copy(frequencies)
        self.spectrum=np.copy(spectrum)
        #the phase should always be kept        
        if phase is not None:
            self.phase=np.copy(phase)
        else:
            self.phase=np.unwrap(np.angle(self.spectrum))
    
    def plotme(self,plotstyle=1):
        if plotstyle==2:
            plt.subplot(2,1,1)
            plt.plot(self.frequencies/1e12,self.getNormalizedSpectrum())
            plt.xlim(0,7)
            plt.ylim(-90,0)
            plt.xlabel('Frequency (THz)')
            plt.ylabel('Normalized Power (dB)')
            plt.subplot(2,1,2)
            plt.plot(self.frequencies/1e12,self.getPhasesRef())
            plt.xlim(0,7)
            
            plt.xlabel('Frequency (THz)')
            plt.ylabel('Phase')
        else:
            plt.plot(self.frequencies/1e12,self.getNormalizedSpectrum())
            plt.xlim(0,7)
            plt.ylim(-90,0)
            plt.xlabel('Frequency (THz)')
            plt.ylabel('Normalized Power (dB)')
        
    def getFrequencies(self):
        return np.copy(self.frequencies)
    
    def getFrequenciesRef(self):
        return self.frequencies
        
    def getSpectrum(self):
        return np.copy(self.spectrum)
    
    def getNormalizedSpectrum(self):
        return 20*np.log10(abs(self.spectrum)/max(abs(self.spectrum)))
    
    def getSpectrumRef(self):
        return self.spectrum
        
    def getPhases(self):
        return np.copy(self.phase)
        
    def getPhasesRef(self):
        return self.phase     
        
    def getMaxFreq(self):
        '''returns the Frequency with the highest Power'''
        return self.getFrequenciesRef()[np.argmax(abs(self.spectrum))]
        
    def getfbins(self):
        #the frequency spacing
        return abs(self.getFrequenciesRef()[1]-self.getFrequenciesRef()[0])

    def getCroppedData(self,startfreq=FMIN,endfreq=FMAX):
        #this function returns the array data cropped from startfreq to endfreq
        ix=np.all([self.getFrequenciesRef()>=startfreq,self.getFrequenciesRef()<=endfreq],axis=0)
        return FrequencyDomainData(self.getFrequenciesRef()[ix],self.getSpectrumRef()[ix],self.getPhasesRef()[ix])

    def removePhaseOffset(self,startfreq=200e9,endfreq=1e12):
        '''interpolates the phase linearly inbetween startfreq and endfreq
        removes afterwards the constant phase shift (zero phase at zero frequency)
        and returns a new FrequencyDomain instance with the new phase
        '''
        #cut phase to reasonable range:
        cropped_data=self.getCroppedData(startfreq,endfreq)
        #determine the slope and the offset      
        p=np.polyfit(cropped_data.getFrequenciesRef(),cropped_data.getPhasesRef(),1)
        #return full phase-offset(croppedPhase)
        
        return FrequencyDomainData(self.getFrequenciesRef(),self.getSpectrumRef(),self.getPhasesRef()-p[1])
    
#    def getFilteredData(self,windowlen=100e9,maxorder=3):
#        #windowlength should be smaller the etalon frequency, in order to emphasize them
#        #for the length algorithm
#        #it should be of the order of typical absorption line widthes for their analysis
#        #remove etalon minima
#        
#        #what happens if the savgol filter is applied to real and imaginary part of the fft?!        
#        #maybe nice to check?!
#        fbins=self.getfbins()
#        #100e9 should be replaced by the Etalon frequency 
#        N_min=int(windowlen/fbins)
#        order=min(N_min-1,maxorder)
#        
#        absdata=signal.savgol_filter(self.getFAbs(),N_min-N_min%2+1,order)
#        phdata=signal.savgol_filter(self.getFAbs(),N_min-N_min%2+1,order)
#       
#        return py.column_stack((self.fdData[:,:3],absdata,phdata,self.fdData[:,5:]))

    def getInterpolatedFDData(self,newfbins,strmode='linear'):
        
        oldfreqs=self.getFrequenciesRef()
        newfreqs=np.arange(min(oldfreqs),max(oldfreqs),newfbins)
                
        inter=interp1d(oldfreqs,np.asarray([self.getSpectrumRef(),self.getPhasesRef()]),strmode,axis=1)
        
        return FrequencyDomainData(newfreqs,inter(newfreqs)[0],inter(newfreqs)[1].real)

    def getSamplingPoints(self):
        #the length of the fdData array
        return self.getFrequenciesRef().shape[0]
    
    def getMovingAveragedData(self,window_size_GHz=-1):
        
        if window_size_GHz<0.5e9:
            #window_size=int(self.getEtalonSpacing()/self.getfbins())
            window_size=int(0.5e9/self.getfbins())
        else:
            window_size=int(window_size_GHz/self.getfbins())
              
        window_size+=window_size%2+1
        window=np.ones(int(window_size))/float(window_size)
        
        spectrum=np.convolve(self.getSpectrumRef(), window, 'valid')
        phase=np.convolve(self.getPhasesRef(), window, 'valid')
        one=np.ones((window_size-1)/2,)
        spectrum=np.concatenate((spectrum[0]*one,spectrum,spectrum[-1]*one))
        phase=np.concatenate((phase[0]*one,phase,phase[-1]*one))
        return FrequencyDomainData(self.getFrequenciesRef(),spectrum,phase)
    
    def getBandwidth(self):
        '''returns the bandwidth , evaluated, using the phase'''
        #average for 100 GHz
        N=int(100e9/self.getfbins())+1
        freqs=self.getFrequenciesRef()
        ix_f=freqs>200e9
        freqs=freqs[ix_f]
        phase=np.convolve(self.getPhases()[ix_f], np.ones((N,))/N)[(N-1):]
        jumps=np.where(np.diff(np.sign(np.diff(phase))))[0]
        freqs=freqs[jumps]
        return np.mean(freqs[0])
     

    def getSNR(self):
        #returns the signal to noise ratio fast method
        norm=abs(self.getSpectrumRef())/np.amax(abs(self.getSpectrumRef()))
        freqs=self.getFrequenciesRef()
        ix_freq=np.all([freqs>=self.getBandwidth(),freqs<8e12],axis=0)
        SNR=norm/np.mean(abs(norm[ix_freq]))
        
        return 20*np.log10(SNR)

    def getDynamicRange(self):
        '''Returns the Dynamic Range'''
#        #this should be the noiselevel of the fft
#        noiselevel=py.sqrt(py.mean(abs(py.fft(self._tdData.getAllPrecNoise()[0]))**2))
#        #apply a moving average filter on log
#        window_size=5
#        window=py.ones(int(window_size))/float(window_size)
#        hlog=py.convolve(20*py.log10(self.getFAbs()), window, 'valid')
#        one=py.ones((2,))
#        hlog=py.concatenate((hlog[0]*one,hlog,hlog[-1]*one))
        return -20*np.ones((self.getFrequenciesRef().shape))


    def getEtalonSpacing(self,fmin=200e9,fmax=1e12):
        '''this should return the frequency of the Etalon
        '''
    
        #how to find a stable range!
        rdata=self.getCroppedData(fmin,fmax)  
        
        #need to interpolate data!
        oldfreqs=rdata.getFrequenciesRef()
        intpdata=interp1d(oldfreqs,abs(rdata.getSpectrumRef()),'cubic')
        
        fnew=np.arange(min(oldfreqs),max(oldfreqs),0.1e9)
        absnew=intpdata(fnew)
        #find minimia and maxima        
        ixmaxima=signal.argrelmax(absnew)[0]
        ixminima=signal.argrelmin(absnew)[0]
        
        fmaxima=np.mean(np.diff(fnew[ixmaxima]))
        fminima=np.mean(np.diff(fnew[ixminima]))
        #calculate etalon frequencies
        df=(fmaxima+fminima)*0.5 #the etalon frequencies
        print(str(df/1e9) + " GHz estimated etalon frequency")
        return df



 
    
#    def calculateFDunc(self):
#        #Calculates the uncertainty of the FFT according to:
#        #   - J. M. Fornies-Marquina, J. Letosa, M. Garcia-Garcia, J. M. Artacho, "Error Propagation for the transformation of time domain into frequency domain", IEEE Trans. Magn, Vol. 33, No. 2, March 1997, pp. 1456-1459
#        #return asarray _tdData
#        #Assumes tha the amplitude of each time sample is statistically independent from the amplitude of the other time
#        #samples
#
#        # Calculates uncertainty of the real and imaginary part of the FFT and ther covariance
#        unc_E_real = []
#        unc_E_imag = []
#        cov = []
#        for f in self.getfreqs():
#            unc_E_real.append(py.sum((py.cos(2*py.pi*f*self._tdData.getTimes())*self._tdData.getUncEX())**2))
#            unc_E_imag.append(py.sum((py.sin(2*py.pi*f*self._tdData.getTimes())*self._tdData.getUncEX())**2))
#            cov.append(-0.5*sum(py.sin(4*py.pi*f*self._tdData.getTimes())*self._tdData.getUncEX()**2))
#        
#        unc_E_real = py.sqrt(py.asarray(unc_E_real))
#        unc_E_imag = py.sqrt(py.asarray(unc_E_imag))
#        cov = py.asarray(cov)
#        
#        # Calculates the uncertainty of the modulus and phase of the FFT
#        unc_E_abs = py.sqrt((self.getFReal()**2*unc_E_real**2+self.getFImag()**2*unc_E_imag**2+2*self.getFReal()*self.getFImag()*cov)/self.getFAbs()**2)
#        unc_E_ph = py.sqrt((self.getFImag()**2*unc_E_real**2+self.getFReal()**2*unc_E_imag**2-2*self.getFReal()*self.getFImag()*cov)/self.getFAbs()**4)
#        
#        t=py.column_stack((self.getfreqs(),unc_E_real,unc_E_imag,unc_E_abs,unc_E_ph))
#        return self.getcroppedData(t)  

#    def findAbsorptionLines(self):
#        #no interpolation here, just return the measured data peak index
#        tresh=6 #db treshhold to detect a line
#        #filter first
#    
#        #bring to logarithmic scale
#        movav=self.getmovingAveragedData()        
#        hlog=20*py.log10(movav[:,2])
#        #smooth data
#        #remove etalon minima
#        ixlines_prob=signal.argrelmin(hlog)[0]
#        ixlines=[]
#        ixmax=signal.argrelmax(hlog)[0]        
#        #check depth of minima compared two the two adjacent maxima
#        
#        if ixmax[0]>ixlines_prob[0]:
#            #true if starts with minimum, treat this case separatedly
#            if hlog[ixmax[0]]-hlog[ixlines_prob[0]]>tresh:
#                ixlines.append(ixlines_prob[0])
#                #and remove first minimum from ixlines_prob
#            ixlines_prob=py.delete(ixlines_prob,0)
#                
#        if ixmax[-1]<ixlines_prob[-1]: #treat this case separatedly
#            if hlog[ixmax[-1]]-hlog[ixlines_prob[-1]]>tresh:
#                ixlines.append(ixlines_prob[0])
#            ixlines_prob=py.delete(ixlines_prob,-1)
#
#        #now the remaining elements of ixlines, ixmax should be in the order of max min max ... min max         
#        leftdist=hlog[ixmax[:-1]]-hlog[ixlines_prob]
#        rightdist=hlog[ixmax[1:]]-hlog[ixlines_prob]
#        
#        for i in range(len(ixlines_prob)):
#            #check if distance is higher than treshhold
#            if (leftdist[i]>tresh or rightdist[i]>tresh):
#                ixlines.append(ixlines_prob[i])
#        
#        #remove inappropriate lines (by distance so far, at least 5 datapoints!
#        ixlines=py.asarray(ixlines)
#        ixlines=py.concatenate(([ixlines[0]],ixlines[py.diff(ixlines)>5]))
#        return ixlines
#    
#    def findAbsorptionPeaks_TESTING(self):
#        #this method might lead to better results than the findAbsorptionLines method?
#        #movav=self.getmovingAveragedData()        
#        hlog=-20*py.log10(self.getFAbs())
#        etalon=self.getEtalonSpacing()
#        Ns=int(etalon/self.getfbins())
#        Ns=py.arange(max(1,Ns-10),Ns+10,1)        
#        peaks=signal.find_peaks_cwt((hlog),Ns)
#        return peaks
#        
#    def fitAbsorptionLines(self):
#        #this method should fit the absorption lines and return fwhm, depth, unfinished
#        peaks=self.findAbsorptionLines()
#        #naive idea: use always 1/n of the space between two peaks for cutting:
#        peakfreqs=self.getfreqs()[peaks]
#        #first and last one are special
#        #first one:
#        if len(peaks)==0:
#            #nothing to do 
#            return 0
#        
#        #if len(peaks)==1:
#            #there is only one
#            
#        #if len(peaks)==2:
#            
#        #if len(peaks)>2:
#        #(peakfreqs[1]-peakfreqs[0])

        
if __name__=="__main__":
    fns=glob.glob('Reference*.txt')
    tdData=importMarburgData(fns)
    tdData=tdData.getWindowedData(5e-12)
    fdData=FrequencyDomainData.fromTimeDomainData(tdData)
    plt.figure(1)
    tdData.plotme()
    plt.figure(2)
    fdData.plotme()
    plt.plot([0,10],[-np.amax(fdData.getSNR()),-np.amax(fdData.getSNR())])
    plt.plot([fdData.getBandwidth()/1e12,fdData.getBandwidth()/1e12],[-100,0])
    
    #plt.figure(3)    
    plt.plot(fdData.getFrequencies()/1e12,fdData.getPhases()/10)
    plt.xlim(0,7)
