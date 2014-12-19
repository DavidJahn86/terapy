import pylab as py
import glob
import os
from uncertainties import unumpy
from scipy.interpolate import interp1d
from scipy.constants import c
from scipy.optimize import minimize
from TeraData import *

class HMeas(FdData):
    '''Transfer function
    HMeas is a measured transferfunction, in inherits FdData, because it is only a special
    for of frequency domain data,
    some functions need to get overriden (i.e. calculateSTDunc,...)
    '''
    def __init__(self,FDref,FDsam,disableCut=False):
        #initialize the transferfunction with a reference and a sample measurement, 
        #both, FDref, FDsam are FdData objects
        #disableCut is normaly false, this means that H is calculated only inside the bandwidth
        self.fdref=FDref
        self.fdsam=FDsam

        #sam td_data is more reliable for example noise calculation
        self._tdData=FDsam.getassTDData()
        
        if not disableCut:
            #standard: cut the fdData inbetween the trustable frequency region
            self.resetfdData()
        else:
            #use all frequencies
            self.fdData=self.calculatefdData()
        self.maxDR=max(self.getDR())
    
    def calculatefdData(self):
        #this function overrides the original fdData (that takes normally a tdData and does the fft)
        #here instead we use the fdref and fdsam objects to calculate H
        
        #check if it is possible to divide the Fdsam and Fdref objects
        prob_str=self._checkDataIntegrity()
        if not prob_str=='good':
            print("interpolation required")
            print(prob_str)
            self._commonFreqSamRef(prob_str)
        
        #calculate the uncertainty of H_real and H_imag
        H_unc=self.calculateSTDunc()
        #phase
        H_ph=self.fdref.getFPh()-self.fdsam.getFPh()
        #H itself
        H=(self.fdsam.getFReal()+1j*self.fdsam.getFImag())/(self.fdref.getFReal()+1j*self.fdref.getFImag())
        H=py.column_stack((self.fdref.getfreqs(),H.real,H.imag,abs(H),H_ph,H_unc))
        return H    

    def calculateSTDunc(self):  
        #this function should give the same result as using uncertainty package with the previously
        #correctly calculated uncertainties of reference and sample measurements
        #x=[a,b,c,d]
        dfda=lambda x: x[2]/(x[2]**2+x[3]**2)
        dfdb=lambda x: x[3]/(x[2]**2+x[3]**2)
        dfdc=lambda x: (x[0]*x[3]**2-x[0]*x[2]**2-2*x[1]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
        dfdd=lambda x: (x[1]*x[2]**2-x[1]*x[3]**2-2*x[0]*x[1]*x[2])/(x[2]**2+x[3]**2)**2
        
        dgda=lambda x: -x[3]/(x[2]**2+x[3]**2)
        dgdb=lambda x: x[2]/(x[2]**2+x[3]**2)
        dgdc=lambda x: -(x[1]*x[2]**2-x[1]*x[3]**2-2*x[0]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
        dgdd=lambda x: (x[0]*x[3]**2-x[0]*x[2]**2-2*x[1]*x[2]*x[3])/(x[2]**2+x[3]**2)**2

        
        ta=self.fdsam.getFReal()
        tb=self.fdsam.getFImag()
        tc=self.fdref.getFReal()
        td=self.fdref.getFImag()
        
        H_unc_real2=(self.fdsam.getFRealUnc()*dfda(py.asarray([ta,tb,tc,td])))**2
        H_unc_real2+=(self.fdsam.getFImagUnc()*dfdb(py.asarray([ta,tb,tc,td])))**2
        H_unc_real2+=(self.fdref.getFRealUnc()*dfdc(py.asarray([ta,tb,tc,td])))**2
        H_unc_real2+=(self.fdref.getFImagUnc()*dfdd(py.asarray([ta,tb,tc,td])))**2
        
        H_unc_imag2=(self.fdsam.getFRealUnc()*dgda(py.asarray([ta,tb,tc,td])))**2
        H_unc_imag2+=(self.fdsam.getFImagUnc()*dgdb(py.asarray([ta,tb,tc,td])))**2
        H_unc_imag2+=(self.fdref.getFRealUnc()*dgdc(py.asarray([ta,tb,tc,td])))**2
        H_unc_imag2+=(self.fdref.getFImagUnc()*dgdd(py.asarray([ta,tb,tc,td])))**2
        absEsam=unumpy.uarray(self.fdsam.getFAbs(),self.fdsam.getFAbsUnc())
        Ephsam=unumpy.uarray(self.fdsam.getFPh(),self.fdsam.getFPhUnc())
        absEref=unumpy.uarray(self.fdref.getFAbs(),self.fdref.getFAbsUnc())
        Ephref=unumpy.uarray(self.fdref.getFPh(),self.fdref.getFPhUnc())
  
        H_uncabs=unumpy.std_devs(absEsam/absEref)
        H_uncph=unumpy.std_devs(Ephsam-Ephref)
        return py.column_stack((py.sqrt(H_unc_real2),py.sqrt(H_unc_imag2),H_uncabs,H_uncph))
    
    def doPlot(self):
        #do the H-plots
        freqs=self.getfreqsGHz()
        
        py.figure('H-UNC-Plot')
        py.plot(freqs,self.getFReal())
        py.plot(freqs,self.getFImag())
        py.plot(freqs,self.getFReal()+self.getFRealUnc(),'g--',freqs,self.getFReal()-self.getFRealUnc(),'g--')
        py.plot(freqs,self.getFImag()+self.getFImagUnc(),'g--',freqs,self.getFImag()-self.getFImagUnc(),'g--')
        py.xlabel('Frequency in GHz')
        py.ylabel('Transfer Function')
        py.legend(('H_real','H_imag'))
#        
        py.figure('H-PHASE-Plot')
        py.plot(freqs,self.getFPh())
       
        py.figure('H-ABS-Plot')
        py.plot(freqs,self.getFAbs())
        
    def estimateLDavid(self):
        #some method to estimate the thickness from phase slope + Etalon frequency
        #crop frequency axis
        rdata=self.getcroppedData(self.fdData,200e9,1e12)
        #calculate phase change
        p=py.polyfit(rdata[:,0],rdata[:,4],1)        
        kappa=abs(p[0])
#        py.figure()
#        py.plot(rdata[:,0],rdata[:,-1])
        df=self.getEtalonSpacing()
        #assume air around
        n=1.0/(1.00027-kappa*df/py.pi)
        l=c/2.0*(1.00027/df-kappa/py.pi)
        return [l,n]

    def getDR(self):
        
         return self.fdsam.getDR()
            
    def manipulateFDData(self,fbins,fbnds,mode='interpolate'):
        #this method provides the means to change the underlying fdData objects and recalculate H
        #if just an interpolated H with fbins is needed, use getInterpolatedFdData from FdData class        
        #maybe I will delete it soon, because it shouldn't be used!       
        if fbins>0.5e9:
            if mode=='zeropadd':
                self.zeroPadd(fbins)
            else:
                self.fdref.setFDData(self.fdref.getInterpolatedFDData(fbins))
                self.fdsam.setFDData(self.fdsam.getInterpolatedFDData(fbins))
               
        self.fdref.setFDData(self.fdref.getcroppedData(self.fdref.fdData,fbnds[0],fbnds[1]))        
        self.fdsam.setFDData(self.fdsam.getcroppedData(self.fdsam.fdData,fbnds[0],fbnds[1]))
        
        self.fdData=self.calculatefdData()

    def resetfdData(self,fbins=-1,fbnds=[FdData.FMIN,FdData.FMAX]):
        #restricts the H to fbnds, and also zeropadds to fbins
        minrf,maxrf=self.fdref.getBandwidth()
        minsf,maxsf=self.fdsam.getBandwidth()
        self.manipulateFDData(fbins,[max(minrf,minsf,FdData.FMIN),min(maxrf,maxsf,FdData.FMAX)])
        
    def zeroPadd(self,fbins):
        #zeropadding of H
        self.fdref.zeroPadd(fbins)
        self.fdsam.zeroPadd(fbins)
        self.calculatefdData()

    def _checkDataIntegrity(self):

        tdrefData=self.fdref.getassTDData()
        tdsamData=self.fdsam.getassTDData()

        #begin and end of timeaxis should differ not more than 5 as
        if abs(tdrefData.tdData[0,0]-tdsamData.tdData[0,0])>5e-15:
            return 'td-Problem-phase'

        #length of time axis should be identical
        if tdrefData.getLength()!=tdsamData.getLength():
            return 'td-Problem-len'
        
        
        if abs(tdrefData.tdData[-1,0]-tdsamData.tdData[-1,0])>1e-16 or\
            abs(tdrefData.tdData[0,0]-tdsamData.tdData[0,0])>1e-16:
            return 'td-Problem-interval'
        
        #length in frequency domain should be equal
        if self.fdref.getLength()!= self.fdsam.getLength():
            return 'fd-Problem-len'
            
        #frequency bins should be equal
        if not all(self.fdref.getfreqs()-self.fdsam.getfreqs()<1e6):
            return 'fd-Problem-axis'
         #else return       
        return 'good'
      
    def _commonFreqSamRef(self,prob_str):
        #take care, this actually manipulates the real underlying data!
        #maybe consider to do a deepcopy!
        tdref=self.fdref.getassTDData()
        tdsam=self.fdsam.getassTDData()
  
        #if the absolute time, at which one of the two measurements starts differs significantly
        #from the other, we should put out a warning, and put trailing zeros
        #after that we are ready for interpolating the data
        if prob_str=='td-Problem-phase':
            t_start_ref=tdref.tdData[0,0]
            t_start_sam=tdsam.tdData[0,0]
            N=int((t_start_sam-t_start_ref)/tdsam.dt)
            tdsam.zeroPaddData(N,'zero','start')
      
        #(max min notation for readability?
        cmin=max(min(tdref.getTimes()),min(tdsam.getTimes()))
        cmax=min(max(tdref.getTimes()),max(tdsam.getTimes()))
        clen=min(tdref.getLength(),tdsam.getLength())
        
        #safe also old bnds        
        minrf,maxrf=self.fdref.getBandwidth()
        minsf,maxsf=self.fdsam.getBandwidth()
                
        tdrefnew=THzTdData(tdref.getInterData(tdref.tdData,clen,cmin,cmax),tdref.getfilename(),tdref._thzdata_raw,existing=True)
        tdsamnew=THzTdData(tdsam.getInterData(tdsam.tdData,clen,cmin,cmax),tdsam.getfilename(),tdsam._thzdata_raw,existing=True)
        
        self.fdref=FdData(tdrefnew,-1,[min(minrf,minsf),max(maxrf,maxsf)])
        self.fdsam=FdData(tdsamnew,-1,[min(minrf,minsf),max(maxrf,maxsf)])  
        
class teralyz():
    '''Finds the optical constants
    this class implements the solver, i.e. the length finding algorithm on the measurement
    data
    '''
    #refractive index of air
    n_0=1.00027-0.0000j

    def __init__(self,measurementdata,thickness=None):
        #measurementdata should be an object of type HMeas
        self.H=measurementdata
        #H of only the first pulse        
        self.H_firstPuls=self.getHFirstPuls()
        
        #just to have a first estimate
        l,n=self.H.estimateLDavid()        
        self.n_estimated=n
        self.l_estimated=l

        #its always better to provide a thickness    
        if thickness==None:
            self.userthickness=self.l_estimated
        else:                 
            self.userthickness=thickness
       
        #calculate an estimated n from the time domain shift of the pulses
        self.n_estimated_user=self.nestimatedTD(measurementdata.fdsam._tdData,measurementdata.fdref._tdData)
        #use as the "optimal" thickness the thickness provided by the user, as a first guess
        self.l_opt=self.userthickness
        #set the found n to zero
        self.n=0
        #output some information
        print("A priori estimated Width: " + str(l*1e6))
        print("A priori estimated n: " + str(n))
        #calculate the number of expected fabry-pero echos
        self.no_echos=self.calculate_no_echos(measurementdata.fdsam._tdData)
        #get a suggestion for the filename
        self.outfilename=self.getFilenameSuggestion()
 
    def calculate_no_echos(self,tdsam):
        #use l and n estimated to estimated roughly the number of fabry-pero pulses
        
        #time window
        t_max=tdsam.getTimeWindowLength()
        
        delta=0.5*(t_max*c/self.l_estimated/self.n_estimated)

        return int(delta)
  
    def calculateinits(self,H,l):
        #this uses the approximation that imaginary part of n is small, and solves
        #the equation in this approximation analytically. Fabry-perot reflections are cut
        #take care that the phase has ben offset_removed!
    
        #the frequency axis that will be used for the calculation
        oldfreqaxis=self.H.getfreqs()
        #this could be done once for all to gain speed!
        intpH=interp1d(self.H_firstPuls.getfreqs(),self.H_firstPuls.fdData[:,1:],axis=0)
        newH=intpH(oldfreqaxis)

        #crop the time domain data to the first pulse and apply than the calculation
        n=self.n_0.real-newH[:,3]/(2*py.pi*oldfreqaxis*l)*c
        alpha=-c/(2*py.pi*oldfreqaxis*l)*py.log(newH[:,2]*(n+1)**2/(4*n))

        return n,alpha
   
    def calculaten(self,H,l):
        #this calculates n for a given transferfunction H and a given thickness l
    
        #calculate the initial values
        inits=py.asarray(self.calculateinits(H,l))
        res=[] #the array of refractive index
        vals=[] # the array of error function values
        bnds=((1,None),(0,None)) #not used at the moment, with SLSQP bnds can be introduced
        nums=len(H[:,0])
        #minimize the deviation of H_theory to H_measured for each frequency
        for i in range(nums):
#            t=minimize(self.error_func,[inits[0,i],inits[1,i]],args=(H[i,:2],l), method='SLSQP',\
#            bounds=bnds, options={'ftol':1e-9,'maxiter':2000, 'disp': False})
            t=minimize(self.error_func,[inits[0,i],inits[1,i]],args=(H[i,:3],l), method='Nelder-Mead',
            options={'xtol': 1e-6,'disp':False})
            res.append(t.x[0]-1j*t.x[1])
            vals.append(t.fun)
        n=py.asarray(res)
        #self.n is a 5xlengthf array, frequency,n_real,n_imag,n_smoothed_real,n_smoothed_imag
        self.n=py.column_stack((H[:,0],n,n))
        return n

        
    def doCalculation(self,bool_findl=1,n_SVMAFS=5,bool_silent=0):
        #this function does the complete calculation
        #bandwidth
        bw=self.H.getBandwidth()
        #crop data
        self.H.manipulateFDData(-1,[bw[0]+50e9,bw[1]-200e9])
        #if l_opt shouldn't be calculated used bool_findl=False
        if bool_findl:
            self.l_opt=self.findLintelli()

        print('\033[92m\033[1m' + '  Use Sample Thickness: ' + str(self.l_opt*1e6) + ' micro m ' + '\033[0m')

        #calculate n for the given l_opt
        n=self.calculaten(self.H.fdData,self.l_opt)
        n_smoothed=n
        i=0
        
        #smooth it with SVMAF
        while i<n_SVMAFS:
            n_smoothed=self.SVMAF(self.H.getfreqs(),n_smoothed,self.l_opt)
            i+=1

        self.n=py.column_stack((self.H.getfreqs(),n,n_smoothed))        
        
        return self.n
  
    def errorL(self,l,H):
        #the error function for the length finding minimization
        
        #calculate n for the short transfer function H and the length l
        n_small=[self.calculaten(H,l)]
        #evaluate the quasi space value
        qs=self.QuasiSpace(n_small,H[1,0]-H[0,0],l)
        #evaluate the total variation value
        tv=self.totalVariation(n_small)
        #plot them
        py.plot(l,qs[0],'+')
        py.plot(l,tv[0],'*')        
        print("Currently evaluating length: "+ str(l[0]*1e6) + " TV Value " + str(tv[0]))
        return qs[0]
    
    def error_func(self,n,H,l):
        #the quadratic deviation of Htheroy and Hmeasuered, used for finding n
        H_t=self.H_theory(H[0],n,l)
        return (H_t.real-H[1])**2+(H_t.imag-H[2])**2
    
    def findLintelli(self):
        #find the frequency with the most amplitude
        fmax=self.H.fdref.getmaxfreq()
        
        #overthink once more the cropping length
        n=self.n_estimated
        #oscillation period can be calculated from H data!
        f_span=c/2*1/self.l_estimated*1/n #(this should be the length of the oscillation)
        #restrict H to only 5 oscillations
        H_small=self.H.getcroppedData(self.H.fdData,fmax-f_span*1,fmax+f_span*4)
        py.figure(33)
        #minimize quasispace/totalvariation value
        t=minimize(self.errorL,self.userthickness,args=((py.asarray(H_small),)),\
        method='Nelder-Mead', options={'xtol':1e-6,'disp': False})#, 'disp': False})
        return t.x[0]

    def getHFirstPuls(self):
        #returns the transferfunction of the timedomain data that corresponds to the first
        #pulse only
        origlen=self.H.fdref._tdData.getLength() 
        refdata=THzTdData(self.H.fdref._tdData.getFirstPuls(10e-12),existing=True)
        samdata=THzTdData(self.H.fdsam._tdData.getFirstPuls(5e-12),existing=True)
        #bring it to the old length
        refdata.zeroPaddData(origlen-refdata.getLength())
        samdata.zeroPaddData(origlen-samdata.getLength())
        #calculate the fouriertransform
        firstref=FdData(refdata)
        firstsam=FdData(samdata)
        #calculate the transferfunction for the first puls
        H_firstPuls=HMeas(firstref,firstsam,disableCut=True)
        
        return H_firstPuls

    def getFilenameSuggestion(self):
        #tries to give a valid filename from the sample filenames
        filenames=self.H.fdsam._tdData.getfilename()
        filenames=filenames[0]
        filenames=os.path.split(filenames)[-1]
             
        if '_' in filenames:
            #try to split along _
            filenames=filenames.split('_')
        elif '-' in filenames:
            #if not found try to use -
            filenames=filenames.splot('-')
        elif ' ' in filenames:
            filenames=filenames.split(' ')
        elif '.' in filenames:
            filenames=filenames.split('.')

        if len(filenames)==2:
            filenames=filenames[0]
        else:
            filenames=filenames[0]+'_'+filenames[1]+'_'
        
        return filenames

    def H_theory(self,freq,n,l):
        #calculate the theoretic transfer function, so far only one medium!
        nc=n[0]-1j*n[1]
        r_as=(self.n_0-nc)/(self.n_0+nc)
       
        FPE=0
        P=lambda tn: py.exp(-1j*l*freq*2*py.pi*tn/c)
        for sumk in range(1,self.no_echos):
            FPE+=((r_as**2)*(P(nc)**2))**sumk
        
        H=4*self.n_0*nc/(nc+self.n_0)**2*P((nc-self.n_0))*(1+FPE)
        return H

    def nestimatedTD(self,tdsam,tdref):
        #the time delay between sample and reference pulse is used to estimate n
        tmaxs=tdsam.getPeakPosition()
        tmaxr=tdref.getPeakPosition()
        n=1+c/self.userthickness*abs(tmaxs-tmaxr)
      
        return n

    def plotRefractiveIndex(self,bool_plotsmoothed=1,savefig=0,filename=None):
        #plot the refractive index
    
        if filename==None:
            figname_b=self.getFilenameSuggestion()
        else:
            figname_b=filename
            
        py.figure('Refractive_Real_Plot')
        if bool_plotsmoothed==1:
            py.plot(self.n[:,0].real/1e9,self.n[:,2].real)
        else:
            py.plot(self.n[:,0].real/1e9,self.n[:,1].real)
        
        py.xlabel('Frequency in GHz')
        py.ylabel('n Real')
        if savefig:
            figname=figname_b+'n_real.png'
            py.savefig(figname,dpi=200)
            
        py.figure('Refractive_Imag_Plot')
        if bool_plotsmoothed==1:
            py.plot(self.n[:,0].real/1e9,-self.n[:,2].imag)
        else:
            py.plot(self.n[:,0].real/1e9,-self.n[:,1].imag)

        py.xlabel('Frequency in GHz')
        py.ylabel('n Imag')
        if savefig:
            figname=figname_b+'n_imag.png'
            py.savefig(figname,dpi=200)
    
    def plotErrorFunction(self,l,freq):
        #plots the error function
        py.figure()
        resolution=300
        ix=py.argmin(abs(freq-self.H.getfreqs()))
        n_i=py.linspace(0.0,0.05,int(resolution/2))
        n_r=py.linspace(1,3,resolution)
        N_R,N_I=py.meshgrid(n_r,n_i)
        E_fu=py.zeros((len(n_i),len(n_r)))
        for i in range(len(n_r)):
            for k in range(len(n_i)):
                E_fu[k,i]=self.error_func([n_r[i],n_i[k]],self.H.fdData[ix,:3],l)
#        print(E_fu[:,2])
        py.pcolor(N_R,N_I,py.log10(E_fu))
        py.colorbar()

    def plotInits(self,H,l,figurenumber=200):
        #plots the initial conditions
        inits=self.calculateinits(H,l)
        py.figure(figurenumber)
        py.title('Initial Conditions')
        py.plot(self.H.getfreqsGHz(),inits[0])
        py.plot(self.H.getfreqsGHz(),-inits[1])
        py.xlabel('Frequency in GHz')
        py.ylabel('optical constant value')
        py.legend(('n_real','n_imag'))
    
    def QuasiSpace(self,ns,df,ls):
        #evaluates the quasi space value
        allqs=[]
        #at the moment, the most easy method(everything else failed...)
        #get the mean absolute value of the complete spectrum. works best
        for i in range(len(ns)):
#            xvalues=py.fftfreq(len(ns[i]),df)*c/4

            QSr=py.fft(ns[i].real-py.mean(ns[i].real))
            QSi=py.fft(ns[i].imag-py.mean(ns[i].real))
            #naive cut:
            ix=range(3,int(len(QSr)/2-1))
        
            QSr=QSr[ix]
            QSi=QSi[ix]
            allqs.append(py.mean(abs(QSr))+py.mean(abs(QSi)))
        
        return allqs

    def saveResults(self,filename=None):
        #save the results to a file        
        H_theory=self.H_theory(self.H.getfreqs(),[self.n[:,1].real,self.n[:,1].imag],self.l_opt)        
        #built the variable that should be saved:        
        savetofile=py.column_stack((
        self.H.getfreqs(), #frequencies
        self.n[:,1].real,-self.n[:,1].imag, #the real and imaginary part of n
        self.n[:,2].real,-self.n[:,2].imag, #the real and imaginary part of the smoothed n
        self.H.getFReal(),self.H.getFImag(),#H_measured
        self.H.getFRealUnc(),self.H.getFImagUnc(),#uncertainties
        self.H.getFAbs(),self.H.getFPh(),#absH,ph H measured    
        H_theory.real,H_theory.imag, #theoretic H
        ))

        headerstr=('freq, ' 
        'ref_ind_real, ref_ind_imag, '
        'ref_ind_sreal, ref_ind_simag, '
        'H_measured_real, H_measured_imag, '
        'Hunc_real, Hunc_imag, '
        'abs(H_measured), angle(H_measured), '
        'H_theory_real, H_theory_imag')
        if filename==None:
            fname=self.getFilenameSuggestion()
        else:
            fname=filename+"_"
        fname+='analyzed_' + 'D=' + str(self.l_opt/1e-6) +'.txt'
        py.savetxt(fname,savetofile,delimiter=',',header=headerstr)

    def setOutfilename(self,newname):
        #outfilename
        self.outfilename=newname

    def SVMAF(self,freq,n,l):
        #Apply the SVMAF filter to the material parameters
        runningMean=lambda x,N: py.hstack((x[:N-1],py.convolve(x,py.ones((N,))/N,mode='same')[N-1:-N+1],x[(-N+1):]))
        #calculate the moving average of 3 points
        n_smoothed=runningMean(n,3)
        #evaluate H_smoothed from n_smoothed
        H_smoothed=self.H_theory(freq,[n_smoothed.real,n_smoothed.imag],l)
        
        H_r=H_smoothed.real
        H_i=H_smoothed.imag
        f=1
        #the uncertainty margins
        lb_r=self.H.getFReal()-self.H.getFRealUnc()*f
        lb_i=self.H.getFImag()-self.H.getFImagUnc()*f
        ub_r=self.H.getFReal()+self.H.getFRealUnc()*f
        ub_i=self.H.getFImag()+self.H.getFImagUnc()*f
        
        #ix=all indices for which after smoothening n H is still inbetwen the bounds        
        ix=py.all([H_r>=lb_r,H_r<ub_r,H_i>=lb_i,H_i<ub_i],axis=0)
#        #dont have a goood idea at the moment, so manually:
        for i in range(len(n_smoothed)):
            if ix[i]==0:
                n_smoothed[i]=n[i]
        print("SVMAF changed the refractive index at " + str(sum(ix)) + " frequencies")
        return n_smoothed      

    def totalVariation(self,ns):
        #calculate the total variation value
        allvs=[]
        
        for i in range(len(ns)):
            tv1=0            
            for dm in range(1,len(ns[i])):
                tv1+=abs(ns[i][dm-1].real-ns[i][dm].real)+\
                abs(ns[i][dm-1].imag-ns[i][dm].imag)
            
            allvs.append(tv1)
            
        return allvs