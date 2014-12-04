import pylab as py
import copy
import glob
import os
from uncertainties import unumpy
from scipy.interpolate import interp1d
from scipy.constants import c
from scipy.optimize import minimize
from TeraData import *

class HMeas(FdData):
 
    def __init__(self,FDref,FDsam,disableCut=False):
        self.fdref=FDref
        self.fdsam=FDsam
        #sam td_data is more reliable for example noise calculation
        self._tdData=FDsam.getassTDData()
        
        if not disableCut:
            self.resetfdData()
        else:
            self.fdData=self.calculatefdData()
        self.maxDR=max(self.getDR())
    #functions that need to be overridden
    def getDR(self):
         return self.fdsam.getDR()
            
    def resetfdData(self,fbins=-1,fbnds=[FdData.FMIN,FdData.FMAX]):
            minrf,maxrf=self.fdref.getBandwidth()
            minsf,maxsf=self.fdsam.getBandwidth()
            
            self.manipulateFDData(fbins,[max(minrf,minsf,FdData.FMIN),min(maxrf,maxsf,FdData.FMAX)])
        
    def manipulateFDData(self,fbins,fbnds,mode='interpolate'):
        #this method provides the means to change the underlying fdData objects and recalculate H
        #if just an interpolated H with fbins is needed, use getInterpolatedFdData from FdData class        
        if fbins>0.5e9:
            if mode=='zeropadd':
                self.zeroPadd(fbins)
            else:
                self.fdref.setFDData(self.fdref.getInterpolatedFDData(fbins))
                self.fdsam.setFDData(self.fdsam.getInterpolatedFDData(fbins))
               
        self.fdref.setFDData(self.fdref.getcroppedData(self.fdref.fdData,fbnds[0],fbnds[1]))        
        self.fdsam.setFDData(self.fdsam.getcroppedData(self.fdsam.fdData,fbnds[0],fbnds[1]))
        
        self.fdData=self.calculatefdData()
   
    def zeroPadd(self,fbins):
        self.fdref.zeroPadd(fbins)
        self.fdsam.zeroPadd(fbins)
        self.calculatefdData()
    
    def calculateSTDunc(self):  
        #x=[a,b,c,d]
        dfda=lambda x: x[2]/(x[2]**2+x[3]**2)
        dfdb=lambda x: x[3]/(x[2]**2+x[3]**2)
        dfdc=lambda x: (x[0]*x[3]**2-x[0]*x[2]**2-2*x[1]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
        dfdd=lambda x: (x[1]*x[2]**2-x[1]*x[3]**2-2*x[0]*x[1]*x[2])/(x[2]**2+x[3]**2)**2
        
        dgda=lambda x: -x[3]/(x[2]**2+x[3]**2)
        dgdb=lambda x: x[2]/(x[2]**2+x[3]**2)
        dgdc=lambda x: -(x[1]*x[2]**2-x[1]*x[3]**2-2*x[0]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
        dgdd=lambda x: (x[0]*x[3]**2-x[0]*x[2]**2-2*x[1]*x[2]*x[3])/(x[2]**2+x[3]**2)**2

        
        ta=self.fdsam.fdData[:,1].real
        tb=self.fdsam.fdData[:,1].imag
        tc=self.fdref.fdData[:,1].real
        td=self.fdref.fdData[:,1].imag
        
        H_unc_real2=(self.fdsam.fdData[:,4]*dfda(py.asarray([ta,tb,tc,td])))**2
        H_unc_real2+=(self.fdsam.fdData[:,5]*dfdb(py.asarray([ta,tb,tc,td])))**2
        H_unc_real2+=(self.fdref.fdData[:,4]*dfdc(py.asarray([ta,tb,tc,td])))**2
        H_unc_real2+=(self.fdref.fdData[:,5]*dfdd(py.asarray([ta,tb,tc,td])))**2
        
        H_unc_imag2=(self.fdsam.fdData[:,4]*dgda(py.asarray([ta,tb,tc,td])))**2
        H_unc_imag2+=(self.fdsam.fdData[:,5]*dgdb(py.asarray([ta,tb,tc,td])))**2
        H_unc_imag2+=(self.fdref.fdData[:,4]*dgdc(py.asarray([ta,tb,tc,td])))**2
        H_unc_imag2+=(self.fdref.fdData[:,5]*dgdd(py.asarray([ta,tb,tc,td])))**2
        absEsam=unumpy.uarray(self.fdsam.fdData[:,2].real,self.fdsam.fdData[:,6].real)
        Ephsam=unumpy.uarray(self.fdsam.fdData[:,3].real,self.fdsam.fdData[:,7].real)
        absEref=unumpy.uarray(self.fdref.fdData[:,2].real,self.fdref.fdData[:,6].real)
        Ephref=unumpy.uarray(self.fdref.fdData[:,3].real,self.fdref.fdData[:,7].real)
  
        H_uncabs=unumpy.std_devs(absEsam/absEref)
        H_uncph=unumpy.std_devs(Ephsam-Ephref)
        return py.column_stack((py.sqrt(H_unc_real2),py.sqrt(H_unc_imag2),H_uncabs,H_uncph))
    
    def _checkDataIntegrity(self):
        tdrefData=self.fdref.getassTDData()
        tdsamData=self.fdsam.getassTDData()

        #begin and end should differ not more than 1 as
        if abs(tdrefData.tdData[0,0]-tdsamData.tdData[0,0])>5e-15:
            return 'td-Problem-phase'

        if tdrefData.getLength()!=tdsamData.getLength():
            return 'td-Problem-len'
        
        if abs(tdrefData.tdData[-1,0]-tdsamData.tdData[-1,0])>1e-16 or\
            abs(tdrefData.tdData[0,0]-tdsamData.tdData[0,0])>1e-16:
            return 'td-Problem-interval'
        
        if len(self.fdref.fdData[:,1])!= len(self.fdsam.fdData[:,1]):
            return 'fd-Problem-len'
        
        if not all(self.fdref.fdData[:,0]-self.fdsam.fdData[:,0]<1e6):
            return 'fd-Problem-axis'
         #else return       
        return 'good'
      
    def interpolateData(self,prob_str):
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
        cmin=max(min(tdref.tdData[:,0]),min(tdsam.tdData[:,0]))
        cmax=min(max(tdref.tdData[:,0]),max(tdsam.tdData[:,0]))
        clen=min(tdref.getLength(),tdsam.getLength())
        
        #safe also old bnds        
        minrf,maxrf=self.fdref.getBandwidth()
        minsf,maxsf=self.fdsam.getBandwidth()
        
        tdrefnew=THzTdData(tdref.getInterData(tdref.tdData,clen,cmin,cmax),tdref.getfilename(),tdref._thzdata_raw,existing=True)
        tdsamnew=THzTdData(tdsam.getInterData(tdsam.tdData,clen,cmin,cmax),tdsam.getfilename(),tdsam._thzdata_raw,existing=True)
        
        self.fdref=FdData(tdrefnew,-1,[min(minrf,minsf),max(maxrf,maxsf)])
        self.fdsam=FdData(tdsamnew,-1,[min(minrf,minsf),max(maxrf,maxsf)])  

    def calculatefdData(self):
        #this function overrides the original fdData (that takes normally a tdData and does the fft)
        #here instead we use the fdref and fdsam objects to calculate H
        
        prob_str=self._checkDataIntegrity()
        if not prob_str=='good':
            print("interpolation required")
            print(prob_str)
            self.interpolateData(prob_str)
        
        H_unc=self.calculateSTDunc()
            #take care that abs(H) is always smaller one!
            #try if this increases the data quality!
#            H_ph=py.unwrap(py.angle(self.fdsam.fdData[:,1]/self.fdref.fdData[:,1]))
        H_ph=self.fdref.getUnwrappedPhase()-self.fdsam.getUnwrappedPhase()
        H=self.fdsam.fdData[:,1]/self.fdref.fdData[:,1]
        H=py.column_stack((self.fdref.fdData[:,0],H,abs(H),H_ph,H_unc))
        return H    
    
    def doPlot(self):
        freqs=self.getfreqsGHz()
        
        py.figure('H-UNC-Plot')
        py.plot(freqs,self.fdData[:,1].real)
        py.plot(freqs,self.fdData[:,1].imag)
        py.plot(freqs,self.fdData[:,1].real+self.fdData[:,4].real,'g--',freqs,self.fdData[:,1].real-self.fdData[:,4].real,'g--')
        py.plot(freqs,self.fdData[:,1].imag+self.fdData[:,5].real,'g--',freqs,self.fdData[:,1].imag-self.fdData[:,5].real,'g--')
        py.xlabel('Frequency in GHz')
        py.ylabel('Transfer Function')
        py.legend(('H_real','H_imag'))
#        
        py.figure('H-PHASE-Plot')
        py.plot(freqs,self.fdData[:,3])
       
        py.figure('H-ABS-Plot')
        py.plot(freqs,self.fdData[:,2])
        
    def estimateLDavid(self):
        rdata=self.getcroppedData(self.fdData,200e9,1e12)
        #calculate phase change
        p=py.polyfit(rdata[:,0].real,rdata[:,3].real,1)        
        kappa=abs(p[0])
#        py.figure()
#        py.plot(rdata[:,0],rdata[:,-1])
        df=self.getEtalonSpacing()
        #assume air around
        n=1.0/(1.00027-kappa*df/py.pi)
        l=c/2.0*(1.00027/df-kappa/py.pi)
        return [l,n]
        
        
class teralyz():
   # n_0=1.00027#+#0.001j
    n_0=1.00027-0.0000j
   
    def __init__(self,measurementdata,thickness=None,unc=50e6,steps=10):
        
        self.H=measurementdata
        self.H_firstPuls=self.getHFirstPuls()
        
        l,n=self.H.estimateLDavid()        
        self.n_estimated=n
        self.l_estimated=l
    
        if thickness==None:
            self.userthickness=self.l_estimated
        else:                 
            self.userthickness=thickness
        
        self.userstd=unc
        self.steps=steps
        self.n_estimated_user=self.nestimatedTD(measurementdata.fdsam._tdData,measurementdata.fdref._tdData)
        self.l_opt=self.userthickness
        self.n=0
        print("A priori estimated Width: " + str(l*1e6))
        print("A priori estimated n: " + str(n))
        self.no_echos=self.calculate_no_echos(measurementdata.fdsam._tdData)
        self.outfilename=self.getFilenameSuggestion()
    
    def getHFirstPuls(self):
        origlen=self.H.fdref._tdData.getLength() 
        refdata=THzTdData(self.H.fdref._tdData.getFirstPuls(10e-12),existing=True)
        samdata=THzTdData(self.H.fdsam._tdData.getFirstPuls(5e-12),existing=True)
        refdata.zeroPaddData(origlen-refdata.getLength())
        samdata.zeroPaddData(origlen-samdata.getLength())
        
        firstref=FdData(refdata)
        firstsam=FdData(samdata)
        
        H_firstPuls=HMeas(firstref,firstsam,disableCut=True)
        
        return H_firstPuls
 
    def setOutfilename(self,newname):
        self.outfilename=newname
 
    def H_theory(self,freq,n,l):
        nc=n[0]-1j*n[1]
        r_as=(self.n_0-nc)/(self.n_0+nc)
       
        FPE=0
        P=lambda tn: py.exp(-1j*l*freq*2*py.pi*tn/c)
        for sumk in range(1,self.no_echos):
            FPE+=((r_as**2)*(P(nc)**2))**sumk
        
        H=4*self.n_0*nc/(nc+self.n_0)**2*P((nc-self.n_0))*(1+FPE)
        
        return H

    def error_func(self,n,H,l):
        #we could also introduce constraints on n by artificialy punishing        
        H_t=self.H_theory(H[0],n,l)
        H_m=H[1]

        return (H_t.real-H_m.real)**2+(H_t.imag-H_m.imag)**2
        
    def calculate_no_echos(self,tdsam):
        t_max=tdsam.getTimeWindowLength()
        
        delta=0.5*(t_max*c/self.l_estimated/self.n_estimated)

        return int(delta)
    
    def nestimatedTD(self,tdsam,tdref):
        tmaxs=tdsam.getPeakPosition()
        tmaxr=tdref.getPeakPosition()
        n=1+c/self.userthickness*abs(tmaxs-tmaxr)
      
        return n
    
    def findLintelli(self):
        fmax=self.H.fdref.getmaxfreq()
        
 #       bnds=(self.thicknessestimate-2*self.thicknessstd,self.thicknessestimate+3*self.thicknessstd)
#        print(bnds)
        #overthink once more the cropping length
        n=self.n_estimated
        #oscillation period can be calculated from H data!
        f_span=c/2*1/self.l_estimated*1/n #(this should be the length of the oscillation)
        
        calcdata=copy.deepcopy(self.H)
        calcdata.manipulateFDData(-1,[fmax-f_span*1,fmax+f_span*4])
        
        H_small=calcdata.fdData
        py.figure(33)
        t=minimize(self.errorL,self.userthickness,args=((py.asarray(H_small),)),\
        method='Nelder-Mead', options={'xtol':1e-6,'disp': False})#, 'disp': False})
        return t.x[0]
        
    def errorL(self,l,H):
        n_small=[self.calculaten(H,l)]
    
        qs=self.QuasiSpace(n_small,H[1,0]-H[0,0],l)
        tv=self.totalVariation(n_small)
        py.plot(l,qs[0],'+')
        py.plot(l,tv[0],'*')        
        print("Currently evaluating length: "+ str(l[0]*1e6) + " TV Value " + str(tv[0]))
        return qs[0]
       
    def QuasiSpace(self,ns,df,ls):
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
    
    def totalVariation(self,ns):
        allvs=[]
        
        for i in range(len(ns)):
            tv1=0            
            for dm in range(1,len(ns[i])):
                tv1+=abs(ns[i][dm-1].real-ns[i][dm].real)+\
                abs(ns[i][dm-1].imag-ns[i][dm].imag)
            
            allvs.append(tv1)
            
        return allvs
        
    def SVMAF(self,freq,n,l):
       
        #Apply the SVMAF filter to the material parameters
        #runningMean=lambda x,N: py.hstack((x[:N-1],py.cvolve(x,py.ones((N,))/N,mode='valid')[(N-1):],x[(-N+1):]))
        runningMean=lambda x,N: py.hstack((x[:N-1],py.convolve(x,py.ones((N,))/N,mode='same')[N-1:-N+1],x[(-N+1):]))
       
        n_smoothed=runningMean(n,3) # no problem in doing real and imaginary part toether
        H_smoothed=self.H_theory(freq,[n_smoothed.real,n_smoothed.imag],l)
        
        H_r=H_smoothed.real
        H_i=H_smoothed.imag
        f=1
        lb_r=self.H.fdData[:,1].real-self.H.fdData[:,4]*f
        lb_i=self.H.fdData[:,1].imag-self.H.fdData[:,5]*f
        ub_r=self.H.fdData[:,1].real+self.H.fdData[:,4]*f
        ub_i=self.H.fdData[:,1].imag+self.H.fdData[:,5]*f
        
        #ix=all indices for which after smoothening n H is still inbetwen the bounds        
        ix=py.all([H_r>=lb_r,H_r<ub_r,H_i>=lb_i,H_i<ub_i],axis=0)
#        #dont have a goood idea at the moment, so manually:
        for i in range(len(n_smoothed)):
            if ix[i]==0:
                n_smoothed[i]=n[i]
        print("SVMAF changed the refractive index at " + str(sum(ix)) + " frequencies")
        return n_smoothed      
      
    def calculateinits(self,H,l):
        #the frequency axis that will be used for the calculation
        oldfreqaxis=self.H.getfreqs()
        #this could be done once for all to gain speed!
        intpH=interp1d(self.H_firstPuls.getfreqs(),self.H_firstPuls.fdData[:,1:],axis=0)
        newH=intpH(oldfreqaxis)
  
        #crop the time domain data to the first pulse and apply than the calculation
        n=self.n_0-newH[:,2]/(2*py.pi*oldfreqaxis*l)*c
        alpha=-c/(2*py.pi*oldfreqaxis*l)*py.log(abs(newH[:,0])*(n+1)**2/(4*n))

        return n.real,alpha.real
    
    def plotInits(self,H,l,figurenumber=200):
        inits=self.calculateinits(H,l)
        py.figure(figurenumber)
        py.title('Initial Conditions')
        py.plot(self.H.getfreqsGHz(),inits[0])
        py.plot(self.H.getfreqsGHz(),-inits[1])
        py.xlabel('Frequency in GHz')
        py.ylabel('optical constant value')
        py.legend(('n_real','n_imag'))
    
    def getLengthes(self):
        return py.linspace(self.thicknessestimate-self.thicknessstd,self.thicknessestimate+self.thicknessstd,self.steps)
   
    def calculaten(self,H,l):
        inits=py.asarray(self.calculateinits(H,l))
        res=[]
        vals=[]
        bnds=((1,None),(0,None))
        nums=len(H[:,0])
#        t=minimize(self.error_func,inits[:,0],args=(H[0,:],l), method='SLSQP',\
#        bounds=bnds, options={'ftol':1e-9,'maxiter':2000, 'disp': False})
#        res.append(t.x[0]-1j*t.x[1])
#        vals.append(t.fun)

        for i in range(nums):
            
#            t=minimize(self.error_func,[inits[0,i],inits[1,i]],args=(H[i,:2],l), method='SLSQP',\
#            bounds=bnds, options={'ftol':1e-9,'maxiter':2000, 'disp': False})
            t=minimize(self.error_func,[inits[0,i],inits[1,i]],args=(H[i,:2],l), method='Nelder-Mead',
            options={'xtol': 1e-6,'disp':False})
            res.append(t.x[0]-1j*t.x[1])
            vals.append(t.fun)
        n=py.asarray(res)
        self.n=py.column_stack((H[:,0],n,n))
        return n

        
    def doCalculation(self,bool_findl=1,n_SVMAFS=5,bool_silent=0):
        bw=self.H.getBandwidth()
        self.H.manipulateFDData(-1,[bw[0]+50e9,bw[1]-200e9])
        if bool_findl:
            self.l_opt=self.findLintelli()

        print('\033[92m\033[1m' + '  Use Sample Thickness: ' + str(self.l_opt*1e6) + ' micro m ' + '\033[0m')


        n=self.calculaten(self.H.fdData,self.l_opt)
        n_smoothed=n
        i=0
        
        while i<n_SVMAFS:
            n_smoothed=self.SVMAF(self.H.getfreqs(),n_smoothed,self.l_opt)
            i+=1

        self.n=py.column_stack((self.H.getfreqs(),n,n_smoothed))        
        
        return self.n
        
    def getFilenameSuggestion(self):
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
        
        
    def saveResults(self,filename=None):
        
        H_theory=self.H_theory(self.H.getfreqs(),[self.n[:,1].real,self.n[:,1].imag],self.l_opt)        
        #built the variable that should be saved:        
        savetofile=py.column_stack((
        self.n[:,0].real,
        self.n[:,1].real,self.n[:,1].imag, #the real and imaginary part of n
        self.n[:,2].real,-self.n[:,2].imag, #the real and imaginary part of the smoothed n
        self.H.fdData[:,1].real,self.H.fdData[:,1].imag,#H_measured
        self.H.fdData[:,4].real,self.H.fdData[:,5].real,#uncertainties
        self.H.fdData[:,2].real,self.H.fdData[:,3].real,#absH,ph H measured    
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
        
            
            
    def plotRefractiveIndex(self,bool_plotsmoothed=1,savefig=0,filename=None):
        if filename==None:
            figname_b=self.getFilenameSuggestion()
        else:
            figname_b=filename
            
        py.figure(1)
        if bool_plotsmoothed==1:
            py.plot(self.n[:,0].real/1e9,self.n[:,2].real)
        else:
            py.plot(self.n[:,0].real/1e9,self.n[:,1].real)
        
        py.xlabel('Frequency in GHz')
        py.ylabel('n Real')
        if savefig:
            figname=figname_b+'n_real.png'
            py.savefig(figname,dpi=200)
            
        py.figure(2)
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
        py.figure()
        resolution=300
        ix=py.argmin(abs(freq-self.H.getfreqs()))
        n_i=py.linspace(0.0,0.05,int(resolution/2))
        n_r=py.linspace(1,3,resolution)
        N_R,N_I=py.meshgrid(n_r,n_i)
        E_fu=py.zeros((len(n_i),len(n_r)))
        for i in range(len(n_r)):
            for k in range(len(n_i)):
                E_fu[k,i]=self.error_func([n_r[i],n_i[k]],self.H.fdData[ix,:],l)
#        print(E_fu[:,2])
        py.pcolor(N_R,N_I,py.log10(E_fu))
        py.colorbar()


def getparams(name):
    #this function should return for several datasets to try the programm the following variables:
    #reffiles: list of filenames of the reference files
    #samfiles: list of filenames of the sample files
    #mode: Different file structures, i.e. INRIM, or Marburg or lucastestformat
    #thickness: the measured thickness of the sample

    mode='Marburg'
    teralyzer=py.zeros((3,3))
    path1='/home/jahndav/Dropbox/THz-Analysis/NEW INRIM Measurements/'
    path2='/home/jahndav/Dropbox/THz-Analysis/'
    if name=='siliInrim':
        
    #    the silicon sample measured at inrim
        thickness=330e-6 
        samfiles=glob.glob(path2+'silicon/*sam*')
        reffiles=glob.glob(path2+'silicon/*ref*') 
        mode='INRIM'
        
    elif name=='siliInrim1':
        thickness=330e-6
        samfiles=glob.glob(path2+'/INRIM1/*sam*.dat')
        reffiles=glob.glob(path2+'/INRIM1/*ref*.dat') 
        mode='INRIM'
    elif name=='Lactose3Inrim':
       
        thickness=2376e-6
        samfiles=glob.glob(path+'Lactose/*lactose-3_HalfStep_RH11_sam*.dat')
        reffiles=glob.glob(path+'Lactose/*lactose-3_HalfStep_RH11_ref*.dat')
        mode='INRIM'
    
    elif name=='Lactose3Inrim2':
        thickness=2376e-6
        samfiles=glob.glob(path+'Lactose/*lactose-3_fast_half-step_RH14_sam*_fast*.dat')
        reffiles=glob.glob(path+'Lactose/*lactose-3_fast_half-step_RH14_ref*_fast*.dat')
        mode='INRIM'
        
    elif name=='Teflon2Inrim':
        thickness=2117e-6
        samfiles=glob.glob(path+'*Teflon-II-4_half-step_RH16_sam*_fast.dat')
        reffiles=glob.glob(path+'*Teflon-II-4_half-step_RH16_ref*_fast.dat')
        mode='INRIM'
        
    elif name=='Teflon3':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=3905e-6
        samfiles=glob.glob(path2+'MarburgData/*_TeflonI-3*')
        reffiles=glob.glob(path2+'MarburgData/*ref7*')+glob.glob(path2+'MarburgData/*ref8*')
    elif name=='Teflon2':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=2117e-6
        samfiles=glob.glob(path2+'MarburgData/*_TeflonI-2*')
        reffiles=glob.glob(path2+'MarburgData/*ref6*')+glob.glob(path2+'MarburgData/*ref7*')
    elif name=='Teflon1':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=2470e-6
        samfiles=glob.glob(path2+'MarburgData/*_TeflonI-1*')
        reffiles=glob.glob(path2+'MarburgData/*ref5*')+glob.glob(path2+'MarburgData/*ref6*')
    elif name=='PP1':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=2039e-6
        samfiles=glob.glob(path2+'MarburgData/*_PPI-1*')
        reffiles=glob.glob(path2+'MarburgData/*ref3*')+glob.glob(path2+'MarburgData/*ref4*')
    elif name=='PP2':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=1294e-6
        samfiles=glob.glob(path2+'MarburgData/*_PPI-2*')
        reffiles=glob.glob(path2+'MarburgData/*ref4*')+glob.glob(path2+'MarburgData/*ref5*')
    elif name=='Lactose1':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=2166e-6
        samfiles=glob.glob(path2+'MarburgData/*_Lact1*')
        reffiles=glob.glob(path2+'MarburgData/*ref0*')+glob.glob(path2+'MarburgData/*ref1*')
    elif name=='Lactose2':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=1611e-6
        samfiles=glob.glob(path2+'MarburgData/*_Lact2*')
        reffiles=glob.glob(path2+'MarburgData/*ref1*')+glob.glob(path2+'MarburgData/*ref2*')
    elif name=='Lactose3':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=2376e-6
        samfiles=glob.glob(path2+'MarburgData/*_Lact3*')
        reffiles=glob.glob(path2+'MarburgData/*ref2*')+glob.glob(path2+'MarburgData/*ref3*')        
        teralyzer=py.loadtxt(path2+'MarburgData/L3.txt.csv',delimiter=',',usecols=(0,1,2),skiprows=1)
    elif name=='Ceramic1':
 #   a Teflon sample, flat real(n) and imag(n) expected
        folder='/home/jahndav/Dropbox/laboratory (1)/data-analysis/python_THz/DataSets/ceramic/'
        thickness=236e-6
        samfiles=glob.glob(folder+'*_sample1*')
        reffiles=glob.glob(folder+'*_ref1*')+glob.glob(folder+'*_ref2_*')
    elif name=='Ceramic4':
 #   a Teflon sample, flat real(n) and imag(n) expected
        folder='/home/jahndav/Dropbox/laboratory (1)/data-analysis/python_THz/DataSets/ceramic/'
        thickness=292e-6
        samfiles=glob.glob(folder+'*_sample4*')
        reffiles=glob.glob(folder+'*_ref4*')+glob.glob(folder+'*_ref5_*') 
    else:    
 #   the testdata used to write the code  
        thickness=677e-6 
        samfiles=glob.glob(path2+'rehi/Sample*')
        reffiles=glob.glob(path2+'rehi/Reference*') 
        mode='lucastestformat'
        teralyzer=py.loadtxt(path2+'rehi/Rehi_Teralyzer_OK.txt')
    
    return thickness,samfiles,reffiles,mode,teralyzer


    
if __name__=="__main__":
    #the Data is at the moment not allowed to be zero padded, 
    #it somehow changes the phase of the signal, and though the initial conditions
    #are shit!    
    

    #Load Parameters from getparams
    thickness,samfiles,reffiles,mode,teralyzer=getparams('Lactose2')
    #depending on format use different import module
    if mode=='lucastestformat':
        reftd=THzTdData(reffiles)
        samtd=THzTdData(samfiles)
    
    if mode=='Marburg':
        reftd=ImportMarburgData(reffiles)
        samtd=ImportMarburgData(samfiles)
    if mode=='INRIM':
        reftd=ImportInrimData(reffiles)
        samtd=ImportInrimData(samfiles)

#    #initialize the fd_data objects
    ref_fd=FdData(reftd)        
    sam_fd=FdData(samtd)
#    #initialize the mdata object (H,and so on)
    mdata=HMeas(ref_fd,sam_fd)
    mdata.manipulateFDData(-1e9,[200e9,2.2e12])
#    peaks=mdata.findAbsorptionLines()
#    py.plot(mdata.getfreqsGHz(),20*py.log10(mdata.fdData[:,2]))
#    py.plot(mdata.getfreqsGHz()[peaks],20*py.log10(mdata.fdData[peaks,2]),'o')
    myana=teralyz(mdata,thickness-30e-6,0.5*thickness,30)
#    myana.doCalculation()
#    myana.plotRefractiveIndex(1,1)
#    inits=myana.calculateinits(myana.H,myana.l_opt)
#    py.figure(1)
##    py.plot(myana.n[:,0].real/1e9,myana.n[:,1])
#    py.plot(myana.n[:,0].real/1e9,inits[0])
#    py.legend(('full H','Initial Conditions'))
#    
#    py.figure(2)
##    py.plot(myana.n[:,0].real/1e9,myana.n[:,1])
#    py.plot(myana.n[:,0].real/1e9,inits[1])
#    py.legend(('full H','Initial Conditions'))
    