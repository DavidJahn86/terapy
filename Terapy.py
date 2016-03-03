import os
from uncertainties import unumpy
from scipy.interpolate import interp1d
from scipy.constants import c
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
import TeraData as TD
import glob

n_0=1.00027-0.0000j

def getThicknessEstimateDavid(H,fmin=200e9,fmax=1e12):
    '''
    This method returns an estimate for the thickness of the sample under investigation
    fmin: minimal frequency used for determining the Etalon Frequency
    fmax: maximal frequency used for determining the Etalon Frequency
    Idea:   a) The periodicity of abs(H) is related to the optical thickness
            b) The slope of the phase of H is related to the distance of the pulses
    combining both yields a (good?) estimate for the thickness of the sample
    
    '''
    rdata=H.getCroppedData(fmin,fmax)
    #calculate phase change
    p=np.polyfit(rdata.getFrequenciesRef(),rdata.getPhasesRef(),1)        
    kappa=abs(p[0])
    df=H.getEtalonSpacing()
    #assume air around
    n=1.0/(1.00027-kappa*df/np.pi)
    l=c/2.0*(1.00027/df-kappa/np.pi)
    return [l,n]

def getRefractiveIndexEstimateTimeDomain(tdsam,tdref,thickness):
    '''This method returns an estimate for the "Time Domain Refractive Index" of the sample
    it uses the pulse maxima of reference and sample measurement and a thickness for returning
    an averaged refractive index
    '''

    #the time delay between sample and reference pulse is used to estimate n
    tmaxs=tdsam.getPeakPosition()
    tmaxr=tdref.getPeakPosition()
    n=1+c/thickness*abs(tmaxs-tmaxr)
  
    return n

def getNumberofEchos(tdData,n,l):
    '''
    Use l and average n for calculating the maximal number of fabry-pero pulses
    present in the time window'''
    #time window
    t_max=tdData.getTimeWindowLength()
    delta=0.5*(t_max*c/(n*l))
    return int(delta)
  
def calculateInits(H,l):
    '''
    This uses the approximation that imaginary part of n is small, and solves
    the fresnell equation in this approximation analytically. Fabry-perot reflections are cut
    take care that the phase has ben offset_removed!'''
    
    
    omega=2*np.pi*H.getFrequenciesRef()
    
    #crop the time domain data to the first pulse and apply than the calculation
    n=n_0.real-H.getPhasesRef()/(omega*l)*c
    alpha=-c/(omega*l)*np.log(abs(H.getSpectrumRef())*(n+n_0.real)**2/(4*n*n_0.real))

    return n,alpha

def plotInits(H,l):
    #plots the initial conditions
    inits=calculateInits(H,l)
    plt.title('Initial Conditions')
    plt.plot(H.getFrequenciesRef()/1e9,inits[0])
    plt.plot(H.getFrequenciesRef()/1e9,-inits[1])
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('optical constant value')
    plt.legend((r'$n_{real}$',r'$n_{imag}$'))
    
    
def calculaten(H,l,echos):
    #this calculates n for a given transferfunction H and a given thickness l

    #calculate the initial values
    inits=np.asarray(calculateInits(H,l))
    
    res=[] #the array of refractive index
    vals=[] # the array of error function values
    bnds=((1,None),(0,None)) #not used at the moment, with SLSQP bnds can be introduced
    freqs=H.getFrequenciesRef()
    spec=H.getSpectrumRef()
    nums=len(freqs)
    #minimize the deviation of H_theory to H_measured for each frequency
    for i in range(nums):
#            t=minimize(self.error_func,[inits[0,i],inits[1,i]],args=(H[i,:2],l), method='SLSQP',\
#            bounds=bnds, options={'ftol':1e-9,'maxiter':2000, 'disp': False})
        t=minimize(error_func,[inits[0,i],inits[1,i]],args=(freqs[i],spec[i],l,echos), method='Nelder-Mead',
        options={'xtol': 1e-6,'disp':False})
        res.append(t.x[0]-1j*t.x[1])
        vals.append(t.fun)
    n=np.asarray(res)
    #self.n is a 5xlengthf array, frequency,n_real,n_imag,n_smoothed_real,n_smoothed_imag
    return n

def error_func(n,f,spectrum,l,echos):
    #the quadratic deviation of Htheroy and Hmeasuered, used for finding n
    H_t=H_theory(f,n,l,echos)
    return abs(H_t-spectrum)**2
    
def H_theory(freq,n,l,echos):
    #calculate the theoretic transfer function, so far only one medium!
    
    nc=n[0]-1j*n[1]
    r_as=(n_0-nc)/(n_0+nc)
   
    FPE=0
    P=lambda tn: np.exp(-1j*l*freq*2*np.pi*tn/c)
    for sumk in range(1,echos):
        FPE+=((r_as**2)*(P(nc)**2))**sumk
    
    H=4*n_0*nc/(nc+n_0)**2*P((nc-n_0))*(1+FPE)
    return H

def doCalculation(tdref,tdsam,measuredthickness=0,fmin=200e9,fmax=2e12,frequencyResolution=0
                ,bool_findl=1,n_SVMAFS=5,bool_silent=0):
    
    '''This function calculates the optical constants and the thickness of a sample
    pass the preprocessed reference and sample data and a thickness for starting the algorithm
    Select the range for calculating the 
    '''

    
    #H=H.getCroppedData(150e9,2e12)

    #H.plotme()
    #plotInits(H,l)
    
    #calculaten(H,l,noEchos)
    
    
    #if l_opt shouldn't be calculated used bool_findl=False
    #if bool_findl:
     #   l_opt=findLintelli()

    #print('\033[92m\033[1m' + '  Use Sample Thickness: ' + str(l_opt*1e6) + ' micro m ' + '\033[0m')

    #noEchos=getNumberofEchos()
    #calculate n for the given l_opt
    #n=calculaten(H,l_opt)
    #n_smoothed=n
    #i=0
    
    #smooth it with SVMAF
    #while i<n_SVMAFS:
    #    n_smoothed=self.SVMAF(self.H.getfreqs(),n_smoothed,self.l_opt)
    #    i+=1

    #self.n=np.column_stack((self.H.getfreqs(),n,n_smoothed))   
    
    #self.calculateinitsunc(self.H.fdData,self.l_opt)
    
    return n

def findLintelli(H,fmax,echos,n_init,l_init):
    
    #oscillation period can be calculated from H data!
    f_span=c/2*1/l_init*1/n_init #(this should be the length of the oscillation)
    #restrict H to only 5 oscillations
    H_small=H.getCroppedData(fmax-f_span*1,fmax+f_span*4)
    plt.figure(33)
    #minimize quasispace/totalvariation value
    t=minimize(errorL,l_init,args=((H_small,echos,)),\
    method='Nelder-Mead', options={'xtol':1e-6,'disp': False})#, 'disp': False})
    return t.x[0]

def errorL(l,H,echos):
    #the error function for the length finding minimization
    
    #calculate n for the short transfer function H and the length l
    n_small=[calculaten(H,l,echos)]
    #evaluate the quasi space value
    qs=quasiSpace(n_small,H.getfbins(),l)
    #evaluate the total variation value
    tv=totalVariation(n_small)
    #plot them
    plt.plot(l,qs[0],'+')
    plt.plot(l,tv[0],'*')        
    print("Currently evaluating length: "+ str(l[0]*1e6) + " TV Value " + str(tv[0]))
    return tv[0]
    
def quasiSpace(ns,df,ls):
    #evaluates the quasi space value
    allqs=[]
    #at the moment, the most easy method(everything else failed...)
    #get the mean absolute value of the complete spectrum. works best
    for i in range(len(ns)):
#            xvalues=py.fftfreq(len(ns[i]),df)*c/4

        QSr=np.fft.fft(ns[i].real-np.mean(ns[i].real))
        QSi=np.fft.fft(ns[i].imag-np.mean(ns[i].real))
        #naive cut:
        ix=list(range(3,int(len(QSr)/2-1)))
    
        QSr=QSr[ix]
        QSi=QSi[ix]
        allqs.append(np.mean(abs(QSr))+np.mean(abs(QSi)))
    
    return allqs
    
def totalVariation(ns):
    #calculate the total variation value
    allvs=[]
    
    for i in range(len(ns)):
        tv1=0            
        for dm in range(1,len(ns[i])):
            tv1+=abs(ns[i][dm-1].real-ns[i][dm].real)+\
            abs(ns[i][dm-1].imag-ns[i][dm].imag)
        
        allvs.append(tv1)
        
    return allvs

def getHFirstPuls(tdRef,tdSam):
    #returns the transferfunction of the timedomain data that corresponds to the first
    #pulse only
    origlen=tdRef.getSamplingPoints()
    tdRef=tdRef.getFirstPuls(10e-12)
    tdSam=tdSam.getFirstPuls(5e-12)
    #bring it to the old length
    tdRef=tdRef.zeroPaddData(origlen-tdRef.getSamplingPoints())
    tdSam=tdSam.zeroPaddData(origlen-tdSam.getSamplingPoints())
    #calculate the fouriertransform
    firstref=TD.FrequencyDomainData.fromTimeDomainData(tdRef)
    firstsam=TD.FrequencyDomainData.fromTimeDomainData(tdSam)
    
    #calculate the transferfunction for the first puls
    H_firstPuls=TD.FrequencyDomainData.divideTwoSpectra(firstsam,firstref)
    
    return H_firstPuls

def getFilenameSuggestion(tdSam):
    #tries to give a valid filename from the sample filenames
    filenames=tdSam.getDataSetName()
    av=filenames.split(',')
    if len(av)>1:
        sug='Average'+os.path.splitext(av[1])[0]
    else:
        sug=os.path.splitext(av[0])[0]
   
    return sug



def SVMAF(self,freq,n,l):
    #Apply the SVMAF filter to the material parameters
    runningMean=lambda x,N: np.hstack((x[:N-1],np.convolve(x,np.ones((N,))/N,mode='same')[N-1:-N+1],x[(-N+1):]))
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
    ix=np.all([H_r>=lb_r,H_r<ub_r,H_i>=lb_i,H_i<ub_i],axis=0)
#        #dont have a goood idea at the moment, so manually:
    for i in range(len(n_smoothed)):
        if ix[i]==0:
            n_smoothed[i]=n[i]
    print("SVMAF changed the refractive index at " + str(sum(ix)) + " frequencies")
    return n_smoothed      



def saveResults(filename):
    #save the results to a file        
    H_theory=self.H_theory(self.H.getfreqs(),[self.n[:,1].real,self.n[:,1].imag],self.l_opt)        
    #built the variable that should be saved:        
    savetofile=np.column_stack((
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
    np.savetxt(fname,savetofile,delimiter=',',header=headerstr)
    
    # Save in a separate file the results of the simplified calculation of n, k, alpha end epsilon along with their unc
    headerstr=('freq, ' 
    'ref_ind_real, ref_ind_imag, '
    'u(ref_ind_real), u(ref_ind_imag), '
    'std(l)_n, std(Esam)_n, std(Eref)_n, f(Theta)_n, f(k<<<)_n, f(FP)_n, f(n0)_n, '
    'std(l)_k, std(Esam)_k, std(Eref)_k, f(Theta)_k, f(k<<<)_k, f(FP)_k, f(n0)_k, '
    'Epsilon_1, Epsilon_2, '
    'u(Epsilon_1), u(Epsilon_2), '
    'alpha, u(alpha), alpha_max, ')
    
    if filename==None:
        fname=self.getFilenameSuggestion()
    else:
        fname=filename+"_"
        
    fname+='SimplifiedAnalysis_' + 'D=' + str(self.l_opt/1e-6) +'.txt'
    np.savetxt(fname,self.n_with_unc,delimiter=',',header=headerstr)

def plotErrorFunction(self,l,freq):
    #plots the error function
    plt.figure()
    resolution=300
    ix=np.argmin(abs(freq-self.H.getfreqs()))
    n_i=np.linspace(0.0,0.05,int(resolution/2))
    n_r=np.linspace(1,3,resolution)
    N_R,N_I=plt.meshgrid(n_r,n_i)
    E_fu=np.zeros((len(n_i),len(n_r)))
    for i in range(len(n_r)):
        for k in range(len(n_i)):
            E_fu[k,i]=self.error_func([n_r[i],n_i[k]],self.H.fdData[ix,:3],l)
#        print(E_fu[:,2])
    plt.pcolor(N_R,N_I,np.log10(E_fu))
    plt.colorbar()


if __name__=='__main__':
    
    fns=glob.glob('*ref*.txt')
    tdRef=TD.importMarburgData(fns)
    fns=glob.glob('*_2_*.txt')
    tdSam=TD.importMarburgData(fns)
    
    tdRef=tdRef.getWindowedData(5e-12)
    tdSam=tdSam.getWindowedData(5e-12)
    
  #  doCalculation(ref,sam)
   
    fmin=300e9
    fmax=2e12
    
    bool_findl=True
    
    #resolution=5e9
   
    fRef=TD.FrequencyDomainData.fromTimeDomainData(tdRef)
    fSam=TD.FrequencyDomainData.fromTimeDomainData(tdSam)
    
    H=TD.FrequencyDomainData.divideTwoSpectra(fSam,fRef)
    l,n=getThicknessEstimateDavid(H,fmin,fmax)
    
    noEchos=getNumberofEchos(tdRef,n,l)

    H_firstpulse=getHFirstPuls(tdRef,tdSam)
    H_firstpulse=H_firstpulse.getCroppedData(fmin,fmax)
    H=H.getCroppedData(fmin,fmax)
    
    #inits=calculateInits(H_firstpulse,l)
    #plotInits(H_firstpulse,l)
    
    fmax=fRef.getMaxFreq()
    
    #if l_opt shouldn't be calculated used bool_findl=False
    if bool_findl:
        l_opt=findLintelli(H,fmax,echos,n,l):

    print('\033[92m\033[1m' + '  Use Sample Thickness: ' + str(l_opt*1e6) + ' micro m ' + '\033[0m')

    #calculate n for the given l_opt
    n=calculaten(H,l_opt)
    #n_smoothed=n
    #i=0
    
    #smooth it with SVMAF
    #while i<n_SVMAFS:
    #    n_smoothed=self.SVMAF(self.H.getfreqs(),n_smoothed,self.l_opt)
    #    i+=1

    #self.n=np.column_stack((self.H.getfreqs(),n,n_smoothed))   
    
    #self.calculateinitsunc(self.H.fdData,self.l_opt)    
    
    
   #tera=teralyz(ref,sam)
    #l,n=tera.getThicknessEstimateDavid()
    
    #tera.calculaten(l)

'''    
   def plotRefractiveIndex(self,bool_plotsmoothed=1,savefig=0,filename=None):
        #plot the refractive index
    
        if filename==None:
            figname_b=self.getFilenameSuggestion()
        else:
            figname_b=filename
            
        plt.figure('Refractive_Real_Plot')
        if bool_plotsmoothed==1:
            plt.plot(self.n[:,0].real/1e9,self.n[:,2].real)
        else:
            plt.plot(self.n[:,0].real/1e9,self.n[:,1].real)
        
        plt.xlabel('Frequency in GHz')
        plt.ylabel('n Real')
        if savefig:
            figname=figname_b+'n_real.png'
            plt.savefig(figname,dpi=200)
            
        plt.figure('Refractive_Imag_Plot')
        if bool_plotsmoothed==1:
            plt.plot(self.n[:,0].real/1e9,-self.n[:,2].imag)
        else:
            plt.plot(self.n[:,0].real/1e9,-self.n[:,1].imag)

        plt.xlabel('Frequency in GHz')
        plt.ylabel('n Imag')
        if savefig:
            figname=figname_b+'n_imag.png'
            plt.savefig(figname,dpi=200)
            
        plt.figure('Refractive_Simplified_Real_Plot')
        plt.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1], color='red')
        plt.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1]-self.n_with_unc[:,3],self.n_with_unc[:,1]+self.n_with_unc[:,3], alpha=0.5, label='k = 1', facecolor='blue')
        plt.xlabel('Frequency in GHz')
        plt.ylabel('n Real')
        if savefig:
            figname=figname_b+'n_simplified_real.png'
            plt.savefig(figname,dpi=200)
            
        plt.figure('Refractive_Simplified_Imag_Plot')
        plt.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,2], color='red')
        plt.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,2]-self.n_with_unc[:,4],self.n_with_unc[:,2]+self.n_with_unc[:,4], alpha=0.5, label='k = 1', facecolor='blue')
        plt.xlabel('Frequency in GHz')
        plt.ylabel('n Imag')
        if savefig:
            figname=figname_b+'n_simplified_imag.png'
            plt.savefig(figname,dpi=200)
        
        plt.figure('Eqs_Comparison_Real_Plot')
        if bool_plotsmoothed==1:
            plt.plot(self.n[:,0].real/1e9,self.n[:,2].real, color='green', label='full eq')
        else:
            plt.plot(self.n[:,0].real/1e9,self.n[:,1].real, color='green', label='full eq')
        plt.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1], color='red', label='simplified eq')
        plt.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1]-self.n_with_unc[:,3],self.n_with_unc[:,1]+self.n_with_unc[:,3], alpha=0.5, facecolor='blue')
        plt.legend(loc='upper left')        
        plt.xlabel('Frequency in GHz')
        plt.ylabel('n Real')
        
        plt.figure('Eqs_Comparison_Imag_Plot')
        if bool_plotsmoothed==1:
            plt.plot(self.n[:,0].real/1e9,-self.n[:,2].imag, color='green', label='full eq')
        else:
            plt.plot(self.n[:,0].real/1e9,-self.n[:,1].imag, color='green', label='full eq')
        plt.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,2], color='red', label='simplified eq')
        plt.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,2]-self.n_with_unc[:,4],self.n_with_unc[:,2]+self.n_with_unc[:,4], alpha=0.5, facecolor='blue')
        plt.legend(loc='upper left')        
        plt.xlabel('Frequency in GHz')
        plt.ylabel('n Imag')
        
        # Plot uncertainty components of n
        plt.figure('Refractive_Simplified_Real_Unc_Components_Plot')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,3], 'r', label='combined unc')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,5], 'g', label='s(l)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,6], 'b', label='s(Esam)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,7], 'y', label='s(Eref)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,8], '--r', label='f($\Theta$)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,9], '--g', label='f(k<<<)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,10], '--b', label='f(FP)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,11], '--y', label='f(n_0)')
        plt.legend(loc='upper left')
        plt.xlabel('Frequency in GHz')
        plt.ylabel('Uncertainty components')
        if savefig:
            figname=figname_b+'n_simplified_real_unc_components.png'
            plt.savefig(figname,dpi=200)
        
        # Plot uncertianty components of k
        plt.figure('Refractive_Simplified_Imag_Unc_Components_Plot')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,4], 'r', label='combined unc')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,12], 'g', label='s(l)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,13], 'b', label='s(Esam)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,14], 'y', label='s(Eref)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,15], '--r', label='f($\Theta$)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,16], '--g', label='f(k<<<)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,17], '--b', label='f(FP)')
        plt.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,18], '--y', label='f(n_0)')
        plt.legend(loc='upper left')
        plt.xlabel('Frequency in GHz')
        plt.ylabel('Uncertainty components')
        if savefig:
            figname=figname_b+'n_simplified_imag_unc_components.png'
            plt.savefig(figname,dpi=200)
            
        #plt.figure('Comparison_with_NPL_Real_Plot')
        #if bool_plotsmoothed==1:
        #    plt.plot(self.n[:,0].real/1e9,self.n[:,2].real, color='green', label='full eq')
        #else:
        #    plt.plot(self.n[:,0].real/1e9,self.n[:,1].real, color='green', label='full eq')
        #plt.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1], color='red', label='simplified eq')
        #plt.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1]-self.n_with_unc[:,3],self.n_with_unc[:,1]+self.n_with_unc[:,3], alpha=0.5, facecolor='blue')
        #plt.plot(self.n_with_unc[:,0]/1e9,[1.95]*len(self.n_with_unc[:,0]), color='black', label='NPL')
        #plt.fill_between(self.n_with_unc[:,0]/1e9,plt.asarray([1.95]*len(self.n_with_unc[:,0]))-0.05,plt.asarray([1.95]*len(self.n_with_unc[:,0]))+0.05, alpha=0.5, facecolor='yellow')
        #plt.legend(loc='upper left')        
        #plt.xlabel('Frequency in GHz')
        #plt.ylabel('n Real')
        
        plt.figure('Epsilon_Real_Plot')
        plt.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,19], color='red')
        plt.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,19]-self.n_with_unc[:,21],self.n_with_unc[:,19]+self.n_with_unc[:,21], alpha=0.5, facecolor='blue')        
        plt.xlabel('Frequency in GHz')
        plt.ylabel('Epsilon Real')
        if savefig:
            figname=figname_b+'Epsilon_real.png'
            plt.savefig(figname,dpi=200)
        
        plt.figure('Epsilon_Imag_Plot')
        plt.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,20], color='red')
        plt.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,20]-self.n_with_unc[:,22],self.n_with_unc[:,20]+self.n_with_unc[:,22], alpha=0.5, facecolor='blue')        
        plt.xlabel('Frequency in GHz')
        plt.ylabel('Epsilon Imag')
        if savefig:
            figname=figname_b+'Epsilon_imag.png'
            plt.savefig(figname,dpi=200)
        
        plt.figure('Alpha_Plot')
        plt.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,23], color='red', label='alpha / cm^-1')
        plt.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,23]-self.n_with_unc[:,24],self.n_with_unc[:,23]+self.n_with_unc[:,24], alpha=0.5, facecolor='blue')
        plt.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,25], '--r', label='alpha max / cm^-1')
        plt.xlabel('Frequency in GHz')
        plt.ylabel('Alpha')
        if savefig:
            figname=figname_b+'Alpha.png'
            plt.savefig(figname,dpi=200)
    
    def calculateinitsunc(self,H,l,sigma_L = 1e-6,sigma_Theta = 1,n_exact = 1,filename=None):
        # Calculates the uncertianty on n and k according to:
        # W. Withayachumnankul, B. M. Fisher, H. Lin, D. Abbott, "Uncertainty in terahertz time-domain spectroscopy measurement", J. Opt. Soc. Am. B., Vol. 25, No. 6, June 2008, pp. 1059-1072
        #
        # sigma_L = standard uncertainty on sample thickness in meter
        # sigma_Theta = interval of sample misallignment in degree
        # n_exact = exact value of the refractive index of air during the measurements
        
        n, k = self.calculateinits(H,l)
        n = py.asarray(n)
        k = py.asarray(k)      
        
        Asam = []
        Aref = []
        Bsam = []
        Bref = []
        for i in range(len(self.H.getfreqs())):
            Asam.append((py.sum(py.imag(self.H.fdsam.getFAbs().tolist()[i]*py.exp(1j*2*py.pi*self.H.getfreqs().tolist()[i]*self.H.fdsam._tdData.getTimes()))*self.H.fdsam._tdData.getUncEX())**2))
            Aref.append((py.sum(py.imag(self.H.fdref.getFAbs().tolist()[i]*py.exp(1j*2*py.pi*self.H.getfreqs().tolist()[i]*self.H.fdref._tdData.getTimes()))*self.H.fdref._tdData.getUncEX())**2))
            Bsam.append((py.sum(py.real(self.H.fdsam.getFAbs().tolist()[i]*py.exp(1j*2*py.pi*self.H.getfreqs().tolist()[i]*self.H.fdsam._tdData.getTimes()))*self.H.fdsam._tdData.getUncEX())**2))
            Bref.append((py.sum(py.real(self.H.fdref.getFAbs().tolist()[i]*py.exp(1j*2*py.pi*self.H.getfreqs().tolist()[i]*self.H.fdref._tdData.getTimes()))*self.H.fdref._tdData.getUncEX())**2))
        
        # Uncertainty on n
        sn_Esam_2 = ((c/(2*py.pi*self.H.getfreqs()*l))**2 * py.asarray(Asam)/self.H.fdsam.getFAbs()**4)/self.H.fdsam._tdData.numberOfDataSets
        sn_Eref_2 = ((c/(2*py.pi*self.H.getfreqs()*l))**2 * py.asarray(Aref)/self.H.fdref.getFAbs()**4)/self.H.fdref._tdData.numberOfDataSets
        sn_l_2 = ((n-self.n_0)*sigma_L/l)**2
        #sn_l_2_1 = (c*self.H.getFPh()/(2*py.pi*self.H.getfreqs()*l*l))**2 * sigma_L**2
        #sn_H_2 = (c/(2*py.pi*self.H.getfreqs()*l))**2 * self.H.getFPhUnc()**2
        fn_Theta = (n-self.n_0.real)*(1/py.cos(sigma_Theta*py.pi/180)-1)
        fn_H = (c/(2*py.pi*self.H.getfreqs()*l))*py.absolute(-py.angle(4*(n-1j*k)*self.n_0/(n-1j*k+self.n_0)**2))
        fn_FP = (c/(2*py.pi*self.H.getfreqs()*l))*py.absolute(-py.angle(1/(1-((n-1j*k-self.n_0)/(n-1j*k+self.n_0))**2*py.exp(-2*1j*(n-1j*k)*2*py.pi*self.H.getfreqs()*l/c))))
        fn_n0 = abs(self.n_0.real - n_exact)*py.ones(len(self.H.getFPh()))
        u_n = py.sqrt(sn_l_2+sn_Esam_2+sn_Eref_2)+fn_Theta+fn_H+fn_FP+fn_n0

        # Uncertianty on k
        sk_Esam_2 = ((c/(2*py.pi*self.H.getfreqs()*l))**2 *(Bsam/self.H.fdsam.getFAbs()**4 + ((n-self.n_0)/(n+self.n_0))**2 * sn_Esam_2/n**2))/self.H.fdsam._tdData.numberOfDataSets
        sk_Eref_2 = ((c/(2*py.pi*self.H.getfreqs()*l))**2 *(Bref/self.H.fdref.getFAbs()**4 + ((n-self.n_0)/(n+self.n_0))**2 * sn_Eref_2/n**2))/self.H.fdref._tdData.numberOfDataSets
        sk_l_2 = (k*sigma_L/l)**2 + (c*(n-self.n_0)/((n+self.n_0)*n*2*py.pi*self.H.getfreqs()*l))**2*sn_l_2
        #sk_l_2_1 = ((c/(2*py.pi*self.H.getfreqs()*l*l))*py.log(self.H.getFAbs()*(n+self.n_0.real)**2/(4*n*self.n_0.real)))**2 * sigma_L**2
        #sk_H_2 = (-c/(2*py.pi*self.H.getfreqs()*l*self.H.getFAbs()))**2 * self.H.getFAbsUnc()**2
        fk_Theta = k*(1/py.cos(sigma_Theta*py.pi/180)-1)+c*(n-self.n_0.real)*fn_Theta/(n*2*py.pi*self.H.getfreqs()*l*(n+self.n_0.real))
        fk_H = (c/(2*py.pi*self.H.getfreqs()*l))*(py.log(py.absolute(n/(n-1j*k)*((n-1j*k+self.n_0.real)/(n+self.n_0.real))**2))+py.absolute(fn_H)*(n-self.n_0.real)/(n*(n+self.n_0.real)))
        fk_FP = (c/(2*py.pi*self.H.getfreqs()*l))*(py.absolute(-py.log(py.absolute(1/(1-((n-1j*k-self.n_0)/(n-1j*k+self.n_0))**2*py.exp(-2*1j*(n-1j*k)*2*py.pi*self.H.getfreqs()*l/c)))))+py.absolute(fn_FP)*(n-self.n_0.real)/(n*(n+self.n_0.real)))
        fk_n0 = (c/(2*py.pi*self.H.getfreqs()*l))*(n-self.n_0.real)*(self.n_0.real - n_exact)/(n*self.n_0.real)
        u_k = py.sqrt(sk_l_2+sk_Esam_2+sk_Eref_2)+fk_Theta+fk_H+fk_FP+fk_n0
        
        # Convert n in epsilon Epsilon = Epsilon_1 + j Epsilon_2 = (n+jk)**2
        # Epsilon_1 = n**2 - k**2
        # Epsilon_2 = -2nk
        Epsilon_1 = n**2 - k **2
        Epsilon_2 = -2 * n * k
        u_Epsilon_1 = py.sqrt((2*n*u_n)**2 + (-2*k*u_k)**2)
        u_Epsilon_2 = py.sqrt((-2*k*u_n)**2 + (-2*n*u_k)**2)
        
        # Calculate absorption coefficient
        # alpha = 4 * pi * k * f / c1
        alpha = 4 * py.pi * k * self.H.getfreqs() / (100 * c)      # in cm^-1
        u_alpha = 4 * py.pi * u_k * self.H.getfreqs() / (100 * c)  # in cm^-1
        
        # Calculate maximum measurable absorption coefficient according to
        # P. U. Jepsen and B. M. Fisher: "Dynamic Range in terahertz time-domain transmission and reflection spectroscopy", Optics Letters, Vol. 30, n. 1, pp. 29-31, Jan 2005
        
        alpha_max = 2 * py.log((self.H.fdsam.getDR() * 4 * n)/(n + 1)**2) / (100 * l) # in cm^-1
        
        # Save results into a table accessible from outside
        self.n_with_unc=py.real(py.column_stack((
        self.H.getfreqs(),                            # frequencies
        n, k,                                         # real and imaginary part of n
        u_n, u_k,                                     # k=1 combined uncertainty on n and k
        py.sqrt(sn_l_2), py.sqrt(sn_Esam_2), py.sqrt(sn_Eref_2), fn_Theta, fn_H, fn_FP, fn_n0, # Uncertainty components of n due to thickness, H, sample misallignment, k<<<, Neglect FP, ref ind of air
        py.sqrt(sk_l_2), py.sqrt(sk_Esam_2), py.sqrt(sk_Eref_2), fk_Theta, fk_H, fk_FP, fk_n0, # Uncertainty components of k due to thickness, H, sample misallignment, k<<<, Neglect FP, ref ind of air
        Epsilon_1, Epsilon_2,                         # Real and imaginary part of Epsilon
        u_Epsilon_1, u_Epsilon_2,                     # k = 1 uncertainty on the real and imaginary part of Epsilon
        alpha, u_alpha,                               # Absorption coefficient and its k = 1 uncertainty
        alpha_max,                                     # Maximum measurable absorption coefficient
        )))
        return

#    def calculateFDunc(self):  
#        #Calculates the uncertainty of the FFT according to:
#        #   - J. M. Fornies-Marquina, J. Letosa, M. Garcia-Garcia, J. M. Artacho, "Error Propagation for the transformation of time domain into frequency domain", IEEE Trans. Magn, Vol. 33, No. 2, March 1997, pp. 1456-1459
#        
#        # Calculate the uncertainty of the real and imaginary part of H
#        #x=[a,b,c,d]
#        dfda=lambda x: x[2]/(x[2]**2+x[3]**2)
#        dfdb=lambda x: x[3]/(x[2]**2+x[3]**2)
#        dfdc=lambda x: (x[0]*x[3]**2-x[0]*x[2]**2-2*x[1]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
#        dfdd=lambda x: (x[1]*x[2]**2-x[1]*x[3]**2-2*x[0]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
#        
#        dgda=lambda x: -x[3]/(x[2]**2+x[3]**2)
#        dgdb=lambda x: x[2]/(x[2]**2+x[3]**2)
#        dgdc=lambda x: -(x[1]*x[2]**2-x[1]*x[3]**2+2*x[0]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
#        dgdd=lambda x: (x[0]*x[3]**2-x[0]*x[2]**2-2*x[1]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
#        
#        ta=self.fdsam.getFReal()
#        tb=self.fdsam.getFImag()
#        tc=self.fdref.getFReal()
#        td=self.fdref.getFImag()
#        
#        H_unc_real2=(self.fdsam.getFRealUnc()*dfda(py.asarray([ta,tb,tc,td])))**2
#        H_unc_real2+=(self.fdsam.getFImagUnc()*dfdb(py.asarray([ta,tb,tc,td])))**2
#        H_unc_real2+=(self.fdref.getFRealUnc()*dfdc(py.asarray([ta,tb,tc,td])))**2
#        H_unc_real2+=(self.fdref.getFImagUnc()*dfdd(py.asarray([ta,tb,tc,td])))**2
#        
#        H_unc_imag2=(self.fdsam.getFRealUnc()*dgda(py.asarray([ta,tb,tc,td])))**2
#        H_unc_imag2+=(self.fdsam.getFImagUnc()*dgdb(py.asarray([ta,tb,tc,td])))**2
#        H_unc_imag2+=(self.fdref.getFRealUnc()*dgdc(py.asarray([ta,tb,tc,td])))**2
#        H_unc_imag2+=(self.fdref.getFImagUnc()*dgdd(py.asarray([ta,tb,tc,td])))**2
#        
#        # Calculate the uncertainty of the modulus and phase of H
#        H_uncabs = py.sqrt((self.fdsam.getFAbsUnc()/self.fdref.getFAbs())**2 + (self.fdsam.getFAbs()*self.fdref.getFAbsUnc())**2/self.fdref.getFAbs()**4)
#        H_uncph = py.sqrt(self.fdsam.getFPhUnc()**2 + self.fdref.getFPhUnc()**2)
#        
#        return py.column_stack((py.sqrt(H_unc_real2),py.sqrt(H_unc_imag2),H_uncabs,H_uncph))
#        
#
'''
    