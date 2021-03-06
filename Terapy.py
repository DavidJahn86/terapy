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
        H_unc=self.calculateFDunc()

        #phase
        H_ph=self.fdref.getFPh()-self.fdsam.getFPh()
        #H itself
        H=(self.fdsam.getFReal()+1j*self.fdsam.getFImag())/(self.fdref.getFReal()+1j*self.fdref.getFImag())
        H=py.column_stack((self.fdref.getfreqs(),H.real,H.imag,abs(H),H_ph,H_unc))
        return H    

    def calculateFDunc(self):  
        #Calculates the uncertainty of the FFT according to:
        #   - J. M. Fornies-Marquina, J. Letosa, M. Garcia-Garcia, J. M. Artacho, "Error Propagation for the transformation of time domain into frequency domain", IEEE Trans. Magn, Vol. 33, No. 2, March 1997, pp. 1456-1459
        
        # Calculate the uncertainty of the real and imaginary part of H
        #x=[a,b,c,d]
        dfda=lambda x: x[2]/(x[2]**2+x[3]**2)
        dfdb=lambda x: x[3]/(x[2]**2+x[3]**2)
        dfdc=lambda x: (x[0]*x[3]**2-x[0]*x[2]**2-2*x[1]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
        dfdd=lambda x: (x[1]*x[2]**2-x[1]*x[3]**2-2*x[0]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
        
        dgda=lambda x: -x[3]/(x[2]**2+x[3]**2)
        dgdb=lambda x: x[2]/(x[2]**2+x[3]**2)
        dgdc=lambda x: -(x[1]*x[2]**2-x[1]*x[3]**2+2*x[0]*x[2]*x[3])/(x[2]**2+x[3]**2)**2
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
        
        # Calculate the uncertainty of the modulus and phase of H
        H_uncabs = py.sqrt((self.fdsam.getFAbsUnc()/self.fdref.getFAbs())**2 + (self.fdsam.getFAbs()*self.fdref.getFAbsUnc())**2/self.fdref.getFAbs()**4)
        H_uncph = py.sqrt(self.fdsam.getFPhUnc()**2 + self.fdref.getFPhUnc()**2)
        
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
        alpha=-c/(2*py.pi*oldfreqaxis*l)*py.log(newH[:,2]*(n+self.n_0.real)**2/(4*n*self.n_0.real))

        return n,alpha
        
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
        
        self.calculateinitsunc(self.H.fdData,self.l_opt)
        
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
            
        py.figure('Refractive_Simplified_Real_Plot')
        py.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1], color='red')
        py.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1]-self.n_with_unc[:,3],self.n_with_unc[:,1]+self.n_with_unc[:,3], alpha=0.5, label='k = 1', facecolor='blue')
        py.xlabel('Frequency in GHz')
        py.ylabel('n Real')
        if savefig:
            figname=figname_b+'n_simplified_real.png'
            py.savefig(figname,dpi=200)
            
        py.figure('Refractive_Simplified_Imag_Plot')
        py.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,2], color='red')
        py.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,2]-self.n_with_unc[:,4],self.n_with_unc[:,2]+self.n_with_unc[:,4], alpha=0.5, label='k = 1', facecolor='blue')
        py.xlabel('Frequency in GHz')
        py.ylabel('n Imag')
        if savefig:
            figname=figname_b+'n_simplified_imag.png'
            py.savefig(figname,dpi=200)
        
        py.figure('Eqs_Comparison_Real_Plot')
        if bool_plotsmoothed==1:
            py.plot(self.n[:,0].real/1e9,self.n[:,2].real, color='green', label='full eq')
        else:
            py.plot(self.n[:,0].real/1e9,self.n[:,1].real, color='green', label='full eq')
        py.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1], color='red', label='simplified eq')
        py.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1]-self.n_with_unc[:,3],self.n_with_unc[:,1]+self.n_with_unc[:,3], alpha=0.5, facecolor='blue')
        py.legend(loc='upper left')        
        py.xlabel('Frequency in GHz')
        py.ylabel('n Real')
        
        py.figure('Eqs_Comparison_Imag_Plot')
        if bool_plotsmoothed==1:
            py.plot(self.n[:,0].real/1e9,-self.n[:,2].imag, color='green', label='full eq')
        else:
            py.plot(self.n[:,0].real/1e9,-self.n[:,1].imag, color='green', label='full eq')
        py.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,2], color='red', label='simplified eq')
        py.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,2]-self.n_with_unc[:,4],self.n_with_unc[:,2]+self.n_with_unc[:,4], alpha=0.5, facecolor='blue')
        py.legend(loc='upper left')        
        py.xlabel('Frequency in GHz')
        py.ylabel('n Imag')
        
        # Plot uncertainty components of n
        py.figure('Refractive_Simplified_Real_Unc_Components_Plot')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,3], 'r', label='combined unc')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,5], 'g', label='s(l)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,6], 'b', label='s(Esam)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,7], 'y', label='s(Eref)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,8], '--r', label='f($\Theta$)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,9], '--g', label='f(k<<<)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,10], '--b', label='f(FP)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,11], '--y', label='f(n_0)')
        py.legend(loc='upper left')
        py.xlabel('Frequency in GHz')
        py.ylabel('Uncertainty components')
        if savefig:
            figname=figname_b+'n_simplified_real_unc_components.png'
            py.savefig(figname,dpi=200)
        
        # Plot uncertianty components of k
        py.figure('Refractive_Simplified_Imag_Unc_Components_Plot')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,4], 'r', label='combined unc')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,12], 'g', label='s(l)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,13], 'b', label='s(Esam)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,14], 'y', label='s(Eref)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,15], '--r', label='f($\Theta$)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,16], '--g', label='f(k<<<)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,17], '--b', label='f(FP)')
        py.semilogy(self.n_with_unc[:,0]/1e9, self.n_with_unc[:,18], '--y', label='f(n_0)')
        py.legend(loc='upper left')
        py.xlabel('Frequency in GHz')
        py.ylabel('Uncertainty components')
        if savefig:
            figname=figname_b+'n_simplified_imag_unc_components.png'
            py.savefig(figname,dpi=200)
            
        #py.figure('Comparison_with_NPL_Real_Plot')
        #if bool_plotsmoothed==1:
        #    py.plot(self.n[:,0].real/1e9,self.n[:,2].real, color='green', label='full eq')
        #else:
        #    py.plot(self.n[:,0].real/1e9,self.n[:,1].real, color='green', label='full eq')
        #py.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1], color='red', label='simplified eq')
        #py.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,1]-self.n_with_unc[:,3],self.n_with_unc[:,1]+self.n_with_unc[:,3], alpha=0.5, facecolor='blue')
        #py.plot(self.n_with_unc[:,0]/1e9,[1.95]*len(self.n_with_unc[:,0]), color='black', label='NPL')
        #py.fill_between(self.n_with_unc[:,0]/1e9,py.asarray([1.95]*len(self.n_with_unc[:,0]))-0.05,py.asarray([1.95]*len(self.n_with_unc[:,0]))+0.05, alpha=0.5, facecolor='yellow')
        #py.legend(loc='upper left')        
        #py.xlabel('Frequency in GHz')
        #py.ylabel('n Real')
        
        py.figure('Epsilon_Real_Plot')
        py.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,19], color='red')
        py.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,19]-self.n_with_unc[:,21],self.n_with_unc[:,19]+self.n_with_unc[:,21], alpha=0.5, facecolor='blue')        
        py.xlabel('Frequency in GHz')
        py.ylabel('Epsilon Real')
        if savefig:
            figname=figname_b+'Epsilon_real.png'
            py.savefig(figname,dpi=200)
        
        py.figure('Epsilon_Imag_Plot')
        py.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,20], color='red')
        py.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,20]-self.n_with_unc[:,22],self.n_with_unc[:,20]+self.n_with_unc[:,22], alpha=0.5, facecolor='blue')        
        py.xlabel('Frequency in GHz')
        py.ylabel('Epsilon Imag')
        if savefig:
            figname=figname_b+'Epsilon_imag.png'
            py.savefig(figname,dpi=200)
        
        py.figure('Alpha_Plot')
        py.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,23], color='red', label='alpha / cm^-1')
        py.fill_between(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,23]-self.n_with_unc[:,24],self.n_with_unc[:,23]+self.n_with_unc[:,24], alpha=0.5, facecolor='blue')
        py.plot(self.n_with_unc[:,0]/1e9,self.n_with_unc[:,25], '--r', label='alpha max / cm^-1')
        py.xlabel('Frequency in GHz')
        py.ylabel('Alpha')
        if savefig:
            figname=figname_b+'Alpha.png'
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
        py.savetxt(fname,self.n_with_unc,delimiter=',',header=headerstr)

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
