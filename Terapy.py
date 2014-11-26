import pylab as py
import copy
import glob
import os
from scipy.signal import argrelmax,argrelmin
from scipy.interpolate import interp1d
from scipy.constants import c
from scipy.optimize import minimize
from TeraData import *

class teradata():
 
    def __init__(self,FDref,FDsam):
        self.fdref=FDref
        self.fdsam=FDsam
        minrf,maxrf=self.fdref.getBandwidth()
        minsf,maxsf=self.fdsam.getBandwidth()
        
#        self.manipulateFDData(-1,[max(minrf,minsf,140e9),min(maxrf,maxsf)])
        self.H=self.calculateH()
        
    def manipulateFDData(self,fbins,fbnds):
        if fbins<self.fdref.getfbins() and fbins>0:
            self.fdref.setFDData(self.fdref.getInterpolatedFDData(fbins))
        if fbins<self.fdsam.getfbins() and fbins>0:
            self.fdsam.setFDData(self.fdsam.getInterpolatedFDData(fbins))
        
        self.fdref.setFDData(self.fdref.cropData(self.fdref.fdData,fbnds[0],fbnds[1]))        
        self.fdsam.setFDData(self.fdsam.cropData(self.fdsam.fdData,fbnds[0],fbnds[1]))
        self.H=self.calculateH()
   
    def calculateConfidenceInterval(self):  
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
                
        return py.sqrt(H_unc_real2),py.sqrt(H_unc_imag2)
    
    def checkDataIntegrity(self):
        if len(self.fdref.getfreqs())!=len(self.fdsam.getfreqs()) or\
            all(self.fdref.getfreqs()-self.fdsam.getfreqs()>0.1e9):
            return 0
        else:
            return 1         
      
    def interpolateData(self):
        #take care, this actually manipulates the real underlying data!
        #maybe consider to do a deepcopy!
    
        tdref=self.fdref.getassTDData()
        tdsam=self.fdsam.getassTDData()
        
        #naive approach: Interpolate on time axis basis, better: Interpolate in Frequency space!
        cmin=max(min(tdref.tdData[:,0]),min(tdsam.tdData[:,0]))
        cmax=min(max(tdref.tdData[:,0]),max(tdsam.tdData[:,0]))
        clen=min(tdref.getlength(),tdsam.getlength())
        
        #safe also old bnds        
        minrf,maxrf=self.fdref.getBandwidth()
        minsf,maxsf=self.fdsam.getBandwidth()
        
        tdrefnew=THzTdData(tdref.getInterData(tdref.tdData,clen,cmin,cmax),tdref.getfilename(),tdref._thzdata_raw,existing=True)
        tdsamnew=THzTdData(tdsam.getInterData(tdsam.tdData,clen,cmin,cmax),tdsam.getfilename(),tdsam._thzdata_raw,existing=True)
      
        
        self.fdref=FdData(tdrefnew,-1,[min(minrf,minsf),max(maxrf,maxsf)])
        self.fdsam=FdData(tdsamnew,-1,[min(minrf,minsf),max(maxrf,maxsf)])

    def calculateH(self):

        if not self.checkDataIntegrity():
            print "interpolation required"
            self.interpolateData()
        
        H_unc_real,H_unc_imag=self.calculateConfidenceInterval()
            #take care that abs(H) is always smaller one!
            #try if this increases the data quality!
#            H_ph=py.unwrap(py.angle(self.fdsam.fdData[:,1]/self.fdref.fdData[:,1]))
        H_ph=-self.fdsam.getUnwrappedPhase()+self.fdref.getUnwrappedPhase()
        H=py.column_stack((self.fdref.fdData[:,0],self.fdsam.fdData[:,1]/self.fdref.fdData[:,1],H_unc_real,H_unc_imag,H_ph))
        return H
      
    def getcroppedData(self,data,startfreq=400e9,endfreq=5e12):
        ix=py.all([data[:,0]>=startfreq,data[:,0]<endfreq],axis=0)
        return data[ix,:]
        
    def doPlots(self):
        freqs=self.getfreqsGHz()
        
        py.figure('H-UNC-Plot')
        py.plot(freqs,self.H[:,1].real)
        py.plot(freqs,self.H[:,1].imag)
        py.plot(freqs,self.H[:,1].real+self.H[:,2].real,'g--',freqs,self.H[:,1].real-self.H[:,2].real,'g--')
        py.plot(freqs,self.H[:,1].imag+self.H[:,3].real,'g--',freqs,self.H[:,1].imag-self.H[:,3].real,'g--')
        py.xlabel('Frequency in GHz')
        py.ylabel('Transfer Function')
        py.legend(('H_real','H_imag'))
#        
#        py.figure(figurenumber+1)
#        py.plot(freqs,abs(self.H[:,1]))
#        py.figure(figurenumber+2)
#        py.plot(freqs,self.H[:,-1])

    def getfreqsGHz(self):
        return self.H[:,0].real*1e9
        
    def findAbsorptionLines(self):
        
        #define treshhold
        tresh=6 #db treshhold to detect a line
        #bring to logarithmic scale
        hlog=20*py.log10(abs(self.H[:,1]))
        #remove etalon minima
        oldfreqs=self.H[:,0].real
        intpdata=interp1d(oldfreqs,hlog,'cubic')
        #take care that the moving average window length is always ap
        fnew=py.arange(min(oldfreqs),max(oldfreqs),1e9)
        hlog=intpdata(fnew)
#        window_size=30
#        window= py.ones(int(window_size))/float(window_size)
#        hlog=py.convolve(hlog, window, 'same')
        ixlines_prob=argrelmin(hlog)[0]
        ixlines=[]
        ixmax=argrelmax(hlog)[0]        
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
            if (leftdist[i]>tresh or rightdist[i]>tresh) and \
            hlog[ixlines_prob[i]]<-10:
                ixlines.append(ixlines_prob[i])        
        py.plot(fnew,hlog)
        py.plot(fnew[ixlines],hlog[ixlines],'*')        
        return fnew[ixlines],ixlines
    
    def getEtalonSpacing(self):
        #how to find a stable range! ? 
        bw=self.fdsam.getBandwidth()
        rdata=self.getcroppedData(self.H,max(bw[0],250e9),min(bw[1],3e12))  
        #get etalon frequencies:
        #need to interpolate data!
        oldfreqs=rdata[:,0].real       
        intpdata=interp1d(oldfreqs,abs(rdata[:,1]),'cubic')
        
        fnew=py.arange(min(oldfreqs),max(oldfreqs),0.1e9)
        absnew=intpdata(fnew)
        ixmaxima=argrelmax(absnew)[0]
        ixminima=argrelmin(absnew)[0]
        
        fmaxima=py.mean(py.diff(fnew[ixmaxima]))
        fminima=py.mean(py.diff(fnew[ixminima]))
        
        df=(fmaxima+fminima)*0.5 #the etalon frequencies
#        py.figure(17)
#        py.plot(fnew,absnew)
#        py.plot(fnew[ixmaxima],absnew[ixmaxima],'+')
#        py.plot(fnew[ixminima],absnew[ixminima],'+')
        print str(df/1e9) + " GHz estimated etalon frequency"
        return df
        
    def estimateLDavid(self):
        rdata=self.getcroppedData(self.H,200e9,2e12)
        #calculate phase change        
        kappa=abs(py.mean(py.diff(rdata[:,-1]))/(self.H[1,0]-self.H[0,0]))
        
        df=self.getEtalonSpacing()
        #assume air around
        n=1.0/(1.00027-kappa*df/py.pi)
        l=c/2.0*(1.00027/df-kappa/py.pi)
        return [l,n]
        
        
class teralyz():
   # n_0=1.00027#+#0.001j
    n_0=1.00027-0.0000j
   
    def __init__(self,measurementdata,thickness=None,unc=50e6,steps=10):
        self.mdata=measurementdata
        l,n=self.mdata.estimateLDavid()        
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
        print "A priori estimated Width: " + str(l*1e6)
        print "A priori estimated n: " + str(n)
        self.no_echos=self.calculate_no_echos(measurementdata.fdsam._tdData)
        self.outfilename=self.getFilenameSuggestion()
 
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
        fmax=self.mdata.fdref.getmaxfreq()
        
 #       bnds=(self.thicknessestimate-2*self.thicknessstd,self.thicknessestimate+3*self.thicknessstd)
#        print bnds        
        #overthink once more the cropping length
        n=self.n_estimated
        #oscillation period can be calculated from H data!
        f_span=c/2*1/self.l_estimated*1/n #(this should be the length of the oscillation)
        
        calcdata=copy.deepcopy(self.mdata)      
        calcdata.manipulateFDData(-1,[fmax-f_span*1,fmax+f_span*7])
        
        H_small=calcdata.H    
        py.figure(33)
        t=minimize(self.errorL,self.userthickness,args=((py.asarray(H_small),)),\
        method='Nelder-Mead', options={'xtol':1e-6,'disp': False})#, 'disp': False})
        return t.x[0]
        
    def errorL(self,l,H):
        n_small=[]
        n_small.append(self.calculaten(H,l))
        qs=self.QuasiSpace(n_small,H[1,0]-H[0,0],l)
        tv=self.totalVariation(n_small)
        py.plot(l,qs[0],'+')
        py.plot(l,tv[0],'*')        
        print "Currently evaluating length: "+ str(l[0]*1e6) + " TV Value " + str(tv[0])
        return qs[0]
    
#    def findL(self):
#        #1) cut frequency bands to narrow range around maxima
#        #take care to dont get the dc component
#        i=1     
#        fmax=0
#        while fmax<100e9:
#            fmax=self.mdata.H[py.argmax(self.mdata.fdref[i:,1]),0]
#            i+=1
#        print fmax/1e9
#        #overthink once more the cropping length
#        ls=self.getLengthes()
#        f_span=c/2*1/self.thicknessestimate #(this should be the length of the oscillation)
#        print str(f_span/1e9) + " GHz oscillation expected" 
#        #f_span=100e9 #try to set it fix first        
#        H_small=self.mdata.getcroppedData(self.mdata.H,fmax-f_span/2,fmax+f_span*3)
#    
#        #2) calculate n for all ls of the smalldata
#        calcs=[]
#        for l in ls:
#            calcs.append(self.calculaten(H_small,l))
#        py.figure(22)
#        allqs=self.QuasiSpace(calcs,H_small[1,0]-H_small[0,0],ls)  
#        alltvs=self.totalVariation(calcs)
#        py.plot(ls,allqs-min(allqs),'*')
#        py.plot(ls,alltvs-min(alltvs))
#        py.legend(('Quasi space method','TV method'))
#        print "optimal l quasi space: " + str(ls[py.argmin(allqs)]*1e6)
#        print "optimal l Total Variation: " + str(ls[py.argmin(alltvs)]*1e6)
#        
#        return ls[py.argmin(alltvs)],calcs
        
    def QuasiSpace(self,ns,df,ls):
        allqs=[]
        #at the moment, the most easy method(everything else failed...)
        #get the mean absolute value of the complete spectrum. works best
        

        for i in range(len(ns)):
#            xvalues=py.fftfreq(len(ns[i]),df)*c/4

            QSr=py.fft(ns[i].real-py.mean(ns[i].real))
            QSi=py.fft(ns[i].imag-py.mean(ns[i].real))
            #naive cut:
            ix=range(3,(len(QSr)/2-1))
        
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
        lb_r=self.mdata.H[:,1].real-self.mdata.H[:,2]*f
        lb_i=self.mdata.H[:,1].imag-self.mdata.H[:,3]*f
        ub_r=self.mdata.H[:,1].real+self.mdata.H[:,2]*f
        ub_i=self.mdata.H[:,1].imag+self.mdata.H[:,3]*f
        
        #ix=all indices for which after smoothening n H is still inbetwen the bounds        
        ix=py.all([H_r>=lb_r,H_r<ub_r,H_i>=lb_i,H_i<ub_i],axis=0)
#        #dont have a goood idea at the moment, so manually:
        for i in range(len(n_smoothed)):
            if ix[i]==0:
                n_smoothed[i]=n[i]
        print "SVMAF changed the refractive index at " + str(sum(ix)) + " frequencies"
        return n_smoothed      
      
    def calculateinits(self,H,l):
        #crop the time domain data to the first pulse and apply than the calculation
#        refdata=THzTdData(tdData=self.mdata.fdref._tdData.getFirstPuls(5e-12,10e-12))
#        samdata=THzTdData(tdData=self.mdata.fdsam._tdData.getFirstPuls(10e-12,5e-12))
#        refdata.zeroPaddData(2000)
#        samdata.zeroPaddData(2000)
        
#        firstref=FdData(refdata)
#        firstsam=FdData(samdata)
#        
#        initTeraData=teradata(firstref,firstsam)
#        oldfreqaxis=H[:,0]
#      
#        intpH=interp1d(initTeraData.H[:,0],initTeraData.H[:,1:],axis=0)
#        newH=intpH(oldfreqaxis)
#        n=self.n_0-newH[:,3]/(2*py.pi*oldfreqaxis*l)*c
#        alpha=-c/(2*py.pi*oldfreqaxis*l)*py.log(abs(newH[:,0])*(n+1)**2/(4*n))
####  
        #why is this NOT working
        tester=0
        #naive approach:
        n=self.n_0-(H[:,4]+py.pi*tester)/(2*py.pi*H[:,0].real*l)*c
        
        alpha=-c/(2*py.pi*H[:,0]*l)*py.log(abs(H[:,1])*(n+1)**2/(4*n))
##  
        return n.real,alpha.real
    
    def plotInits(self,H,l,figurenumber=200):
        inits=self.calculateinits(H,l)
        py.figure(figurenumber)
        py.title('Initial Conditions')
        py.plot(H[:,0].real/1e9,inits[0])
        py.plot(H[:,0].real/1e9,-inits[1])
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
            
            t=minimize(self.error_func,[inits[0,i],inits[1,i]],args=(H[i,:],l), method='SLSQP',\
            bounds=bnds, options={'ftol':1e-9,'maxiter':2000, 'disp': False})
#            t=minimize(self.error_func,[inits[0,i],inits[1,i]],args=(H[i,:],l), method='Nelder-Mead',
#            options={'xtol': 1e-6,'disp':False})
            res.append(t.x[0]-1j*t.x[1])
            vals.append(t.fun)
        return py.asarray(res)
    
    def doCalculation(self,bool_findl=1,n_SVMAFS=5,bool_silent=0):
        if bool_findl:
            self.l_opt=self.findLintelli()

        print '\033[92m\033[1m' + '  Use Sample Thickness: ' + str(self.l_opt*1e6) + ' micro m ' + '\033[0m'


        n=self.calculaten(self.mdata.H,self.l_opt)
        n_smoothed=n
        i=0
        
        while i<n_SVMAFS:
            n_smoothed=self.SVMAF(self.mdata.H[:,0],n_smoothed,self.l_opt)
            i+=1

        self.n=py.column_stack((self.mdata.H[:,0],n,n_smoothed))        
        
        return self.n
        
    def getFilenameSuggestion(self):
        filenames=self.mdata.fdsam._tdData.getfilename()
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
        
        H_theory=self.H_theory(self.mdata.H[:,0].real,[self.n[:,1].real,self.n[:,1].imag],self.l_opt)        
        #built the variable that should be saved:        
        savetofile=py.column_stack((
        self.n[:,0].real,
        self.n[:,1].real,self.n[:,1].imag, #the real and imaginary part of n
        self.n[:,2].real,self.n[:,2].imag, #the real and imaginary part of the smoothed n
        self.mdata.H[:,1].real,self.mdata.H[:,1].imag,#H_measured
        self.mdata.H[:,2].real,self.mdata.H[:,3].real,#uncertainties
        abs(self.mdata.H[:,1]),self.mdata.H[:,4].real,#absH,ph H measured    
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
        ix=py.argmin(abs(freq-self.mdata.H[:,0]))
        n_i=py.linspace(0.0,0.05,int(resolution/2))
        n_r=py.linspace(1,3,resolution)
        N_R,N_I=py.meshgrid(n_r,n_i)
        E_fu=py.zeros((len(n_i),len(n_r)))
        for i in range(len(n_r)):
            for k in range(len(n_i)):
                E_fu[k,i]=self.error_func([n_r[i],n_i[k]],self.mdata.H[ix,:],l)
#        print E_fu[:,2]
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
    path='/home/jahndav/Dropbox/THz-Analysis/NEW INRIM Measurements/'
     
    if name=='siliInrim':
        
    #    the silicon sample measured at inrim
        thickness=330e-6 
        samfiles=glob.glob('./silicon/*sam*')
        reffiles=glob.glob('./silicon/*ref*') 
        mode='INRIM'
        
    elif name=='siliInrim1':
        thickness=330e-6
        samfiles=glob.glob('./INRIM1/*sam*.dat')
        reffiles=glob.glob('./INRIM1/*ref*.dat') 
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
        samfiles=glob.glob('./MarburgData/*_TeflonI-3*')
        reffiles=glob.glob('./MarburgData/*ref7*')+glob.glob('./MarburgData/*ref8*')
    elif name=='Teflon2':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=2117e-6
        samfiles=glob.glob('./MarburgData/*_TeflonI-2*')
        reffiles=glob.glob('./MarburgData/*ref6*')+glob.glob('./MarburgData/*ref7*')
    elif name=='Teflon1':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=2470e-6
        samfiles=glob.glob('./MarburgData/*_TeflonI-1*')
        reffiles=glob.glob('./MarburgData/*ref5*')+glob.glob('./MarburgData/*ref6*')
    elif name=='PP1':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=2039e-6
        samfiles=glob.glob('./MarburgData/*_PPI-1*')
        reffiles=glob.glob('./MarburgData/*ref3*')+glob.glob('./MarburgData/*ref4*')
    elif name=='PP2':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=1294e-6
        samfiles=glob.glob('./MarburgData/*_PPI-2*')
        reffiles=glob.glob('./MarburgData/*ref4*')+glob.glob('./MarburgData/*ref5*')
    elif name=='Lactose1':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=2166e-6
        samfiles=glob.glob('./MarburgData/*_Lact1*')
        reffiles=glob.glob('./MarburgData/*ref0*')+glob.glob('./MarburgData/*ref1*')
    elif name=='Lactose2':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=1611e-6
        samfiles=glob.glob('./MarburgData/*_Lact2*')
        reffiles=glob.glob('./MarburgData/*ref1*')+glob.glob('./MarburgData/*ref2*')
    elif name=='Lactose3':
 #   a Teflon sample, flat real(n) and imag(n) expected
        thickness=2376e-6
        samfiles=glob.glob('./MarburgData/*_Lact3*')
        reffiles=glob.glob('./MarburgData/*ref2*')+glob.glob('./MarburgData/*ref3*')        
        teralyzer=py.loadtxt('./MarburgData/L3.txt.csv',delimiter=',',usecols=(0,1,2),skiprows=1)
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
        samfiles=glob.glob('./rehi/Sample*')
        reffiles=glob.glob('./rehi/Reference*') 
        mode='lucastestformat'
        teralyzer=py.loadtxt('./rehi/Rehi_Teralyzer_OK.txt')
    
    return thickness,samfiles,reffiles,mode,teralyzer


    
if __name__=="__main__":
    #the Data is at the moment not allowed to be zero padded, 
    #it somehow changes the phase of the signal, and though the initial conditions
    #are shit!    
    

    #Load Parameters from getparams
    thickness,samfiles,reffiles,mode,teralyzer=getparams('Lactose3')

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

#    reftd.doPlotWithunc()
#    samtd.doPlotWithunc()

#    #initialize the fd_data objects        
    ref_fd=FdData(reftd)
    sam_fd=FdData(samtd)

    ref_fd.doPlot()
    sam_fd.doPlot()
##    #initialize the mdata object (H,and so on)
#    mdata=teradata(ref_fd,sam_fd)
#    mdata.doPlots()
#    mdata.manipulateFDData(-11e9,[200e9,2.2e12])
#    l3=mdata.findAbsorptionLines()
 
#    myana=teralyz(mdata,thickness-30e-6,0.5*thickness,30)
#    myana.doCalculation()
#    myana.plotRefractiveIndex(1,1)
#    myana.saveResults()
   