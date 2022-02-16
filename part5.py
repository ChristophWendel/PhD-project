## ------------------------------------------------------------------------------------     o   o
## Electromagnetic cascades                                                                   I
##                                                                                          \.../
## Numerical solution of the kinetic equation 
## ------------------------------------------------------------------------------------

import os
import multiprocessing

from part4 import * # Import of the file.


## Considering the first term of equation 1

# The first term on the right-hand side of equation 1 is the above UsedDotNi(ga).


## Considering the second term of equation 1

def SecondTermOfEq1(ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,RelativeError=integrate.quad.__defaults__[3]):
    # The absolute value of the second term on the right-hand side is NElectrons*COfA21 with COfA21 from part 1. According to my analysis, it describes the spectral decrease-rate of electrons, i. e. the number-density of electrons per time-interval and per ga-interval. In the case of IncludeSynchrotron==True, the SyncPs are included as soft photons, in the case of IncludeSynchrotron==False, only the soft background photons are included (and one has to set nSyncPs=None and ValuesForxSyncSuperIt=None). It is in units of 1/(s*m^3).
    return NElectrons(ga)*COfA21(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def ElectronSpectralLossDensityRate(ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,ValuesForga=None,RelativeError=integrate.quad.__defaults__[3]):
    # With inclusion of electron-escape and / or synchrotron-energy-losses, the absolute value of the second term on the right-hand side is generalised to NElectrons*ElectronLossRate. It describes the spectral loss(-density)-rate of electrons, i. e. the number-density of electrons per time-interval and per ga-interval. Again, in the case of IncludeSynchrotron==True, the SyncPs are included as soft photons for IC-scattering and in the case of IncludeSynchrotron==False, only the soft background photons are included (and one has to set nSyncPs=None and ValuesForxSyncSuperIt=None). It is in units of 1/(s*m^3).
    return NElectrons(ga)*ElectronLossRate(ga,n0,nSyncPs,ValuesForxSyncSuperIt,ValuesForga,RelativeError)

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    def EvaluateSecondTermOfEq1():
        ValuesForgammaLog = np.logspace(2.0,5.0,21)  # A sampling-range for ga. ga ranges over reasonable values. Logarithmic division of the interval is used.
        print('Sample values for gamma:\n', ValuesForgammaLog)
        if not(IncludeSynchrotron):
            ValuesForSecondTermOfEq1 = [SecondTermOfEq1(i,Usedn0,None,None,CurrentNElectrons) for i in ValuesForgammaLog]
            ValuesForElectronSpectralLossDensityRate = [ElectronSpectralLossDensityRate(i,Usedn0,None,None,CurrentNElectrons) for i in ValuesForgammaLog]
        elif IncludeSynchrotron:
            ValuesForSecondTermOfEq1 = [SecondTermOfEq1(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons) for i in ValuesForgammaLog]
            ValuesForElectronSpectralLossDensityRate = [ElectronSpectralLossDensityRate(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons) for i in ValuesForgammaLog]
        print('Corresponding values for SecondTermOfEq1:\n', ValuesForSecondTermOfEq1)
        pl.figure(figsize=(12, 9), num="Cascade-equation: Second term (electrons' spectral loss-density-rate) versus original electron energy")
        pl.loglog(ValuesForgammaLog,ValuesForSecondTermOfEq1, label = 'Only via IC-down-scattering (second term)')
        pl.loglog(ValuesForgammaLog,ValuesForElectronSpectralLossDensityRate, label = 'Total loss (IC-down-scattering + escape + synchrotron)')
        pl.legend(loc="best", fontsize=14)
        pl.xlabel("Original electron energy $\gamma$")
        pl.ylabel("Absolute value of spectral number-density-rate in s$^{-1} \cdot $m$^{-3}$")


## Considering the third term of equation 1

# The third term on the right-hand side is similar to COfA21 from part 1 but has three features. Firstly, inside the integral, NElectrons is included. Secondly, it is integrated over the first argument of COfA1. Thus, one mustn't pull NElectrons out of the integral. Thirdly, the integration-range is going from ga to infinity. 
# According to my analysis, it describes the spectral density of the production-rate of electrons, i. e. the number-density of electrons per time-interval and per ga-interval. It is in units of 1/(s*m^3). In the following, the integrand is defined first, and afterwards the integration is done. 

def IntegrandOfThirdTermOfEq1(ga,gaP,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,RelativeError=integrate.quad.__defaults__[3]):      # In the following integral it is to be integrated over ga, so this has to be the first argument. 
    return NElectrons(ga)*COfA1(ga,gaP,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def BiggestIntegrationBorderOfThirdTermOfEq1(ValuesForxSyncSuperIt,gaP,NElectronsga0):
    '''Define the upper integration boundary of the following integral. According to Zdz. one should integrate from gaP to infinity. But integration up to infinity is not necessary as was shown. So, the upper border can be reduced. NElectronsga0 is an obvious cut-off. Moreover, for the case gaP<1/(4*x0), GammaLimit(x0,gaP) is used as an additional upper boundary, where GammaLimit is specified in part 1. In the synchrotron-case, there is an additional GammaLimit(xSync0Used,gaP) and the maximum of both these GammaLimits has to be used.
    Obsolete comment: This upper border could be an input value. This might be a critical parameter! It should be as big as possible, if NElectrons=1, because in this case, the high-ga part of the integrand is non-vanishing and approximately asymptotic to ga**(-1), which means that it is not negligible during the integration.  For real NElectrons, the decrease of the integrand is strengthened by the decrease of NElectrons. Thus it decreases more rapidly than ga**(-1) at high enough ga and hence it is negligible in the integration. This was confirmed via NElectronsGaussian at about 3*10**4 and similar width. For n0Planck, n0PL and n0Exp with x0=0.001, integration warnings are raised, if this value is bigger than about 10**17. For x0=0.005 or smaller, 10**18 is okay.
    Comment from version 5: I could show that for StartingNElectrons=NElectronsRecommended, the iterated NElectrons always has an upper cut-off at NElectronsga0. Hence, this can be taken as biggest integration-border.
    Comment of v1_1: The min was added. Also in this case, the upper cut-off of NElectrons is constraining the integration range.
    Comment of v2_4: AuxiliaryBiggestIntegrationBorderOfThirdTermOfEq1 was removed. The argument NElectronsga0 was introduced as well as ValuesForxSyncSuperIt'''
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        if gaP>=1.0/(4.0*ValuesForxSyncSuperIt[-1]): # ValuesForxSyncSuperIt[-1] is the current xSync0Used.
            return NElectronsga0
        else:
            return min(NElectronsga0,max(GammaLimit(x0,gaP),GammaLimit(ValuesForxSyncSuperIt[-1],gaP))) # Actually it might be max(min(NElectronsga0,GammaLimit(x0,gaP)),min(NElectronsga0,GammaLimit(xSync0Used,gaP))). This max() can be used because IntegrandOfA1 and IntegrandOfA1WithSynchrotron are defined such that they are =0 outside of their definition-range. Hence, from min() to max() one does not integrate over an unwanted contribution. It is however max(min(A,B),min(A,C)) = min(A,max(B,C)), which simplifies it.
    else: # The case without synchrotron.
        if gaP>=1.0/(4.0*x0):
            return NElectronsga0
        else:
            return min(NElectronsga0,GammaLimit(x0,gaP))

# Now integrate over the integrand of the third term of equation 1:
def ThirdTermOfEq1(gaP,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is in units of 1/(s*m^3). It describes the spectral density of the probability-rate for IC-scattering events of an electron-population with any energy bigger than gaP to final energy gaP off a soft external photon background of spectral number-density n0 (in case IncludeSynchrotron==False, one has to set nSyncPs=None and ValuesForxSyncSuperIt=None) or off the external background photons plus the SyncPs (in the case IncludeSynchrotron==True). The spectral number-density of the electrons is NElectrons. Thus it is a probability per unit time dt and per unit volume.'''
    LowestIntegrationBorder=gaP # Obsolete comment: For the lower border, gaP+0.0001*gaP is used instead of just gaP, as the integrand is not a number at ga=gaP. Now, the used lower border is always 0.01% above the value, that should actually be used. Thus, the error is small. Comment of v1_1: It was reduced to gaP, because in the Thomson-regime the upper integration border gets smaller and smaller and thus the integration range gets smaller and smaller. The problem with the integrand being nan at ga=gaP isn't a problem any more, due to the separate handling of this case in IntegrandOfA1.
    BiggestIntegrationBorder = BiggestIntegrationBorderOfThirdTermOfEq1(ValuesForxSyncSuperIt,gaP,LastValuesForgaIt[-1]) # Actually, ValuesForxSyncSuperIt[-1] is the current xSync0Used but the slice is taken only inside BiggestIntegrationBorderOfThirdTermOfEq1. LastValuesForgaIt[-1] is the current NElectronsga0.
    if LowestIntegrationBorder>=BiggestIntegrationBorder:
        return 0.0
    # Now the integration: 
    else:
        # Consider the case of the zeroth iteration, i. e. the case where CurrentNElectrons=NElectronsRecommended:
        if IterationCounter==0:
            NumberOfBorders = max(3,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
            ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
            # If a function with a Heaviside-function is used as NElectrons, the integration algorithm has problems at the jumping points. Especially in the case NElectronsRecommended, there are jumping points at Injectedga0 and at Injectedga1. Thus, these points are included into the list of integration-borders if they are situated in the interior of the integration-range.
            if LowestIntegrationBorder<Injectedga0<BiggestIntegrationBorder:
                ListOfIntegrationBorders = np.concatenate((ListOfIntegrationBorders,[Injectedga0]))
            if LowestIntegrationBorder<Injectedga1<BiggestIntegrationBorder:
                ListOfIntegrationBorders = np.concatenate((ListOfIntegrationBorders,[Injectedga1]))
        # Consider the cases after the zeroth iteration. Here CurrentNElectrons is an interpolated object, which has kinks at each item of that LastValuesForgaIt, that was used for the interpolation. So, take every kink as an element of ListOfIntegrationBorders if this element is inside the integration-range. As LastValuesForgaIt can change during the iteration, LastValuesForgaIt has to be an argument.
        elif IterationCounter:
            ListOfIntegrationBorders=np.asarray([LowestIntegrationBorder])
            for PotentialBorder in LastValuesForgaIt:
                if LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder:
                    ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
            if BiggestIntegrationBorder not in ListOfIntegrationBorders: # This if-test was added in v1_1. It is necessary for the case gaP<1/(4*x0) because here BiggestIntegrationBorder can be below NElectronsga0 and in the extreme case it is only very slightly above LowestIntegrationBorder.
                ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,BiggestIntegrationBorder)
        PotentialAdditionalBorders = [] # The integrand (the function COfA1) has kinks, due to additional contributions to COfA1. These additional contributions are caused by discontinuities along x. They should be included into the list of borders.
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[GammaLimit(ValuesForxSyncSuperIt[-1],gaP)])
        if Usedn0==n0MultiDelta:
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[GammaLimit(i,gaP) for i in x0MultiDelta])
        else:
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[GammaLimit(x0,gaP)])
        for PotentialBorder in PotentialAdditionalBorders: # Include these points into the list of integration-borders, if they are in the interior of the range.
            if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
                ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)        
        ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
        #print('    Integration borders of ThirdTermOfEq1: ', ListOfIntegrationBorders)
        ThirdTermOfEq1Result = 0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            ThirdTermOfEq1Result += integrate.quad(IntegrandOfThirdTermOfEq1, LeftBorder, RightBorder, args=(gaP,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,RelativeError), epsrel=RelativeError)[0]
        return ThirdTermOfEq1Result
    # Comment of version 9: IterationCounter was added to the function arguments.
    # Comment of version 10: LastValuesForgaIt was added to the function arguments.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    def EvaluateIntegrandOfThirdTermOfEq1(): # Look at the integrand:
        gammaP=0.1/x0
        LowestIntegrationBorder=gammaP
        BiggestIntegrationBorder=BiggestIntegrationBorderOfThirdTermOfEq1(ValuesForxSyncSuperIt,gammaP,CurrentLastValuesForgaIt[-1])*0.9999999999 # IntegrandOfThirdTermOfEq1 is = 0 at ga=BiggestIntegrationBorderOfThirdTermOfEq1(x0,gammaP). Hence, in logarithmic plotting this should be excluded.
        ValuesForgammaLog = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),70000)  # A sampling-range for ga. ga ranges from gaP to the upper integration border. Logarithmic division of the interval is used.
        pl.figure(figsize=(12, 9), num="Cascade-equation: Integrand of third term versus original electron energy")
        if not(IncludeSynchrotron):
            ValuesForIntegrandOfThirdTermOfEq1 = [IntegrandOfThirdTermOfEq1(i,gammaP,Usedn0,None,None,CurrentNElectrons) for i in ValuesForgammaLog]
            pl.scatter(LowestIntegrationBorder*x0, IntegrandOfThirdTermOfEq1(LowestIntegrationBorder,gammaP,Usedn0,None,None,CurrentNElectrons), 20, color='red')
            pl.scatter(BiggestIntegrationBorder*x0, IntegrandOfThirdTermOfEq1(BiggestIntegrationBorder,gammaP,Usedn0,None,None,CurrentNElectrons), 20, color='red')
        elif IncludeSynchrotron:
            ValuesForIntegrandOfThirdTermOfEq1 = [IntegrandOfThirdTermOfEq1(i,gammaP,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons) for i in ValuesForgammaLog]
            pl.scatter(LowestIntegrationBorder*x0, IntegrandOfThirdTermOfEq1(LowestIntegrationBorder,gammaP,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons), 20, color='red')
            pl.scatter(BiggestIntegrationBorder*x0, IntegrandOfThirdTermOfEq1(BiggestIntegrationBorder,gammaP,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons), 20, color='red')
        pl.loglog(ValuesForgammaLog*x0,ValuesForIntegrandOfThirdTermOfEq1, label="$\gamma' \cdot x_0 = %s$" % (gammaP*x0))# An error is raised here for reasons of scaling, if gammaP is slightly near to a power of 10.
        pl.legend(loc="best", fontsize=13)
        pl.xlabel("Original electron energy $\gamma$ times $x_0$")
        pl.ylabel("Integrand of third term") 
        pl.savefig("Cascade-equation - Integrand of third term versus original electron energy from %s.svg" % GetCurrentDate())
    # For values of ga smaller than gammaP, the integrand is not a number (which comes from COfA1 being not a number here). This is logical because for a certain gaP, ga has to be bigger than gaP, due to energy conservation. Thus, in the evaluation-function ValuesForgammaLog was prevented to contain values smaller than gaP.
    # Additionally, I have analysed the following: For the case gaP<1/(4*x0), it must be ga<GammaLimit(x0,gaP), where GammaLimit is specified in part 1. (In this case, the integrand is of value -0.0 at ga above GammaLimit(x0,gaP).) Attention has to be paid to this during the integration. For gaP>=1/(4*x0), ga has no upper boundary.
    # For the four different n0, the course of the integrand is qualitatively very similar. For ga>>gaP, it always decreases approximately like x**(-1). For gaP approx 1/(4*x0) it changes its course. For gaP <= 1/(4*x0) the value at ga=gaP is not a number (probably infinity). 
    # If one now pays attention to the fact that it was NElectrons(ga)=1 here but that in reality NElectrons(ga) is vanishing for big enough ga, one can infer that the integrand is vanishing for big enough ga, too. This fact makes it possible not having to integrate up to infinity in the integral of the third term.
    
    def EvaluateThirdTermOfEq1():
        #ValuesForgammaPLog = np.logspace(2,6,31)  # A sampling-range for gaP for fast evaluation. gaP ranges over reasonable values.
        FirstDivisionPoint=np.log10(0.1*ThomsonKNBoundary)
        SecondDivisionPoint=np.log10(100*ThomsonKNBoundary)
        ValuesForgammaPLog = np.concatenate((np.logspace(np.log10(Lowestga),FirstDivisionPoint,300,endpoint=False),np.logspace(FirstDivisionPoint,SecondDivisionPoint,200,endpoint=False),np.logspace(SecondDivisionPoint,np.log10(NElectronsga0),300)))  # A sampling-range for gaP for smooth plots. Logarithmic division is used and around 1/(4*x0) the spacing between sampling points is reduced, because the function changes rapidly here.
        pl.figure(figsize=(12, 9), num="Cascade-equation: Third term versus final electron energy in the case n_0(x)=%s" % Usedn0.__name__)
        if not(IncludeSynchrotron):
            ValuesForThirdTermOfEq1 = [ThirdTermOfEq1(i,Usedn0,None,None,CurrentNElectrons,1,CurrentLastValuesForgaIt) for i in ValuesForgammaPLog]
            pl.scatter(1/(8), ThirdTermOfEq1(1/(8*x0),Usedn0,None,None,CurrentNElectrons,1,CurrentLastValuesForgaIt), 20, color='red')
        elif IncludeSynchrotron:
            ValuesForThirdTermOfEq1 = [ThirdTermOfEq1(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons,1,CurrentLastValuesForgaIt) for i in ValuesForgammaPLog]
            pl.scatter(1/(8), ThirdTermOfEq1(1/(8*x0),Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons,1,CurrentLastValuesForgaIt), 20, color='red')
        pl.loglog(ValuesForgammaPLog*x0,ValuesForThirdTermOfEq1)        
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        #pl.legend(loc="best", fontsize=13)
        pl.xlabel("Final electron energy $\gamma' \cdot x_0$", fontsize=16)
        pl.ylabel("Third term in s$^{-1} \cdot $m$^{-3}$", fontsize=16)
        pl.savefig("Cascade-equation - Third term versus final electron energy from %s.svg" % GetCurrentDate())
    # At gaP=1/(4*x0) the third term always has a kind of step, which might come from the change of the shape of the integrand at this value (and perhaps from the change of the integration-range). For reasonable NElectrons, the third term is vanishing for high enough gaP.


## Considering the brackets in the fourth term of equation 1:

def IntegrandOfBracketsOfEq1(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,RelativeError=integrate.quad.__defaults__[3]):      # This is the integrand of the integral within the brackets. In the integral it is to be integrated over ga, so this has to be the first argument. The units of this quantity are 1/(s*m^3).
    return NElectrons(ga)*COfA1DependingOnxga(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def IntegralOfBracketsOfEq1(xga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is in units of 1/(s*m^3). It describes the spectral density of the production-rate of photons of final energy xga via IC-scattering events of an electron-population off soft photons. The spectral number-density of the electrons is NElectrons. Thus it is a probability per unit time dt and per unit volume. Again, in the case of IncludeSynchrotron==True, the SyncPs and the soft external background photons are included as soft photons for IC-scattering and in the case of IncludeSynchrotron==False, only the external photons are included (and one has to set nSyncPs=None and ValuesForxSyncSuperIt=None).'''
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        LowestIntegrationBorder=min(GammaMinOnxga(x0,xga),GammaMinOnxga(ValuesForxSyncSuperIt[-1],xga)) # ValuesForxSyncSuperIt[-1] is the current xSync0Used. Similarly to ThirdTermOfEq1, the min() may be used here, because the function SyncP-contribution as well as the ExPs-contribution to COfA1DependingOnxga are defined such that they vanish outside of their respective definition range.
    else: # This case has to be treated separately, because ValuesForxSyncSuperIt has the value None, which cannot be sliced.
        LowestIntegrationBorder=GammaMinOnxga(x0,xga)
    BiggestIntegrationBorder=LastValuesForgaIt[-1] # LastValuesForgaIt[-1] is the current NElectronsga0. Obsolete comment: This could be an input value. The upper integration border of this integral actually should be infinity. But for the case NElectrons=1 the integrand is already substantially steeper than ga**(-1). Thus, integration up to infinity is not necessary, all the more if NElectrons is decreasing with high ga, too. Comment from version 5: I could show that for StartingNElectrons=NElectronsRecommended, the iterated NElectrons always has an upper cut-off at NElectronsga0. Hence, this can be taken as biggest integration-border.
    if LowestIntegrationBorder>=BiggestIntegrationBorder:
        return 0.0
    # Now the integration: 
    else:
        # Consider the case of the zeroth iteration, i. e. the case where CurrentNElectrons=NElectronsRecommended:
        if IterationCounter==0:
            NumberOfBorders = max(3,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
            ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
            # If a function with a Heaviside-function is used as NElectrons, the integration algorithm may get problems at the jumping points. Especially in the case NElectronsRecommended, there are jumping points at Injectedga0 and at Injectedga1. Thus, these points are included into the list of integration-borders if they are situated in the interior of the integration-range.
            if LowestIntegrationBorder<Injectedga0<BiggestIntegrationBorder:
                ListOfIntegrationBorders = np.concatenate((ListOfIntegrationBorders,[Injectedga0]))
            if LowestIntegrationBorder<Injectedga1<BiggestIntegrationBorder:
                ListOfIntegrationBorders = np.concatenate((ListOfIntegrationBorders,[Injectedga1]))
        # Consider the cases after the zeroth iteration. Here CurrentNElectrons is an interpolated object, which has kinks at each item of the LastValuesForgaIt, that was used for the interpolation. So, take every kink as an element of ListOfIntegrationBorders if this element is inside the integration-range. As LastValuesForgaIt can change during the iteration, LastValuesForgaIt has to be an argument.
        elif IterationCounter:
            ListOfIntegrationBorders=[LowestIntegrationBorder]
            for PotentialBorder in LastValuesForgaIt:
                if LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder:
                    ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
        # for PotentialBorder in [1.0*xga+0.001*10**7,1.0*xga+0.003*10**7,1.0*xga+0.004*10**7,1.0*xga+0.0045*10**7,1.0*xga+0.006*10**7]: # It was found that, for n0Planck and high values of xga, COfA1DependingOnxga and hence IntegrandOfBracketsOfEq1, too, has a sharp maximum here. Comment of v2_5: It is not sure whether these additional borders are reasonable any more at all. Therefore they are commented out.
        #     if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
        #         ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
        ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
        PotentialAdditionalBorders = [] # The integrand can have kinks, due to additional contributions to COfA1. These additional contributions are caused by discontinuities along x.
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[GammaMinOnxga(ValuesForxSyncSuperIt[-1],xga)])
        if Usedn0==n0MultiDelta:
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[GammaMinOnxga(i,xga) for i in x0MultiDelta])
        else:
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[GammaMinOnxga(x0,xga)])
        
        for PotentialBorder in PotentialAdditionalBorders: # Include these points into the list of integration-borders, if they are in the interior of the range.
            if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
                ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)        
        ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
        # print('Integration borders of IntegralOfBracketsOfEq1: ', ListOfIntegrationBorders)
        IntegralOfBracketsOfEq1Result = 0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            IntegralOfBracketsOfEq1Result += integrate.quad(IntegrandOfBracketsOfEq1, LeftBorder, RightBorder, args=(xga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,RelativeError), epsrel=RelativeError)[0]
        return IntegralOfBracketsOfEq1Result
    # Comment of version 9: IterationCounter was added to the function arguments.
    # Comment of version 10: LastValuesForgaIt was added to the function arguments.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    
    # IntegrandOfBracketsOfEq1 can be plotted in dependence on ga for various values of xga:
    def EvaluateIntegrandOfBracketsOfEq1ForCertainxga():
        ValuesForgaLog = np.logspace(np.log10(GammaMinOnxga(x0,xgamma)+0.1),8.0,2000)  # A sampling-range for ga. Actually ga ranges from GammaMinOnxga(x0,xga) to infinity. Logarithmic division of the interval is used.
        ValuesForgaLin = np.linspace(GammaMinOnxga(x0,xgamma)+0.1,10.0**6.0,2000)  # A sampling-range for ga. Actually ga ranges from GammaMinOnxga(x0,xga) to infinity. Linear division of the interval is used here.
        if not(IncludeSynchrotron):
            ValuesForIntegrandOfBracketsOfEq1ForCertainxga = [IntegrandOfBracketsOfEq1(i,xgamma,Usedn0,None,None,CurrentNElectrons) for i in ValuesForgaLog]
        elif IncludeSynchrotron:
            ValuesForIntegrandOfBracketsOfEq1ForCertainxga = [IntegrandOfBracketsOfEq1(i,xgamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons) for i in ValuesForgaLog]
        pl.figure(figsize=(12, 9), num="Cascade-equation: Integrand within the brackets versus gamma")
        pl.loglog(ValuesForgaLog,ValuesForIntegrandOfBracketsOfEq1ForCertainxga, label=LabelForFixedxgamma)
        pl.legend(loc="best")
        pl.xlabel("Electron energy $\gamma$")
        pl.ylabel("Integrand within the brackets $C(\gamma,\gamma-x_{\gamma}) \cdot N(\gamma)$ in s$^{-1}$m$^{-3}$")
    
    def EvaluateIntegralOfBracketsOfEq1():
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            MinimumxgaLog=np.log10(min(xgammaLimitClean(x0,gamma),xgammaLimitClean(ValuesForxSyncSuperIt[-1],gamma))) # This is the lower integration border of 4th term.
        else: # This case has to be treated separately, because ValuesForxSyncSuperIt has the value None, which cannot be sliced.
            MinimumxgaLog=np.log10(xgammaLimitClean(x0,gamma))
        MaximumxgaLog = np.log10(DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,CurrentLastValuesForgaIt[-1]))
        FirstDivisionxgaLog = 1.0005*MinimumxgaLog
        SecondDivisionxgaLog = 0.9980*MaximumxgaLog
        ValuesForxgaLog = np.concatenate((np.logspace(MinimumxgaLog,FirstDivisionxgaLog,100,endpoint=False),np.logspace(FirstDivisionxgaLog,SecondDivisionxgaLog,300,endpoint=False),np.logspace(SecondDivisionxgaLog,MaximumxgaLog,100)))  # A sampling-range for xga.
        if not(IncludeSynchrotron):
            ValuesForIntegralOfBracketsOfEq1 = [IntegralOfBracketsOfEq1(i,Usedn0,None,None,CurrentNElectrons,1,CurrentLastValuesForgaIt) for i in ValuesForxgaLog]
        elif IncludeSynchrotron:
            ValuesForIntegralOfBracketsOfEq1 = [IntegralOfBracketsOfEq1(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons,1,CurrentLastValuesForgaIt) for i in ValuesForxgaLog]
        pl.figure(figsize=(12, 9), num="Cascade-equation: Integral within the brackets versus final photon energy")
        pl.loglog(ValuesForxgaLog,ValuesForIntegralOfBracketsOfEq1)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        pl.xlabel("Final photon energy $x_\gamma$", fontsize=16)
        pl.ylabel("Integral within brackets in s$^{-1} \cdot $m$^{-3}$", fontsize=16)
        pl.savefig("Cascade-equation - Integral within brackets from %s.svg" % GetCurrentDate())
    # Concerning the case NElectrons=NElectronsRecommended: At a certain xga the lower integration border is equal to the upper jumping point of NElectronsRecommended, i. e. GammaMinOnxga(xga,x0)=Injectedga0. I could show that this happens at xga=Injectedga0/(1+1/(Injectedga0*4*x0)). For bigger xga, the lower integration border is above the upper jumping point. In this case, the integral is exactly vanishing. Thus, at xga=Injectedga0/(1+1/(Injectedga0*4*x0)) the integration of the fourth term may have problems.


## Considering the fourth term of equation 1:

# The fourth term of equation 1 of Zdz can now be divided into two summands. The first summand's integrand is pOfB18 times IntegralOfBracketsOfEq1 while the second summand's integrand is pOfB18 times UsedDotni. In the case where HE photon escape is included, pOfB18 is substituted for by NormalisedSpectralPPProbability. For the case IncludeSynchrotron==True, NormalisedSpectralPPProbability also includes the pair-production on SyncPs.

# According to my analysis this describes the double-spectral rate of pair-production events, that produce particles of energy ga via HE-photons, that were previously produced via IC-scattering. Thus, it is a rate per unit volume, per unit final particle energy dga and per HE-photon energy interval dxga. It is in units of 1/(s*m^3).
if IncludeHEPhotonEscape==False:
    def IntegrandOfFirstSummandOf4thTerm(xga,ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
        return pOfB18(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)*IntegralOfBracketsOfEq1(xga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError)
elif IncludeHEPhotonEscape==True:
    def IntegrandOfFirstSummandOf4thTerm(xga,ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
        return NormalisedSpectralPPProbability(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)*IntegralOfBracketsOfEq1(xga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError)
# Comment of version 9: IterationCounter was added to the function arguments.
# Comment of version 10: LastValuesForgaIt was added to the function arguments.

# According to my analysis this describes the double-spectral rate of pair-production events, that produce particles of energy ga via HE-photons, that were injected. Thus it is a rate per unit volume, per unit final particle energy dga and per HE-photon energy interval dxga. It is in units of 1/(s*m^3).
if IncludeHEPhotonEscape==False:
    def IntegrandOfSecondSummandOf4thTerm(xga,ga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError=integrate.quad.__defaults__[3]):
        return pOfB18(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)*Dotni(xga)
elif IncludeHEPhotonEscape==True:
    def IntegrandOfSecondSummandOf4thTerm(xga,ga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError=integrate.quad.__defaults__[3]):
        return NormalisedSpectralPPProbability(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)*Dotni(xga)

# According to my analysis this describes the double-spectral rate of pair-production events, that produce particles of energy ga via all HE-photons. Thus it is a rate per unit volume, per unit final particle energy dga and per HE-photon energy interval dxga. It is in units of 1/(s*m^3).
# In order to treat this integrand in only one integral, pOfB18 (NormalisedSpectralPPProbability) is drawn before the parentheses. Thus, pOfB18 (NormalisedSpectralPPProbability) has to be evaluated only once, which is an advantage for time-efficiency. But the second summand is supposed to be integrated only in the range between HEPhotonsxga1 and HEPhotonsxga0 (as Dotni is non-vanishing only here). Therefore, the Heaviside-function is included in the second term.
if IncludeHEPhotonEscape==False:
    def CompleteIntegrandOf4thTerm(xga,ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):  
        return pOfB18(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)*(IntegralOfBracketsOfEq1(xga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError)+Dotni(xga)*Heavi(HEPhotonsxga1,xga,HEPhotonsxga0))
elif IncludeHEPhotonEscape==True:
    def CompleteIntegrandOf4thTerm(xga,ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):  
        return NormalisedSpectralPPProbability(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)*(IntegralOfBracketsOfEq1(xga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError)+Dotni(xga)*Heavi(HEPhotonsxga1,xga,HEPhotonsxga0))
# Comment of version 9: IterationCounter was added to the function arguments.
# Comment of version 10: LastValuesForgaIt was added to the function arguments. 

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    def EvaluateSampleIntegrationBorders():
        BiggestIntegrationBorder = 4.0*10**10#BiggestIntegrationBorderOf4thTermOfEq1
        LowestIntegrationBorder = np.logspace(np.log10(5.0*10**4), np.log10(BiggestIntegrationBorder),600, endpoint=False)
        ListsOfIntegrationBorders = np.asarray([UsedSampleIntBordersIn4thTerm(i,BiggestIntegrationBorder)[0] for i in LowestIntegrationBorder])
        pl.figure(figsize=(12, 10), num="Integration borders according to %s" % UsedSampleIntBordersIn4thTerm.__name__)
        for i in range(len(LowestIntegrationBorder)):
            ListOfIntegrationBorders = ListsOfIntegrationBorders[i]
            for j in range(len(ListOfIntegrationBorders)):
                pl.loglog(LowestIntegrationBorder[i],ListOfIntegrationBorders[j], marker='.', color='b', markersize=2)
        #pl.text(BiggestIntegrationBorder/100,10*ListsOfIntegrationBorders[0][0],r'%s' % UsedSampleIntBordersIn4thTerm.__name__, fontsize=20, color='b')
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=28, direction='inout', width=3, length=10, pad=7)
        ax.yaxis.set_tick_params(labelsize=28, direction='inout', width=3, length=10)
        pl.xlabel("Lowest integration border $b_{\mathrm{l}}$", fontsize=34)
        pl.ylabel("Sampled integration borders $\{b_{\mathrm{l}},\,\ldots,\,b_{n+1}\}\,\,\,\,\,$  ", fontsize=34)
        pl.ylim(ListsOfIntegrationBorders[0][0],2*BiggestIntegrationBorder)
        pl.xlim(LowestIntegrationBorder[0],2*BiggestIntegrationBorder)
        pl.subplots_adjust(top=0.97, bottom=0.14, left=0.15, right=0.97)
        pl.savefig("Run %s, %s from %s.pdf" % (RunIdentifier,UsedSampleIntBordersIn4thTerm.__name__,GetCurrentDate()))

def AuxiliaryIntegrateIntegrandOf4thTerm(Integrand, LeftBorder, RightBorder, ga, n0, ValuesForxSyncSuperIt, ValuesFornSyncPs, Dotni, IterationCounter, LastValuesForgaIt, RelativeError, NElectrons=None, ValuesForNextNElectrons=None):
    '''This will be used, if the FourthTerm is determined with help of multiprocessing. Then, this function just determines the integral of the respective Integrand from LeftBorder to RightBorder. If NElectrons is an interpolated object then it will be retrieved from LastValuesForgaIt and ValuesForNextNElectrons due to the pickling issue. Integrand is intended to be IntegrandOfFirstSummandOf4thTerm, IntegrandOfSecondSummandOf4thTerm or CompleteIntegrandOf4thTerm. ValuesForxSyncSuperIt and ValuesFornSyncPs have to be given. In the case without synchrotron-radiation, they should be None, in the case with synchrotron, they have to contain the value-arrays. From that, nSyncPs can be reconstructed via interpolation.'''
    ProcessID = multiprocessing.current_process().name+' AuxiliaryIntegrator'
    if NElectrons is not None and ValuesForNextNElectrons is None: # This case is tailored for NElectrons being a common function. To call this case, only specify NElectrons, but do not specify ValuesForNextNElectrons.
        print('%-47s:             No internal interpolation necessary.' % ProcessID)
    elif NElectrons is None and ValuesForNextNElectrons is not None: # This case is tailored for NElectrons being a non-picklable object, especially an interpolated object. In this case the interpolated object is reconstructed via LastValuesForgaIt and ValuesForNextNElectrons. To call this case, only specify ValuesForNextNElectrons, but do not specify NElectrons.
        #print('%-47s:             Doing internal interpolation.' % ProcessID)
        NElectrons = interp1d(LastValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
    else:
        raise WrongFunctionCallError(AuxiliaryIntegrateIntegrandOf4thTerm.__name__)
    if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray: # If these two objects (arrays) are given, then nSyncPs is reconstructed. This case might be equal to the case PerformanceMode=='SuperIteration' and IncludeSynchrotron==True.
        #print('%-47s:             Doing interpolation of nSyncPs.' % ProcessID)
        nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
    else:
        print('%-47s:             Using nSyncPs=None.' % ProcessID)
        nSyncPs = None
    StartingTimeOfAuxiliaryIntegration = time.time() # Measure the spent time of the auxiliary integration.
    if Integrand == IntegrandOfFirstSummandOf4thTerm:
        AuxiliaryResult = integrate.quad(Integrand, LeftBorder, RightBorder, args=(ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError), epsrel=RelativeError)[0]
    elif Integrand == IntegrandOfSecondSummandOf4thTerm:
        AuxiliaryResult = integrate.quad(Integrand, LeftBorder, RightBorder, args=(ga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError), epsrel=RelativeError)[0]
    elif Integrand == CompleteIntegrandOf4thTerm:
        AuxiliaryResult = integrate.quad(Integrand, LeftBorder, RightBorder, args=(ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,RelativeError), epsrel=RelativeError)[0]
    EndingTimeOfAuxiliaryIntegration = time.time() # Now, determine the point of time after the end of the auxiliary integration.
    TimeIntervalOfAuxiliaryIntegration = (EndingTimeOfAuxiliaryIntegration-StartingTimeOfAuxiliaryIntegration) # Determine the time in s, that was spent by the evaluation. 
    print('%-47s:             Time for aux. integr. from %12.1f to %12.1f is %5.1f s. Result: %s /(s*m^3)' % (ProcessID,LeftBorder,RightBorder,TimeIntervalOfAuxiliaryIntegration,AuxiliaryResult))
    return AuxiliaryResult

def FourthTermOfEq1Alternative1(ga,n0,ValuesForxSyncSuperIt,ValuesFornSyncPs,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,Include1stIntRangeOf4thTerm=True,ValuesForNextNElectrons=None,RelativeError=integrate.quad.__defaults__[3]):
    '''An auxiliary function for the computation of the fourth term. The most obvious computation of the integral in the case UsedDotni != DotniDelta or DotniZero: Here, the complete integrand is integrated at once.'''
    ProcessID = multiprocessing.current_process().name+' Computer of 4. term (alt. 1)'
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        LowestIntegrationBorder=min(xgammaLimitClean(x0,ga),xgammaLimitClean(ValuesForxSyncSuperIt[-1],ga)) # ValuesForxSyncSuperIt[-1] is the current xSync0Used. Again, the min() may be used here, because the SyncP-contribution as well as the ExPs-contribution to COfA1DependingOnxga are defined such that they vanish outside of their respective definition range.
    else: # This case has to be treated separately, because ValuesForxSyncSuperIt has the value None, which cannot be sliced.
        LowestIntegrationBorder=xgammaLimitClean(x0,ga)
    BiggestIntegrationBorder=DetermineBiggestIntegrationBorderOf4thTermOfEq1(ValuesForxSyncSuperIt,LastValuesForgaIt[-1]) # The slicing of ValuesForxSyncSuperIt is performed internally. LastValuesForgaIt[-1] is the current NElectronsga0.
    if LowestIntegrationBorder>=BiggestIntegrationBorder:
        return 0.0
    # Now the integration: 
    else:
        ListOfIntegrationBorders = np.asarray([LowestIntegrationBorder,BiggestIntegrationBorder]) # Initialise the sampling of the integration borders into subranges.
        # Additional borders can be necessary: There are discontinuities of the integrand at HEPhotonsxga1 and at HEPhotonsxga0, where the Heaviside-function jumps. So, choose these two points as additional borders if they are inside the integration interval. Furthermore, it was inferred from EvaluateIntegralOfBracketsOfEq1() and analytically in "1988ApJ...335..786Z - Determining an upper cut-off.png", too, that IntegralOfBracketsOfEq1 in dependence on xga has an upper cut-off at DetermineUpperCutOffOfEq5(xSync0Used,NElectronsga0).
        PotentialAdditionalBorders = [HEPhotonsxga1,HEPhotonsxga0,DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,LastValuesForgaIt[-1])] # Again, the slicing of ValuesForxSyncSuperIt is performed internally.
        # Moreover, the numerator (POfB1) of the function pOfB18, respectively NormalisedSpectralPPProbability for the inclusion of escape, has jumps (additional contributing components) at xga=xgammaLimitClean(i,ga) for i being a jump of the LEP-distributions. It might be useful to include these jump points as additional borders.
        # Moreover, the denominator (POfB11, respectively HEPhotonLossRate for the inclusion of escape) of the function pOfB18, respectively NormalisedSpectralPPProbability for the inclusion of escape, has jumps at 1/i, so include these, too:
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[xgammaLimitClean(ValuesForxSyncSuperIt[-1],ga),1/ValuesForxSyncSuperIt[-1]]) # ValuesForxSyncSuperIt[-1] is xSync0Used.
        if Usedn0==n0MultiDelta:
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[xgammaLimitClean(i,ga) for i in x0MultiDelta])
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[1/i for i in x0MultiDelta])
        else:
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[xgammaLimitClean(x0,ga),1/x0])
        if InjectionType=='ResultsOfOldIt': # In this case, UsedDotni has kinks, which have to be used as integration borders.
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,OldValuesForxgammanHEPhotons)
        for PotentialBorder in PotentialAdditionalBorders: # Include these points into the list of integration-borders, if they are in the interior of the range.
            if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
                ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)        
        ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
        # Now, subdivide each of the arising ranges and take the divisions as additional borders:
        AdditionalBorders = []
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            AdditionalBorders=np.append(AdditionalBorders,UsedSampleIntBordersIn4thTerm(LeftBorder,RightBorder)[0][1:-1])
        for PotentialBorder in AdditionalBorders: # Include these borders into the list of integration-borders, if they are not included yet.
            if PotentialBorder not in ListOfIntegrationBorders:
                ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
        ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
        print('%-47s:         Integration borders of 4. term:' % ProcessID, ListOfIntegrationBorders)
        if not(Include1stIntRangeOf4thTerm) and len(ListOfIntegrationBorders)>3: # The case where both the 1st integration range can be neglected and there are at least 3 ranges.
            ListOfIntegrationBorders=ListOfIntegrationBorders[1:] # Delete the first border, thus neglect the first range.
            print('%-47s:         Actually used integration borders:' % ProcessID, ListOfIntegrationBorders)
        NumberOfSubtasks = len(ListOfIntegrationBorders)-1 # The integration along one subrange is considered as one subtask.
        if UseMultiprocessingIn4thTermOfEq1==False: # No multiprocessing.
            print('%-47s:         Determine 4. term without multiprocessing.' % ProcessID)
            if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray: # If these two objects (arrays) are given, then nSyncPs is reconstructed. This case might be equal to the case PerformanceMode=='SuperIteration' and IncludeSynchrotron==True.
                #print('%-47s:         Doing interpolation of nSyncPs.' % ProcessID)
                nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
            else:
                print('%-47s:         Using nSyncPs=None.' % ProcessID)
                nSyncPs = None
            FourthTermOfEq1Result = 0.0
            for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                FourthTermOfEq1Result += integrate.quad(CompleteIntegrandOf4thTerm, LeftBorder, RightBorder, args=(ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,RelativeError), epsrel=RelativeError)[0]
        elif UseMultiprocessingIn4thTermOfEq1==True: # Multiprocessing is used, added in v1_4.
            NumberOfSubprocesses=min(NumberOfSubtasks,int(np.floor(multiprocessing.cpu_count()/(2*NumberOfIteratingProcesses)))) # multiprocessing.cpu_count()/2 gives the number of physical CPUs available on this machine. Division by NumberOfIteratingProcesses, gives the number of physical CPUs that are maximally available for iteration at one point of gamma.
            PoolOfWorkers = multiprocessing.Pool(processes=NumberOfSubprocesses) # Initialise multiple processes.
            if isinstance(NElectrons, types.FunctionType): # The case where NElectrons is a common picklable function.
                print('%-47s:         Determine 4. term for function-object with multiprocessing.' % ProcessID)
                FourthTermOfEq1ResultQueue = [PoolOfWorkers.apply_async(AuxiliaryIntegrateIntegrandOf4thTerm, args=[CompleteIntegrandOf4thTerm, LeftBorder, RightBorder, ga, n0, ValuesForxSyncSuperIt, ValuesFornSyncPs, Dotni, IterationCounter, LastValuesForgaIt, RelativeError], kwds={'NElectrons': NElectrons, 'ValuesForNextNElectrons': None}) for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:])]
            elif isinstance(NElectrons, interp1d): # Here, NElectrons has to be transferred to the child processes as arrays and interpolated in the child processes.
                print('%-47s:         Determine 4. term for interp1d-object with multiprocessing.' % ProcessID)
                FourthTermOfEq1ResultQueue = [PoolOfWorkers.apply_async(AuxiliaryIntegrateIntegrandOf4thTerm, args=[CompleteIntegrandOf4thTerm, LeftBorder, RightBorder, ga, n0, ValuesForxSyncSuperIt, ValuesFornSyncPs, Dotni, IterationCounter, LastValuesForgaIt, RelativeError], kwds={'NElectrons': None, 'ValuesForNextNElectrons': ValuesForNextNElectrons}) for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:])]
            FourthTermOfEq1Result = np.sum([r.get() for r in FourthTermOfEq1ResultQueue]) # Sum over the part-integrals
            PoolOfWorkers.close() # Clear the Pool-object to exit the processes and to...
            PoolOfWorkers.join() # ...retrieve the used RAM.
        print('%-47s:         Finished determining 4. term. Result: %s /(s*m^3)' % (ProcessID,FourthTermOfEq1Result))
        return FourthTermOfEq1Result
    # Comment of version 9: IterationCounter was added to the function arguments.
    # Comment of version 10: LastValuesForgaIt was added to the function arguments. 

def FourthTermOfEq1Alternative2(ga,n0,ValuesForxSyncSuperIt,ValuesFornSyncPs,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,Include1stIntRangeOf4thTerm=True,ValuesForNextNElectrons=None,RelativeError=integrate.quad.__defaults__[3]):
    '''An auxiliary function for the computation of the fourth term. An alternative for the computation of the integral in the case UsedDotni != DotniDelta or DotniZero: Here, the integral is split into two summands.'''
    ProcessID = multiprocessing.current_process().name+' Computer of 4. term (alt. 2)'
    if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray: # If these two objects (arrays) are given, then nSyncPs is reconstructed. This case might be equal to the case PerformanceMode=='SuperIteration' and IncludeSynchrotron==True.
        #print('%-47s:         Doing interpolation of nSyncPs.' % ProcessID)
        nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
    else:
        print('%-47s:         Using nSyncPs=None.' % ProcessID)
        nSyncPs = None
    # Integration of the first summand:
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        Part1LowestIntegrationBorder=min(xgammaLimitClean(x0,ga),xgammaLimitClean(ValuesForxSyncSuperIt[-1],ga)) # ValuesForxSyncSuperIt[-1] is the current xSync0Used. This is the physical lower border, which is also valid in part2 of this definition.
    else: # This case has to be treated separately, because ValuesForxSyncSuperIt has the value None, which cannot be sliced.
        Part1LowestIntegrationBorder=xgammaLimitClean(x0,ga)
    Part1BiggestIntegrationBorder=DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,LastValuesForgaIt[-1]) # This is the upper cut-off of the integral in brackets. The slicing of ValuesForxSyncSuperIt is performed internally.
    if Part1LowestIntegrationBorder>=Part1BiggestIntegrationBorder:
        Part14thTermOfEq1Result = 0.0
    # Now the integration: 
    else:
        Part1ListOfIntegrationBorders = np.asarray([Part1LowestIntegrationBorder,Part1BiggestIntegrationBorder]) # Initialise the sampling of the integration borders.
        Part1PotentialAdditionalBorders = [] # Additional borders can be necessary.
        #The numerator (POfB1) of the function pOfB18, respectively NormalisedSpectralPPProbability for the inclusion of escape, has jumps (additional contributing components) at xga=xgammaLimitClean(i,ga) for i being a jump of the LEP-distributions. It might be useful to include these jump points as additional borders.
        # Moreover, the denominator (POfB11, respectively HEPhotonLossRate for the inclusion of escape) of the function pOfB18, respectively NormalisedSpectralPPProbability for the inclusion of escape, has jumps at 1/i, so include these, too:        
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            Part1PotentialAdditionalBorders = np.append(Part1PotentialAdditionalBorders,[xgammaLimitClean(ValuesForxSyncSuperIt[-1],ga),1/ValuesForxSyncSuperIt[-1]]) # ValuesForxSyncSuperIt[-1] is xSync0Used.
        if Usedn0==n0MultiDelta:
            Part1PotentialAdditionalBorders = np.append(Part1PotentialAdditionalBorders,[xgammaLimitClean(i,ga) for i in x0MultiDelta])
            Part1PotentialAdditionalBorders = np.append(Part1PotentialAdditionalBorders,[1/i for i in x0MultiDelta])
        else:
            Part1PotentialAdditionalBorders = np.append(Part1PotentialAdditionalBorders,[xgammaLimitClean(x0,ga),1/x0])
        for PotentialBorder in Part1PotentialAdditionalBorders: # Include these points into the list of integration-borders, if they are in the interior of the range.
            if (Part1LowestIntegrationBorder<PotentialBorder<Part1BiggestIntegrationBorder) and (PotentialBorder not in Part1ListOfIntegrationBorders):
                Part1ListOfIntegrationBorders=np.append(Part1ListOfIntegrationBorders,PotentialBorder)
        Part1ListOfIntegrationBorders = np.sort(Part1ListOfIntegrationBorders)
        # Now, subdivide each of the arising ranges and take the divisions as additional borders:
        Part1AdditionalBorders = []
        for LeftBorder, RightBorder in zip(Part1ListOfIntegrationBorders[:-1], Part1ListOfIntegrationBorders[1:]):
            Part1AdditionalBorders=np.append(Part1AdditionalBorders,UsedSampleIntBordersIn4thTerm(LeftBorder,RightBorder)[0][1:-1])
        for PotentialBorder in Part1AdditionalBorders: # Include these borders into the list of integration-borders, if they are not included yet.
            if PotentialBorder not in Part1ListOfIntegrationBorders:
                Part1ListOfIntegrationBorders=np.append(Part1ListOfIntegrationBorders,PotentialBorder)
        Part1ListOfIntegrationBorders = np.sort(Part1ListOfIntegrationBorders)
        print('%-47s:         Integration borders of 1. summand of 4. term:' % ProcessID, Part1ListOfIntegrationBorders)
        if not(Include1stIntRangeOf4thTerm) and len(Part1ListOfIntegrationBorders)>3: # The case where both the 1st integration range can be neglected and there are at least 3 ranges.
            Part1ListOfIntegrationBorders=Part1ListOfIntegrationBorders[1:] # Delete the first border, thus neglect the first range.
            print('%-47s:         Actually used integration borders:' % ProcessID, Part1ListOfIntegrationBorders)
        Part1NumberOfSubtasks = len(Part1ListOfIntegrationBorders)-1
        if UseMultiprocessingIn4thTermOfEq1==False: # No multiprocessing.
            print('%-47s:         Determine 1. summand of 4. term without multiprocessing.' % ProcessID)
            Part14thTermOfEq1Result = 0.0
            for LeftBorder, RightBorder in zip(Part1ListOfIntegrationBorders[:-1], Part1ListOfIntegrationBorders[1:]):
                Part14thTermOfEq1Result += integrate.quad(IntegrandOfFirstSummandOf4thTerm, LeftBorder, RightBorder, args=(ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError), epsrel=RelativeError)[0]
        elif UseMultiprocessingIn4thTermOfEq1==True: # Multiprocessing is used, added in v1_4.
            Part1NumberOfSubprocesses=min(Part1NumberOfSubtasks,int(np.floor(multiprocessing.cpu_count()/(2*NumberOfIteratingProcesses)))) # multiprocessing.cpu_count()/2 gives the number of physical CPUs available on this machine. Division by NumberOfIteratingProcesses, gives the number of physical CPUs that are maximally available for iteration at one point of gamma.
            Part1PoolOfWorkers = multiprocessing.Pool(processes=Part1NumberOfSubprocesses) # Initialise multiple processes.
            if isinstance(NElectrons, types.FunctionType): # The case where NElectrons is a common picklable function.
                print('%-47s:         Determine 1. summand of 4. term for function-object with multiprocessing.' % ProcessID)
                Part14thTermOfEq1ResultQueue = [Part1PoolOfWorkers.apply_async(AuxiliaryIntegrateIntegrandOf4thTerm, args=[IntegrandOfFirstSummandOf4thTerm, LeftBorder, RightBorder, ga, n0, ValuesForxSyncSuperIt, ValuesFornSyncPs, Dotni, IterationCounter, LastValuesForgaIt, RelativeError], kwds={'NElectrons': NElectrons, 'ValuesForNextNElectrons': None}) for LeftBorder, RightBorder in zip(Part1ListOfIntegrationBorders[:-1], Part1ListOfIntegrationBorders[1:])]
            elif isinstance(NElectrons, interp1d): # Here, NElectrons has to be transferred to the child processes as arrays and interpolated in the child processes.
                print('%-47s:         Determine 1. summand of 4. term for interp1d-object with multiprocessing.' % ProcessID)
                Part14thTermOfEq1ResultQueue = [Part1PoolOfWorkers.apply_async(AuxiliaryIntegrateIntegrandOf4thTerm, args=[IntegrandOfFirstSummandOf4thTerm, LeftBorder, RightBorder, ga, n0, ValuesForxSyncSuperIt, ValuesFornSyncPs, Dotni, IterationCounter, LastValuesForgaIt, RelativeError], kwds={'NElectrons': None, 'ValuesForNextNElectrons': ValuesForNextNElectrons}) for LeftBorder, RightBorder in zip(Part1ListOfIntegrationBorders[:-1], Part1ListOfIntegrationBorders[1:])]
            Part14thTermOfEq1Result = np.sum([r.get() for r in Part14thTermOfEq1ResultQueue]) # Sum over the part-integrals
            Part1PoolOfWorkers.close() # Clear the Pool-object to exit the processes and to...
            Part1PoolOfWorkers.join() # ...retrieve the used RAM.
        print('%-47s:         Finished determining 1. summand of 4. term. Result: %s /(s*m^3)' % (ProcessID,Part14thTermOfEq1Result))
    # Integration of the second summand: 
    if HEPhotonsxga0<=Part1LowestIntegrationBorder: # If the HE-photon population is situated outside the integration-range, then the second summand vanishes. 
        return Part14thTermOfEq1Result
    else: # Consider the case, where the HE-photon population (at least partly) intersects with the integration-range. 
        Part2LowestIntegrationBorder=max(Part1LowestIntegrationBorder,HEPhotonsxga1)
        Part2BiggestIntegrationBorder=HEPhotonsxga0
        Part2ListOfIntegrationBorders = np.asarray([Part2LowestIntegrationBorder,Part2BiggestIntegrationBorder]) # Initialise the sampling of the integration borders.
        Part2PotentialAdditionalBorders = [] # Analogously to above, additional borders can be necessary:
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            Part2PotentialAdditionalBorders = np.append(Part2PotentialAdditionalBorders,[xgammaLimitClean(ValuesForxSyncSuperIt[-1],ga),1/ValuesForxSyncSuperIt[-1]]) # ValuesForxSyncSuperIt[-1] is xSync0Used.
        if Usedn0==n0MultiDelta:
            Part2PotentialAdditionalBorders = np.append(Part2PotentialAdditionalBorders,[xgammaLimitClean(i,ga) for i in x0MultiDelta])
            Part2PotentialAdditionalBorders = np.append(Part2PotentialAdditionalBorders,[1/i for i in x0MultiDelta])
        else:
            Part2PotentialAdditionalBorders = np.append(Part2PotentialAdditionalBorders,[xgammaLimitClean(x0,ga),1/x0])
        if InjectionType=='ResultsOfOldIt': # In this case, UsedDotni has kinks, which have to be used as integration borders.
            Part2PotentialAdditionalBorders = np.append(Part2PotentialAdditionalBorders,OldValuesForxgammanHEPhotons)
        for PotentialBorder in Part2PotentialAdditionalBorders: # Include these points into the list of integration-borders, if they are in the interior of the range.
            if (Part2LowestIntegrationBorder<PotentialBorder<Part2BiggestIntegrationBorder) and (PotentialBorder not in Part2ListOfIntegrationBorders):
                Part2ListOfIntegrationBorders=np.append(Part2ListOfIntegrationBorders,PotentialBorder)
        Part2ListOfIntegrationBorders = np.sort(Part2ListOfIntegrationBorders)
        # Now, subdivide each of the arising ranges and take the divisions as additional borders:
        Part2AdditionalBorders = []
        for LeftBorder, RightBorder in zip(Part2ListOfIntegrationBorders[:-1], Part2ListOfIntegrationBorders[1:]):
            Part2AdditionalBorders=np.append(Part2AdditionalBorders,UsedSampleIntBordersIn4thTerm(LeftBorder,RightBorder)[0][1:-1])
        for PotentialBorder in Part2AdditionalBorders: # Include these borders into the list of integration-borders, if they are not included yet.
            if PotentialBorder not in Part2ListOfIntegrationBorders:
                Part2ListOfIntegrationBorders=np.append(Part2ListOfIntegrationBorders,PotentialBorder)
        Part2ListOfIntegrationBorders = np.sort(Part2ListOfIntegrationBorders)
        print('%-47s:         Integration borders of 2. summand of 4. term:' % ProcessID, Part2ListOfIntegrationBorders)
        if not(Include1stIntRangeOf4thTerm) and len(Part2ListOfIntegrationBorders)>3: # The case where both the 1st integration range can be neglected and there are at least 3 ranges.
            Part2ListOfIntegrationBorders=Part2ListOfIntegrationBorders[1:] # Delete the first border, thus neglect the first range.
            print('%-47s:         Actually used integration borders:' % ProcessID, Part2ListOfIntegrationBorders)
        Part2NumberOfSubtasks = len(Part2ListOfIntegrationBorders)-1 # The integration along one subrange is considered as one subtask.
        if UseMultiprocessingIn4thTermOfEq1==False: # No multiprocessing.
            print('%-47s:         Determine 2. summand of 4. term without multiprocessing.' % ProcessID)
            Part24thTermOfEq1Result = 0.0
            for LeftBorder, RightBorder in zip(Part2ListOfIntegrationBorders[:-1], Part2ListOfIntegrationBorders[1:]):
                Part24thTermOfEq1Result += integrate.quad(IntegrandOfSecondSummandOf4thTerm, LeftBorder, RightBorder, args=(ga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError), epsrel=RelativeError)[0]
        elif UseMultiprocessingIn4thTermOfEq1==True: # Multiprocessing is used, added in v1_4.
            Part2NumberOfSubprocesses=min(Part2NumberOfSubtasks,int(np.floor(multiprocessing.cpu_count()/(2*NumberOfIteratingProcesses)))) # multiprocessing.cpu_count()/2 gives the number of physical CPUs available on this machine. Division by NumberOfIteratingProcesses, gives the number of physical CPUs that are maximally available for iteration at one point of gamma.
            Part2PoolOfWorkers = multiprocessing.Pool(processes=Part2NumberOfSubprocesses) # Initialise multiple processes.
            if isinstance(NElectrons, types.FunctionType): # The case where NElectrons is a common picklable function.
                print('%-47s:         Determine 2. summand of 4. term for function-object with multiprocessing.' % ProcessID)
                Part24thTermOfEq1ResultQueue = [Part2PoolOfWorkers.apply_async(AuxiliaryIntegrateIntegrandOf4thTerm, args=[IntegrandOfSecondSummandOf4thTerm, LeftBorder, RightBorder, ga, n0, ValuesForxSyncSuperIt, ValuesFornSyncPs, Dotni, IterationCounter, LastValuesForgaIt, RelativeError], kwds={'NElectrons': NElectrons, 'ValuesForNextNElectrons': None}) for LeftBorder, RightBorder in zip(Part2ListOfIntegrationBorders[:-1], Part2ListOfIntegrationBorders[1:])]
            elif isinstance(NElectrons, interp1d): # Here, NElectrons has to be transferred to the child processes as arrays and interpolated in the child processes.
                print('%-47s:         Determine 2. summand of 4. term for interp1d-object with multiprocessing.' % ProcessID)
                Part24thTermOfEq1ResultQueue = [Part2PoolOfWorkers.apply_async(AuxiliaryIntegrateIntegrandOf4thTerm, args=[IntegrandOfSecondSummandOf4thTerm, LeftBorder, RightBorder, ga, n0, ValuesForxSyncSuperIt, ValuesFornSyncPs, Dotni, IterationCounter, LastValuesForgaIt, RelativeError], kwds={'NElectrons':None, 'ValuesForNextNElectrons': ValuesForNextNElectrons}) for LeftBorder, RightBorder in zip(Part2ListOfIntegrationBorders[:-1], Part2ListOfIntegrationBorders[1:])]
            Part24thTermOfEq1Result = np.sum([r.get() for r in Part24thTermOfEq1ResultQueue]) # Sum over the part-integrals
            Part2PoolOfWorkers.close() # Clear the Pool-object to exit the processes and to...
            Part2PoolOfWorkers.join() # ...retrieve the used RAM.
        print('%-47s:         Finished determining 2. summand of 4. term. Result: %s /(s*m^3)' % (ProcessID,Part24thTermOfEq1Result))
        return Part14thTermOfEq1Result+Part24thTermOfEq1Result
    # Comment of version 9: IterationCounter was added to the function arguments.
    # Comment of version 10: LastValuesForgaIt was added to the function arguments. 

# Now, the fourth term is computed. According to my analysis this describes the spectral rate of pair-production events, that produce particles of energy ga via HE-photons. (n0 still is the population of external soft photons in both pair-production and IC-scattering. NElectrons is the population of particles. Dotni is the population of high-energy photons, that are injected artificially in the process of pair-production.) Thus it is a rate per unit volume and per unit final particle energy dga. It is in units of 1/(s*m^3). 
# In the case Dotni=DotniDelta and Dotni=DotniZero the integrand is divided into two summands for reasons of time-efficiency. In the other cases, the whole integrand is integrated.
if UsedDotni == DotniDelta or UsedDotni == DotniZero:
    def FourthTermOfEq1(ga,n0,ValuesForxSyncSuperIt,ValuesFornSyncPs,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,Include1stIntRangeOf4thTerm=True,ValuesForNextNElectrons=None,RelativeError=integrate.quad.__defaults__[3]):
        '''In case IncludeSynchrotron==False, and one has to set ValuesForxSyncSuperIt=None and ValuesFornSyncPs=None, whilst in the case IncludeSynchrotron==True, one has to plug in these arrays. nSyncPs is recreated as soon as necessary.'''
        ProcessID = multiprocessing.current_process().name+' Computer of 4. term (Dotni=0/Delta)'
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            PPThreshold=min(1.0/(4.0*x0),1.0/(4.0*ValuesForxSyncSuperIt[-1]))
        else: # This case has to be treated separately, because ValuesForxSyncSuperIt has the value None, which cannot be sliced.
            PPThreshold=1.0/(4.0*x0)
        if ga>PPThreshold:
            if IncludeSynchrotron and PerformanceMode!='CommonIteration':
                LowestIntegrationBorder=min(xgammaLimitClean(x0,ga),xgammaLimitClean(ValuesForxSyncSuperIt[-1],ga))
            else: # This case has to be treated separately, because ValuesForxSyncSuperIt has the value None, which cannot be sliced.
                LowestIntegrationBorder=xgammaLimitClean(x0,ga)
            BiggestIntegrationBorder=DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,LastValuesForgaIt[-1]) # The upper cut-off of the integral in brackets is the upper integration-border for the first part. The slicing of ValuesForxSyncSuperIt to ValuesForxSyncSuperIt[-1] is done internally.
            if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray: # If these two objects (arrays) are given, then nSyncPs is reconstructed. This case might be equal to the case PerformanceMode=='SuperIteration' and IncludeSynchrotron==True.
                #print('%-47s:         Doing interpolation of nSyncPs.' % ProcessID)
                nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
            else:
                print('%-47s:         Using nSyncPs=None.' % ProcessID)
                nSyncPs = None
            # At first, treat the first summand:
            if LowestIntegrationBorder>=BiggestIntegrationBorder:
                FirstSummandOf4thTermOfEq1Result = 0.0
            else:
                # The first summand is integrated.
                ListOfIntegrationBorders = np.asarray([LowestIntegrationBorder,BiggestIntegrationBorder]) # Initialise the sampling of the integration borders.
                PotentialAdditionalBorders = [] # Again, additional borders can be necessary:
                if IncludeSynchrotron and PerformanceMode!='CommonIteration':
                    PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[xgammaLimitClean(ValuesForxSyncSuperIt[-1],ga),1/ValuesForxSyncSuperIt[-1]]) # ValuesForxSyncSuperIt[-1] is xSync0Used.
                if Usedn0==n0MultiDelta:
                    PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[xgammaLimitClean(i,ga) for i in x0MultiDelta])
                    PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[1/i for i in x0MultiDelta])
                else:
                    PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[xgammaLimitClean(x0,ga),1/x0])
                for PotentialBorder in PotentialAdditionalBorders: # Include these points into the list of integration-borders, if they are in the interior of the range.
                    if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
                        ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
                ListOfIntegrationBorders = np.sort(ListOfIntegrationBorders)
                # Now, subdivide each of the arising ranges and take the divisions as additional borders:
                AdditionalBorders = []
                for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                    AdditionalBorders=np.append(AdditionalBorders,UsedSampleIntBordersIn4thTerm(LeftBorder,RightBorder)[0][1:-1])
                for PotentialBorder in AdditionalBorders: # Include these borders into the list of integration-borders, if they are not included yet.
                    if PotentialBorder not in ListOfIntegrationBorders:
                        ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)                
                ListOfIntegrationBorders = np.sort(ListOfIntegrationBorders)
                print('%-47s:         Integration borders of 1. summand of 4. term:' % ProcessID, ListOfIntegrationBorders)
                if not(Include1stIntRangeOf4thTerm) and len(ListOfIntegrationBorders)>3: # The case where both the 1st integration range can be neglected and there are at least 3 ranges.
                    ListOfIntegrationBorders=ListOfIntegrationBorders[1:] # Delete the first border, thus neglect the first range.
                    print('%-47s:         Actually used integration borders:' % ProcessID, ListOfIntegrationBorders)
                NumberOfSubtasks = len(ListOfIntegrationBorders)-1
                if UseMultiprocessingIn4thTermOfEq1==False: # No multiprocessing.
                    print('%-47s:         Determine 1. summand of 4. term without multiprocessing.' % ProcessID)
                    FirstSummandOf4thTermOfEq1Result = 0.0
                    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                        FirstSummandOf4thTermOfEq1Result += integrate.quad(IntegrandOfFirstSummandOf4thTerm, LeftBorder, RightBorder, args=(ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError), epsrel=RelativeError)[0] # I could show, that the bigger is Dotni, the more evaluation time is taken by this integration. However, the evaluation time of IntegrandOfFirstSummandOf4thTerm is not dependent on Dotni. The reason might be that the bigger is Dotni, the bigger is the variation within IntegrandOfFirstSummandOf4thTerm and bigger variations in an integrand usually make the integration more difficult.
                elif UseMultiprocessingIn4thTermOfEq1==True: # Multiprocessing is used, added in v1_4.
                    NumberOfSubprocesses=min(NumberOfSubtasks,int(np.floor(multiprocessing.cpu_count()/(2*NumberOfIteratingProcesses)))) # multiprocessing.cpu_count()/2 gives the number of physical CPUs available on this machine. Division by NumberOfIteratingProcesses, gives the number of physical CPUs that are maximally available for iteration at one point of gamma.
                    PoolOfWorkers = multiprocessing.Pool(processes=NumberOfSubprocesses) # Initialise multiple processes.
                    if isinstance(NElectrons, types.FunctionType): # The case where NElectrons is a common picklable function.
                        print('%-47s:         Determine 1. summand of 4. term for function-object with multiprocessing.' % ProcessID)
                        FirstSummandOf4thTermOfEq1ResultQueue = [PoolOfWorkers.apply_async(AuxiliaryIntegrateIntegrandOf4thTerm, args=[IntegrandOfFirstSummandOf4thTerm, LeftBorder, RightBorder, ga, n0, ValuesForxSyncSuperIt, ValuesFornSyncPs, Dotni, IterationCounter, LastValuesForgaIt, RelativeError], kwds={'NElectrons': NElectrons, 'ValuesForNextNElectrons': None}) for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:])]
                    elif isinstance(NElectrons, interp1d): # Here, NElectrons has to be transferred to the child processes as arrays and interpolated in the child processes.
                        print('%-47s:         Determine 1. summand of 4. term for interp1d-object with multiprocessing.' % ProcessID)
                        FirstSummandOf4thTermOfEq1ResultQueue = [PoolOfWorkers.apply_async(AuxiliaryIntegrateIntegrandOf4thTerm, args=[IntegrandOfFirstSummandOf4thTerm, LeftBorder, RightBorder, ga, n0, ValuesForxSyncSuperIt, ValuesFornSyncPs, Dotni, IterationCounter, LastValuesForgaIt, RelativeError], kwds={'NElectrons': None, 'ValuesForNextNElectrons': ValuesForNextNElectrons}) for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:])]
                    FirstSummandOf4thTermOfEq1Result = np.sum([r.get() for r in FirstSummandOf4thTermOfEq1ResultQueue]) # Sum over the part-integrals
                    PoolOfWorkers.close() # Clear the Pool-object to exit the processes and to...
                    PoolOfWorkers.join() # ...retrieve the used RAM.
                print('%-47s:         Finished determining 1. summand of 4. term. Result: %s /(s*m^3)' % (ProcessID,FirstSummandOf4thTermOfEq1Result))
            # Now, treat the second summand. The second summand is zero in case UsedDotni=DotniZero. In case UsedDotni=DotniDelta, it is determined via using the substitution feature of the Dirac-Delta-function.
            if UsedDotni == DotniZero:
                return FirstSummandOf4thTermOfEq1Result
            else:
                if HEPhotonsxga0 >= LowestIntegrationBorder:
                    if IncludeHEPhotonEscape==False:
                        SecondSummandOf4thTermOfEq1Result=pOfB18(ga,HEPhotonsxga0,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)*Dotni
                    elif IncludeHEPhotonEscape==True:
                        SecondSummandOf4thTermOfEq1Result=NormalisedSpectralPPProbability(ga,HEPhotonsxga0,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)*Dotni
                else:
                    SecondSummandOf4thTermOfEq1Result=0.0
                return FirstSummandOf4thTermOfEq1Result+SecondSummandOf4thTermOfEq1Result
        else: # The case ga<=min(1.0/(4.0*x0),1.0/(4.0*ValuesForxSyncSuperIt[-1])).
            return 0.0 # Zdz. showed that pairs are never created at energies below 1/(4*x0). Therefore this case is set to zero manually.
else:
    def FourthTermOfEq1(ga,n0,ValuesForxSyncSuperIt,ValuesFornSyncPs,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,Include1stIntRangeOf4thTerm=True,ValuesForNextNElectrons=None,RelativeError=integrate.quad.__defaults__[3]):
        '''In case IncludeSynchrotron==False, and one has to set ValuesForxSyncSuperIt=None and ValuesFornSyncPs=None, whilst in the case IncludeSynchrotron==True, one has to plug in these arrays. nSyncPs is recreated as soon as necessary.'''
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            PPThreshold=min(1.0/(4.0*x0),1.0/(4.0*ValuesForxSyncSuperIt[-1]))
        else: # This case has to be treated separately, because ValuesForxSyncSuperIt has the value None, which cannot be sliced.
            PPThreshold=1.0/(4.0*x0)
        if ga>PPThreshold:
            if UsedAlternativeFor4thTermOfEq1==1:
                return FourthTermOfEq1Alternative1(ga,n0,ValuesForxSyncSuperIt,ValuesFornSyncPs,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,Include1stIntRangeOf4thTerm=Include1stIntRangeOf4thTerm,ValuesForNextNElectrons=ValuesForNextNElectrons,RelativeError=RelativeError)
            elif UsedAlternativeFor4thTermOfEq1==2:
                return FourthTermOfEq1Alternative2(ga,n0,ValuesForxSyncSuperIt,ValuesFornSyncPs,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,Include1stIntRangeOf4thTerm=Include1stIntRangeOf4thTerm,ValuesForNextNElectrons=ValuesForNextNElectrons,RelativeError=RelativeError)
            else:
                print('Wrong value for UsedAlternativeFor4thTermOfEq1!')
        else: # The case ga<=min(1.0/(4.0*x0),1.0/(4.0*ValuesForxSyncSuperIt[-1])).
            return 0.0 # Zdz. showed that pairs are never created at energies below 1/(4*x0). Therefore this case is set to zero manually.
# Comment of version 9: IterationCounter was added to the function arguments. 
# Comment of version 10: LastValuesForgaIt was added to the function arguments.
# Comment of version 14: The if-test discriminating between > / <= 1/(4*x0) was introduced. By this, the function can be evaluated for ga smaller than 1/(4*x0), too, without raising errors.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    def EvaluateIntegrandOfFirstSummandOf4thTerm():
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            MinimumxgaLog=np.log10(min(xgammaLimitClean(x0,gamma),xgammaLimitClean(ValuesForxSyncSuperIt[-1],gamma))) # This is the lower integration border of 4th term.
        else: # This case has to be treated separately, because ValuesForxSyncSuperIt has the value None, which cannot be sliced.
            MinimumxgaLog=np.log10(xgammaLimitClean(x0,gamma))
        MaximumxgaLog = np.log10(DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,CurrentLastValuesForgaIt[-1]))
        FirstDivisionxgaLog = 1.0005*MinimumxgaLog
        SecondDivisionxgaLog = 0.9980*MaximumxgaLog
        ValuesForxgaLog = np.concatenate((np.logspace(MinimumxgaLog,FirstDivisionxgaLog,300,endpoint=False),np.logspace(FirstDivisionxgaLog,SecondDivisionxgaLog,800,endpoint=False),np.logspace(SecondDivisionxgaLog,MaximumxgaLog,900)))  # A sampling-range for xga.
        if not(IncludeSynchrotron):
            ValuesForIntegrandOfFirstSummandOf4thTerm = np.array([IntegrandOfFirstSummandOf4thTerm(i,gamma,Usedn0,None,None,CurrentNElectrons,1,CurrentLastValuesForgaIt,0.1) for i in ValuesForxgaLog])
        elif IncludeSynchrotron:
            ValuesForIntegrandOfFirstSummandOf4thTerm = np.array([IntegrandOfFirstSummandOf4thTerm(i,gamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons,1,CurrentLastValuesForgaIt,0.1) for i in ValuesForxgaLog])
        pl.figure(figsize=(12, 9), num="Cascade-equation: Integrand of 1. summand of 4. term versus incident photon energy")
        pl.loglog(ValuesForxgaLog,ValuesForIntegrandOfFirstSummandOf4thTerm)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        pl.xlabel("Incident photon energy $x_{\gamma}$", fontsize=16)
        pl.ylabel("Integrand of first summand of fourth term", fontsize=16)
        pl.savefig("Cascade-equation - Integrand of first summand of fourth term from %s.svg" % GetCurrentDate())
    # For n0=n0Planck and NElectrons=NElectronsGaussian the integrand of the first summand increases at first (due to the increase of pOfB18). After the peak, the integrand decreases (due to the decrease of IntegralOfBracketsOfEq1 which again is the effect of NElectrons).

    def EvaluateFourthTermOfEq1():
        MinimumGammaLog = np.log10(1/(4.0*x0))
        MaximumGammaLog = np.log10(1.0*DetermineUpperCutOffOf4thTermOfEq1(CurrentxSync0Used,CurrentNElectronsga0))
        DivisionPoint = 0.50*MaximumGammaLog
        ValuesForgaLog = np.concatenate((np.logspace(MinimumGammaLog,DivisionPoint,10,endpoint=False),np.logspace(DivisionPoint,MaximumGammaLog,100)))  # A sampling-range for ga. Logarithmic division of the interval is used.
        if not(IncludeSynchrotron):
            ValuesForFourthTermOfEq1 = np.array([FourthTermOfEq1(i,Usedn0,None,None,CurrentNElectrons,UsedDotni,1,CurrentLastValuesForgaIt,CurrentValuesForNElectrons,0.05) for i in ValuesForgaLog])
        elif IncludeSynchrotron:
            ValuesForFourthTermOfEq1 = np.array([FourthTermOfEq1(i,Usedn0,ValuesForxSyncSuperIt,ValuesFornSyncPs,CurrentNElectrons,UsedDotni,1,CurrentLastValuesForgaIt,CurrentValuesForNElectrons,0.05) for i in ValuesForgaLog])
        print('ValuesForgaLog:',ValuesForgaLog)
        print('ValuesForFourthTermOfEq1',ValuesForFourthTermOfEq1)
        pl.figure(figsize=(12, 9), num="Cascade-equation: Fourth term versus particle energy from %s" % GetCurrentDate())
        pl.loglog(ValuesForgaLog*x0,ValuesForFourthTermOfEq1)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        #pl.ylim(10**(1),10**(4))
        #pl.legend(loc="best", fontsize=14)
        pl.xlabel("$\gamma \cdot x_0$", fontsize=16)
        pl.ylabel("Fourth term", fontsize=16)
        #pl.xscale('log')
        pl.savefig("Cascade-equation - Fourth term versus particle energy from %s.svg" % GetCurrentDate())


## Considering all terms and iterating:

def NextNElectrons(ga,n0,ValuesForxSyncSuperIt,ValuesFornSyncPs,Dotni,DotNi,IterationCounter,LastValuesForgaIt,RelativeError,Include1stIntRangeOf4thTerm=True,NElectrons=None,ValuesForNextNElectrons=None):
    '''The electrons' kinetic equation (equivalent to equation 1 by Zdz.) is solved for the N(gamma) in the second term. Thus, N(gamma) = (the first term + the third term + the fourth term + the synchrotron-term)/electron-spectral-loss-rate. Then, this N(gamma) can be computed by evaluating all these terms. The first term, the third term, the fourth term and the synchrotron-term are all spectral injection-rates per unit volume and measured in 1/(s*m^3). The electron-spectral-loss-rate (denominator) is in 1/s and gives the spectral probability-rate for loss-processes (IC-down-scattering, escape and / or synchrotron). Hence, the computed N(gamma) is a spectral number-density and measured in 1/m^3, as it is supposed to be according to the above defined electron distribution.
IC-scattering and pair-production is always incorporated. Alternatively, one can incorporate electron escape and / or synchrotron-losses. If synchrotron-radiation is included, one has to properly insert ValuesForxSyncSuperIt and ValuesFornSyncPs, if it is not included, one has to set ValuesForxSyncSuperIt=None, ValuesFornSyncPs=None.
Actually, ValuesForNextNElectrons has kinks at the elements of LastValuesForgaIt. To execute the underlying integrations properly, one thus has to propagate LastValuesForgaIt as function-arguments. Principally, LastValuesForgaIt is different from ValuesForgaIt because the sampling of ga can change during the iteration.'''
    ProcessID = multiprocessing.current_process().name+' Computer'
    print('\n%-47s:     Computing value at gamma*x_0 = %s:' % (ProcessID,ga*x0))
    if NElectrons is not None and ValuesForNextNElectrons is None: # This case is tailored for NElectrons being a common function. To call this case, only specify NElectrons, but do not specify ValuesForNextNElectrons.
        print('%-47s:         No internal interpolation necessary.' % ProcessID)
    elif NElectrons is None and ValuesForNextNElectrons is not None: # This case is tailored for NElectrons being a non-picklable object, especially an interpolated object. In this case the interpolated object is reconstructed via LastValuesForgaIt and ValuesForNextNElectrons. To call this case, only specify ValuesForNextNElectrons, but do not specify NElectrons.
        #print('%-47s:         Doing internal interpolation.' % ProcessID)
        InternalInterpolatedObjectForValuesForNextNElectrons = interp1d(LastValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
        #print('%-47s:         Doing internal Assignment.' % ProcessID)
        NElectrons = InternalInterpolatedObjectForValuesForNextNElectrons
    else:
        raise WrongFunctionCallError(NextNElectrons.__name__)
    if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray: # If these two objects (arrays) are given, then nSyncPs is reconstructed. This case might be equal to the case PerformanceMode=='SuperIteration' and IncludeSynchrotron==True.
        #print('%-47s:         Doing interpolation of nSyncPs.' % ProcessID)
        nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
    else:
        print('%-47s:         Using nSyncPs=None.' % ProcessID)
        nSyncPs = None
    Term1=DotNi(ga)
    print('%-47s:         Term1 done:' % ProcessID, Term1, '/(s*m^3)')
    Term3=ThirdTermOfEq1(ga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError)
    print('%-47s:         Term3 done:' % ProcessID, Term3, '/(s*m^3)')
    Term4=FourthTermOfEq1(ga,n0,ValuesForxSyncSuperIt,ValuesFornSyncPs,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,Include1stIntRangeOf4thTerm=Include1stIntRangeOf4thTerm,ValuesForNextNElectrons=ValuesForNextNElectrons,RelativeError=RelativeError)
    print('%-47s:         Term4 done:' % ProcessID, Term4, '/(s*m^3)')
    if IncludeSynchrotron==False:
        Numerator = Term1+Term3+Term4
    elif IncludeSynchrotron==True: # The alternative with synchrotron-radiation.
        TermSynchrotron=SynchrotronSpectralProductionRate(ga,NElectrons,LastValuesForgaIt)
        print('%-47s:         TermSynchrotron done:' % ProcessID, TermSynchrotron, '/(s*m^3)')
        Numerator = Term1+Term3+Term4+TermSynchrotron
    Denominator = ElectronLossRate(ga,n0,nSyncPs,ValuesForxSyncSuperIt,LastValuesForgaIt,RelativeError) # ElectronLossRate autonomously chooses which processes are included. 
    print('%-47s:         Denominator:' % ProcessID, Denominator, '/s')
    Result = Numerator/Denominator
    print('%-47s:         Value of NElectrons:' % ProcessID, Result, '/m^3')
    return Result
    # In ThirdTermOfEq1 gaP was substituted for ga.
    # Comment of version 9: IterationCounter was added to the function arguments. 
    # Comment of version 10: LastValuesForgaIt was added to the function arguments.
    # Comment of version 12: The if-elif-else-test was added.
    # Comment of v1_0: COfA21 was replaced by ElectronDisappearanceRate in the denominator.

def AuxiliaryPushUpVersion1(ga,LowestPushedUpga,HighestPushedUpga):
    '''This function is equal to 1 left of LowestPushedUpga and right of HighestPushedUpga. In between, it is linearly decreasing with value 4 at LowestPushedUpga and 1 at HighestPushedUpga.'''
    return 1.0+3.0*Heavi(LowestPushedUpga,ga,HighestPushedUpga)*((ga-HighestPushedUpga)/(LowestPushedUpga-HighestPushedUpga))

def AuxiliaryPushUpVersion2(ga,LowestPushedUpga,HighestPushedUpga):
    '''This function is equal to 1 left of LowestPushedUpga and right of HighestPushedUpga. In between, it is decreasing via a power-law with value HighestPushUp at LowestPushedUpga and 1 at HighestPushedUpga.'''
    HighestPushUp = 6.0
    Coefficient = HighestPushUp**(1.0/(1.0-np.log10(LowestPushedUpga)/np.log10(HighestPushedUpga)))
    Index = np.log10(HighestPushUp)/(np.log10(HighestPushedUpga/LowestPushedUpga))
    return np.where(operator.or_(LowestPushedUpga > ga,ga > HighestPushedUpga), 1.0, 0.0) + Heavi(LowestPushedUpga,ga,HighestPushedUpga)*Coefficient*ga**(-Index)

def MoreIterationsNecessary(ga,NElectronsIterated,IterationCounter,ValuesForgaIt):
    '''This function determines whether more iterations have to be performed at ga (return True) or whether the iteration can be finished (return False).
    At first, and in any case, at least two iteration steps are computed.
    Then it is checked whether all the iterations at higher values of ValuesForgaIt than ga have been finished. (Notice that in principle LastValuesForgaIt mustn't be used here.)
    After that, it is checked whether the last iterated value is =0.
    
    If this is the case, at least one additional iteration step is computed with inclusion of the 1st range (and the boolean 'OneStepDoneAfterTaking1stIntRangeOf4thTermIntoAccount' is set equal to True).
    Then, it is checked whether more iteration steps have to be conducted or whether convergence at 0 is present.

    If the last iterated value is not =0, it is checked whether the relative change of the iterated values of NElectrons (the last and the second last element of StackOfIteratedValues) deceeds a certain threshold RelativeChangeToBeAchievedToTake1stIntRangeOf4thTermIntoAccount or whether the 1st range has already been taken into account ('Take1stIntRangeOf4thTermIntoAccount' is True).
    If this is the case, the 1st integration range of the 4th term is taken into account (and not neglected any more) from now on ('Take1stIntRangeOf4thTermIntoAccount' is set True). With inclusion of the 1st range, at least one additional iteration step is computed (and the boolean 'OneStepDoneAfterTaking1stIntRangeOf4thTermIntoAccount' is set equal to True).
    After this, it is checked for convergence of the values, i.e. it checks whether the relative change is below RelativeChangeToBeAchieved or whether more iterations have to be performed.
    
    If DefaultTake1stIntRangeOf4thTermIntoAccount is True, the 1st integration range is always and in every step taken into account.'''
    ProcessID = multiprocessing.current_process().name+' Checker'
    print('\n%-47s:     Checking convergence:' % ProcessID)
    if IterationCounter in [0,1]:
        print('%-47s:         More iterations necessary (due to iteration is just beginning).' % ProcessID)
        return True # At least two iterations have to be performed always.
    else: # The case IterationCounter > 1.
        ValuesBiggerThanga = ValuesForgaIt[ga<ValuesForgaIt] # This returns an array containing those elements of ValuesForgaIt that are bigger than ga.
        BooleansOfIterationFinished = np.asarray([NElectronsIterated[g]['IterationFinished'] for g in ValuesBiggerThanga]) # This array contains the booleans that state whether the iteration at the respective ga out of ValuesBiggerThanga is already finished or not.
        if False in BooleansOfIterationFinished: # This means that at least one iteration at higher ga has not been finished yet.
            print('%-47s:         More iterations necessary (due to higher iterations not finished yet).' % ProcessID)
            return True
        else: # The case "False not in BooleansOfIterationFinished" which means that all iterations at higher ga have been finished.
            SecondLastElementInStackOfIteratedValues = NElectronsIterated[ga]['StackOfIteratedValues'][-2]
            FirstLastElementInStackOfIteratedValues = NElectronsIterated[ga]['StackOfIteratedValues'][-1]
            if FirstLastElementInStackOfIteratedValues==0: # In this case, RelativeChange would encounter division by zero. So, RelativeChange cannot be determined.
                NElectronsIteratedAtga = NElectronsIterated[ga]
                NElectronsIteratedAtga['Take1stIntRangeOf4thTermIntoAccount']=True # Once this is set = True, the 1st range of the 4th term is computed at this ga.
                NElectronsIterated[ga] = NElectronsIteratedAtga
                if NElectronsIterated[ga]['OneStepDoneAfterTaking1stIntRangeOf4thTermIntoAccount']==False: # The next five lines ensure that at least one additional iteration step is performed after the previous criterion was satisfied (1st range is taken into account). This is necessary to be safe that the relative change is determined at least once for an iteration step after inclusion of the 1st range.
                    NElectronsIteratedAtga = NElectronsIterated[ga]
                    NElectronsIteratedAtga['OneStepDoneAfterTaking1stIntRangeOf4thTermIntoAccount']=True # Once this case is applied for this ga, it automatically ensures that it is applied never again.
                    NElectronsIterated[ga] = NElectronsIteratedAtga
                    print('%-47s:         First step with inclusion of 1st int. range of 4th term (flag was set).' % ProcessID)
                    return True
                else:
                    if SecondLastElementInStackOfIteratedValues!=0: # Nevertheless convergence is not achieved yet, so do one more iteration step.
                        print('%-47s:         More iterations necessary (due to change from !=0 to =0).' % ProcessID)
                        return True
                    else: # ...but (in the case SecondLastElementInStackOfIteratedValues==0) convergence is achieved, so...
                        print('%-47s:         Iteration finished (with convergence at =0).' % ProcessID)
                        NElectronsIteratedAtga = NElectronsIterated[ga]
                        NElectronsIteratedAtga['RelativeChangeInFinalIterationStep']='Convergence at =0' # Store the relative change in this special case.
                        NElectronsIterated[ga] = NElectronsIteratedAtga
                        return False # finish the iteration.
            else: # In the case FirstLastElementInStackOfIteratedValues!=0, RelativeChange can be determined...
                RelativeChange = (np.abs(SecondLastElementInStackOfIteratedValues-FirstLastElementInStackOfIteratedValues))/FirstLastElementInStackOfIteratedValues # Determine the relative change of the value of NElectrons in the last iteration.
                print('%-47s:         Relative change: %s' % (ProcessID,RelativeChange))
                if ga*x0>10: # The Klein-Nishina-regime.
                    RelativeChangeToBeAchievedToTake1stIntRangeOf4thTermIntoAccount = 0.1 # This is the relative change, that has to be deceeded for the first integration range of the 4th term being taken into account. This can be a critical value. Remark that this has to be bigger than RelativeChangeToBeAchieved, in any case.
                else:
                    RelativeChangeToBeAchievedToTake1stIntRangeOf4thTermIntoAccount = 0.1*(ga*x0/10)**0.5 # In the Thomson-regime it decreases.
                if RelativeChange > RelativeChangeToBeAchievedToTake1stIntRangeOf4thTermIntoAccount and NElectronsIterated[ga]['Take1stIntRangeOf4thTermIntoAccount']==False: # The second condition ensures that this test evaluates to False as soon as the 1st range has once been taken into account, as well as the case when DefaultTake1stIntRangeOf4thTermIntoAccount=True.
                    print('%-47s:         More iterations necessary without 1st int. range (due to required precision not met yet).' % ProcessID)
                    return True
                else:
                    NElectronsIteratedAtga = NElectronsIterated[ga]
                    NElectronsIteratedAtga['Take1stIntRangeOf4thTermIntoAccount']=True # Once this is set = True, the 1st range of the 4th term is computed at this ga.
                    NElectronsIterated[ga] = NElectronsIteratedAtga
                    if NElectronsIterated[ga]['OneStepDoneAfterTaking1stIntRangeOf4thTermIntoAccount']==False: # The next five lines ensure that at least one additional iteration step is performed after the previous criterion was satisfied (1st range is taken into account). This is necessary to be safe that the relative change is determined at least once for an iteration step after inclusion of the 1st range.
                        NElectronsIteratedAtga = NElectronsIterated[ga]
                        NElectronsIteratedAtga['OneStepDoneAfterTaking1stIntRangeOf4thTermIntoAccount']=True # Once this case is applied for this ga, it automatically ensures that it is applied never again.
                        NElectronsIterated[ga] = NElectronsIteratedAtga
                        print('%-47s:         First step with inclusion of 1st int. range of 4th term (flag was set).' % ProcessID)
                        return True
                    else: 
                        if ga*x0>10: # The Klein-Nishina-regime.
                            RelativeChangeToBeAchieved = 0.001 # This is the relative change of the last two iterated values of NElectrons, that has to be achieved for an iteration to be finished. This can be a critical value, especially in or near the Thomson-regime, because the energy transfer is small in this regime and the change of NElectrons from one iteration step to the next step can be small, although the iteration is still far from convergence. The deeper one wants to iterate into the Thomson-regime, the smaller this value has to be.
                        else:
                            RelativeChangeToBeAchieved = 0.001*(ga*x0/10)**2 # In the Thomson-regime it decreases.
                        if RelativeChange > RelativeChangeToBeAchieved:
                            print('%-47s:         More iterations necessary with 1st int. range (due to no convergence yet).' % ProcessID)
                            return True
                        else:
                            print('%-47s:         Iteration finished (with relative change %s < %s).' % (ProcessID,RelativeChange,RelativeChangeToBeAchieved))
                            NElectronsIteratedAtga = NElectronsIterated[ga]
                            NElectronsIteratedAtga['RelativeChangeInFinalIterationStep']=RelativeChange # Store the relative change.
                            NElectronsIterated[ga] = NElectronsIteratedAtga
                            return False # If this case applies, an iteration is finished.

def DetermineRelativeError(IterationCounter):
    if IterationCounter==0: # The case of the zeroth iteration.
        RelativeError = 0.01 # This is the relative error, that has to be achieved in the computations of NextNElectrons.
    else: # The case IterationCounter != 0:
        RelativeError = max(0.001,0.01/IterationCounter) # This reduces the RelativeError in the course of the iteration but never deceedes 0.01.
    return RelativeError

def IterateOnePointPointwise(ga,NElectronsIterated,ValuesForxSyncSuperIt,ValuesFornSyncPs,SuperItCounter,ValuesForgaIt):
    '''This performs the iteration of one point of ga, in the case 'PointsFromRToLIterationPointwise' with multiprocessing. Now, the discrimination between ValuesForgaIt and LastValuesForgaIt is abandoned, because within one common iteration, the sampling of ga does not change. Both are equal and are called ValuesForgaIt.'''
    ProcessID = multiprocessing.current_process().name+' Iterator'
    print('\n%-47s: Current value of gamma*x_0: %s' % (ProcessID,ga*x0))
    CurrentIterationCounter = 0 # Initialise the iteration at this value of ga.
    # The iteration is performed until the convergence criterion is satisfied:
    while MoreIterationsNecessary(ga,NElectronsIterated,CurrentIterationCounter,ValuesForgaIt):
        CurrentValuesForNElectrons = np.asarray([NElectronsIterated[ga]['StackOfIteratedValues'][-1] for ga in ValuesForgaIt]) # Determine the list of latest values of NElectrons by catching the last element of each StackOfIteratedValues. The i. element of CurrentValuesForNElectrons corresponds to the i. element of ValuesForgaIt.
        CurrentNElectrons = interp1d(ValuesForgaIt, CurrentValuesForNElectrons, kind='linear', bounds_error=False, fill_value=0.0)   # If fill_value=0.0 is included, the value of the interpolating object outside the interpolation range is 0.0.
        print('\n%-47s:     Current iteration counter: %s' % (ProcessID,CurrentIterationCounter))
        CurrentRelativeError = DetermineRelativeError(CurrentIterationCounter)
        print('%-47s:     Current relative error: %s' % (ProcessID,CurrentRelativeError))
        IteratedValueForNElectrons = NextNElectrons(ga,Usedn0,ValuesForxSyncSuperIt,ValuesFornSyncPs,UsedDotni,UsedDotNi,1,ValuesForgaIt,CurrentRelativeError,Include1stIntRangeOf4thTerm=NElectronsIterated[ga]['Take1stIntRangeOf4thTermIntoAccount'],ValuesForNextNElectrons=CurrentValuesForNElectrons) # Now perform the computation of the next iterated value of the electrons' spectral number-density. IterationCounter=1 is used as NElectrons is an interpolated object in every cases and never an analytical function. Actually, one should plug in LastValuesForgaIt instead of ValuesForgaIt. However, they are equal in our case, because the sampling of ga does not change during a common iteration.
        NElectronsIteratedAtga = NElectronsIterated[ga] # This is necessary for a proper modification of the interior dictionary, cf. the tests in MyMultiprocessingManagerWorking or MyMultiprocessingManagerDictOfDicts.
        NElectronsIteratedAtga['StackOfIteratedValues'] = np.append(NElectronsIteratedAtga['StackOfIteratedValues'],IteratedValueForNElectrons) # Save the computed value at the end of the stack.
        NElectronsIteratedAtga['IterationCounter'] = CurrentIterationCounter # Save the CurrentIterationCounter.
        NElectronsIterated[ga] = NElectronsIteratedAtga # Reassign the modified dictionary.
        CurrentIterationCounter += 1 # Increment the counter.
    NElectronsIteratedAtga = NElectronsIterated[ga] # This is necessary for a proper modification of the interior dictionary.
    NElectronsIteratedAtga['IterationFinished'] = True # This signalises that convergence at this gamma has been achieved.
    NElectronsIterated[ga] = NElectronsIteratedAtga # Reassign the modified dictionary.
    ExportDataOfIteration(ga, NElectronsIterated, "Temporary", SuperItCounter) # Save an intermediate belay point for the case of a sudden abort of the iteration.
   
def WorkerToIteratePointwise(InputQueue,NElectronsIterated,ValuesForxSyncSuperIt,ValuesFornSyncPs,SuperItCounter,ValuesForgaIt):
    """Draw input-values and crunch them as long as there are some in the input-queue."""
    ProcessID = multiprocessing.current_process().name+' Worker'
    while True: # Do this as long as there are points in InputQueue.
        try: # Try to retrieve a ga-value from the queue and to execute the iteration at this value.
            print('\n%-47s: Try to retrieve next Currentga from queue.' % ProcessID)
            Currentga = InputQueue.get(block=False) # Get the next value of ga from the queue.
            print('%-47s: Retrieved Currentga = %s.' % (ProcessID,Currentga))
            IterateOnePointPointwise(Currentga,NElectronsIterated,ValuesForxSyncSuperIt,ValuesFornSyncPs,SuperItCounter,ValuesForgaIt) # Perform a complete iteration at Currentga.
            print('\n%-47s: Finished iteration.\n' % ProcessID)
        except Empty: # If the queue has no more elements, finish the while-loop and quit the worker.
            break

def PrintActuallyAchievedAccuracy(DictionaryNElectronsIterated,ValuesForgaIt):
    '''This retrieves the values of the relative change that occurred from the second last to the first last iterated value and prints them as well as the logarithmic mean.'''
    ValuesForRelativeChangeInFinalIterationStep = np.asarray([DictionaryNElectronsIterated[ga]['RelativeChangeInFinalIterationStep'] for ga in ValuesForgaIt]) # Compile an array which contains the values of the actually achieved accuracy.
    CleanedValuesForRelativeChangeInFinalIterationStep = np.asarray(ValuesForRelativeChangeInFinalIterationStep[np.asarray([i!='Convergence at =0' and i!=0 for i in ValuesForRelativeChangeInFinalIterationStep])],dtype=float) # The same list, just without the items 'Convergence at =0' and without the items =0.
    MeanLogRelativeChangeInFinalIterationStep = np.sum(np.log10(CleanedValuesForRelativeChangeInFinalIterationStep))/len(CleanedValuesForRelativeChangeInFinalIterationStep) # The arithmetic mean of the logarithms of the actually achieved accuracies.
    print("\nList of actually achieved accuracies:", ValuesForRelativeChangeInFinalIterationStep)
    print("\nCorresponding mean logarithmic actually achieved accuracy:", MeanLogRelativeChangeInFinalIterationStep)

if "MainProcess" in multiprocessing.current_process().name and UsedIterationScheme == 'IterationStepByStepPointsFromLToR' and not(IncludeSynchrotron): # The first condition prevents evaluation in child-processes.
    def EvaluateIteration(CurrentIterationCounter, CurrentNElectrons, LastValuesForgaIt, NElectronsIterated, NElectronsga0, LastValuesForNextNElectrons=None):
        '''The iteration-loop (according to the scheme 'IterationStepByStepPointsFromLToR') with time-measurement. If a new iteration is executed, then the keyword argument LastValuesForNextNElectrons mustn't be given. If an imported iteration is resumed, then LastValuesForNextNElectrons has to be given ValuesForNextNElectrons.
        Comment of v2_3: The iteration scheme IterationStepByStepPointsFromLToR seems outdated, obsolete and inferior in comparison to PointsFromRToLIterationPointwise. Therefore, synchrotron-radiation is not implemented here as soft target photons.'''
        if type(LastValuesForNextNElectrons)!=type(None):
            ValuesForNextNElectrons = LastValuesForNextNElectrons # This is necessary in case, when an iteration, that was already executed, is resumed. In this case, the keyword argument LastValuesForNextNElectrons has to be given ValuesForNextNElectrons.
        StartingTime = time.time() # To measure the spent time of the iteration, determine the point of time before the start of the iteration.
        while CurrentIterationCounter<20:
            print('\n%-19s: Iteration number %s:' % (multiprocessing.current_process().name,CurrentIterationCounter))
            if CurrentIterationCounter<24:
                NumberOfSamplingPointsForTotalRange = 30+1*CurrentIterationCounter # Number of sampling points for the range from Lowestga to Highestga.
                NumberOfSamplingPointsForLowerRange = 3.0+np.round(CurrentIterationCounter/2.0) # Number of sampling points for the range slightly above Lowestga, where some additional points are inserted.
                NumberOfSamplingPointsForUpperRange = 5.0+CurrentIterationCounter # Number of sampling points for the range slightly below Highestga, where some additional points are inserted.
                NumberOfSamplingPointsForVeryUpperRange = 10 # Number of sampling points for the range very very slightly below Highestga, where some additional points are inserted.
                Lowestga = DetermineLowestga(x0)
                LowerDivisionPointgaLog = np.log10(10*Lowestga/5)
                IntermediateDivisionPointgaLog = (1.0*np.log10(ThomsonKNBoundaryExPs)+3*np.log10(NElectronsga0))/4
                UpperDivisionPointgaLog = (1.0*np.log10(ThomsonKNBoundaryExPs)+9*np.log10(NElectronsga0))/10
                SecondHighestgaLog = (1.0*np.log10(ThomsonKNBoundaryExPs)+999*np.log10(NElectronsga0))/1000
                ValuesForgaLogLowerRange = np.logspace(np.log10(Lowestga),LowerDivisionPointgaLog,NumberOfSamplingPointsForLowerRange,endpoint=False) # A part of the sampling-range for ga, at low values of ga.
                ValuesForgaLogUpperRange = np.logspace(IntermediateDivisionPointgaLog,SecondHighestgaLog,NumberOfSamplingPointsForUpperRange,endpoint=False) # A part of the sampling-range for ga, at intermediate / high values of ga.
                ValuesForgaLogVeryUpperRange = np.logspace(UpperDivisionPointgaLog,SecondHighestgaLog,NumberOfSamplingPointsForVeryUpperRange,endpoint=False) # A part of the sampling-range for ga, at very high values of ga.
                ValuesForgaIt = np.logspace(np.log10(Lowestga),SecondHighestgaLog,NumberOfSamplingPointsForTotalRange) # A part of the sampling-range for ga. This range is used for the evaluation and interpolation of NextNElectrons.
                for PotentialBorder in ValuesForgaLogLowerRange: # Test for each element of ValuesForgaLogLowerRange ...
                    if PotentialBorder not in ValuesForgaIt: # whether it is not yet in the list ValuesForgaIt...
                        ValuesForgaIt=np.append(ValuesForgaIt,PotentialBorder) # ...and if so, append it to the list.
                for PotentialBorder in ValuesForgaLogUpperRange: # Test for each element of ValuesForgaLogUpperRange ...
                    if PotentialBorder not in ValuesForgaIt: # whether it is not yet in the list ValuesForgaIt...
                        ValuesForgaIt=np.append(ValuesForgaIt,PotentialBorder) # ...and if so, append it to the list.
                for PotentialBorder in ValuesForgaLogVeryUpperRange: # This is a remedy for the sparse sampling problem.
                    if PotentialBorder not in ValuesForgaIt:
                        ValuesForgaIt=np.append(ValuesForgaIt,PotentialBorder)
                if 1.0/(2.0*x0) not in ValuesForgaIt:
                    ValuesForgaIt = np.append(ValuesForgaIt,1.0/(2.0*x0)) # At 1/(2*x0) FourthTermOfEq1 has a kink, so, include this point.
                if 1.01*Lowestga not in ValuesForgaIt:
                    ValuesForgaIt = np.append(ValuesForgaIt,1.01*Lowestga) # Include a point, which is situated very slightly above Lowestga.
                if NElectronsga0 not in ValuesForgaIt:
                    ValuesForgaIt = np.append(ValuesForgaIt,NElectronsga0) # Include Highestga.
                ValuesForgaIt = np.sort(ValuesForgaIt)
            else:
                ValuesForgaIt = LastValuesForgaIt
            if CurrentIterationCounter==0: # The case of the zeroth iteration.
                CurrentRelativeError = 0.01 # This is the relative error, that has to be achieved in the computations of NextNElectrons.
            else: # The case CurrentIterationCounter != 0:
                CurrentRelativeError = max(0.001,0.01/CurrentIterationCounter) # This reduces the RelativeError in the course of the iteration.
            print('%-19s:     Current relative error: %s' % (multiprocessing.current_process().name,CurrentRelativeError))
            # Now perform the computation of the next iterated electrons' spectral number-density. Actually, do the evaluation of NextNElectrons for ga=ValuesForgaIt, n0=Usedn0, NElectrons=CurrentNElectrons, IterationCounter=CurrentIterationCounter, ValuesForgaIt=LastValuesForgaIt, Dotni=UsedDotni and DotNi=UsedDotNi. As this evaluation is quite CPU-intensive, multiprocessing may be used.
            if NumberOfIteratingProcesses==1: # Do not use multiprocessing:
                if isinstance(CurrentNElectrons, types.FunctionType): # In this case, CurrentNElectrons is a usual callable function.
                    ValuesForNextNElectrons = np.asarray([NextNElectrons(i,Usedn0,None,None,UsedDotni,UsedDotNi,CurrentIterationCounter,LastValuesForgaIt,CurrentRelativeError,NElectrons=CurrentNElectrons) for i in ValuesForgaIt])
                elif isinstance(CurrentNElectrons, interp1d): # In this case, CurrentNElectrons is an interpolated object.
                    ValuesForNextNElectrons = np.asarray([NextNElectrons(i,Usedn0,None,None,UsedDotni,UsedDotNi,CurrentIterationCounter,LastValuesForgaIt,CurrentRelativeError,ValuesForNextNElectrons=ValuesForNextNElectrons) for i in ValuesForgaIt])
            else: # Use multiprocessing:
                PoolOfWorkers = multiprocessing.Pool(processes=NumberOfIteratingProcesses) # Initialise the multiple processes.
                if isinstance(CurrentNElectrons, types.FunctionType): # In this case, CurrentNElectrons is a usual callable function.
                    QueueOfResultsOfNextNElectrons = [PoolOfWorkers.apply_async(NextNElectrons, args=[i,Usedn0,None,None,UsedDotni,UsedDotNi,CurrentIterationCounter,LastValuesForgaIt,CurrentRelativeError], kwds={'NElectrons': CurrentNElectrons}) for i in ValuesForgaIt] # Distribute the tasks and start the work.
                elif isinstance(CurrentNElectrons, interp1d): # In this case, CurrentNElectrons is an interpolated object.
                    QueueOfResultsOfNextNElectrons = [PoolOfWorkers.apply_async(NextNElectrons, args=[i,Usedn0,None,None,UsedDotni,UsedDotNi,CurrentIterationCounter,LastValuesForgaIt,CurrentRelativeError], kwds={'ValuesForNextNElectrons': ValuesForNextNElectrons}) for i in ValuesForgaIt] # Distribute the tasks and start the work.
                ListOfResultsOfNextNElectrons = [r.get() for r in QueueOfResultsOfNextNElectrons] # Collect the results into a list.
                PoolOfWorkers.close() # Clear the Pool-object to exit the processes and to...
                PoolOfWorkers.join() # ...retrieve the used RAM.
                ValuesForNextNElectrons = np.asarray(ListOfResultsOfNextNElectrons)
            # Now, a speed-up of the iteration's convergence in the case of n0Planck:
            if Usedn0==n0Planck and CurrentIterationCounter in [5,6,7]:
                ValuesForNextNElectrons = ValuesForNextNElectrons*AuxiliaryPushUpVersion2(ValuesForgaIt,Lowestga,10*Lowestga) # This auxiliary function slightly raises the values of ValuesForNextNElectrons in the range between Lowestga and 50*Lowestga.
            # Now, the interpolation of the computed values: 
            InterpolatedObjectForValuesForNextNElectrons = interp1d(ValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0)   # If fill_value=0.0 is included, the value of the interpolating object outside the interpolation range is 0.0.
            # For each iteration-step, save the results in an own dictionary:
            NElectronsIterated['%s. Iteration' % CurrentIterationCounter] = {}
            NElectronsIterated['%s. Iteration' % CurrentIterationCounter]['Values for gamma'] = ValuesForgaIt # Now, save ValuesForgaIt in this dictionary.
            NElectronsIterated['%s. Iteration' % CurrentIterationCounter]['Interpolated object'] = InterpolatedObjectForValuesForNextNElectrons # Save the interpolated object in the same dictionary.
            NElectronsIterated['%s. Iteration' % CurrentIterationCounter]['Values for next NElectrons'] = ValuesForNextNElectrons # Save ValuesForNextNElectrons in the same dictionary.
            # Now, the iteration: 
            CurrentNElectrons = InterpolatedObjectForValuesForNextNElectrons
            LastValuesForgaIt = ValuesForgaIt
            CurrentIterationCounter +=1
            # Save an intermediate belay point for the case of a sudden abort of the iteration:
            ExportDataOfIteration(CurrentIterationCounter, NElectronsIterated, "Temporary", None)
        EndingTime = time.time() # Now, determine the point of time after the end of the iteration.
        TimeInterval = (EndingTime-StartingTime)/3600.0 # Determine the time in hours, that was spent by the iteration. 
        print("Time spent on the iteration:", TimeInterval, "hours")
        DeleteFile('Run %s, N(gamma) versus gamma (data, temporary' % RunIdentifier) # Delete the temporary belay point savings...
        ExportDataOfIteration(CurrentIterationCounter, NElectronsIterated, "DateAndTime", None) # ...and save all data again with a final filename.
        return CurrentIterationCounter, CurrentNElectrons, LastValuesForgaIt, NElectronsIterated, ValuesForgaIt, ValuesForNextNElectrons
elif "MainProcess" in multiprocessing.current_process().name and UsedIterationScheme == 'PointsFromRToLIterationPointwise': # The first condition prevents evaluation in child-processes.
    def EvaluateIteration(CurrentNElectrons, ValuesForgaIt, CurrentValuesForNElectrons, DictOfNElectrons, ValuesForxSyncSuperIt, ValuesFornSyncPs, SuperItCounter=None, ResumeSuperIt=False):
        '''The iterative walkthrough (according to the scheme 'PointsFromRToLIterationPointwise) with time-measurement. DictOfNElectrons can either be NElectronsIterated (in the case of a pure common iteration) or NElectronsSuperIt (in the case of a common iteration within a super-iteration).'''
        StartingTime = time.time() # To measure the spent time of the iteration, determine the point of time before the start of the iteration.
        LastValuesForgaIt = ValuesForgaIt # Principally, one has to discriminate between LastValuesForgaIt and ValuesForgaIt, because the sampling of ga can change during the iteration. However, actually, the sampling does change from one SuperIt-step to the next one, but not during common iteration-steps within one SuperIt. Thus, this assignment can be set.
        if SuperItCounter: # In this case, DictOfNElectrons is NElectronsSuperIt.
            if ResumeSuperIt==False:
                NElectronsIterated = CreateANewNElectronsIterated(CurrentNElectrons,ValuesForgaIt,SuperItCounter) # Create a new dictionary for this common iteration.
            elif ResumeSuperIt==True:
                NElectronsIterated = DictOfNElectrons[SuperItCounter] # Reimport the iteration that has been started earlier but not finished.
        else:
            NElectronsIterated=DictOfNElectrons
        ValuesForgaItToBeWalkedThrough = ValuesForgaIt[np.logical_not(np.asarray([NElectronsIterated[ga]['IterationFinished'] for ga in ValuesForgaIt]))] # This yields a list of those values of ValuesForgaIt, for which the iteration has not been finished yet and has to be performed.
        if NumberOfIteratingProcesses==1: # Do not use multiprocessing, here the function WorkerToIteratePointwise is not used and instead an explicit while-loop is used:
            for Currentga in np.flipud(ValuesForgaItToBeWalkedThrough):
                print('\n%-47s: Current value of gamma*x_0: %s' % (multiprocessing.current_process().name+' Iterator',Currentga*x0))
                CurrentIterationCounter = 0 # Initialise the iteration at this Currentga.
                # The iteration is performed until the convergence criterion is satisfied:
                while MoreIterationsNecessary(Currentga,NElectronsIterated,CurrentIterationCounter,ValuesForgaIt):
                    CurrentValuesForNElectrons = np.asarray([NElectronsIterated[ga]['StackOfIteratedValues'][-1] for ga in ValuesForgaIt]) # Determine the list of latest values of NElectrons by catching the last element of each StackOfIteratedValues. The i. element of CurrentValuesForNElectrons corresponds to the i. element of ValuesForgaIt.
                    CurrentNElectrons = interp1d(ValuesForgaIt, CurrentValuesForNElectrons, kind='linear', bounds_error=False, fill_value=0.0) # If fill_value=0.0 is included, the value of the interpolating object outside the interpolation range is 0.0.
                    print('\n%-47s:     Current iteration counter: %s' % (multiprocessing.current_process().name+' Iterator',CurrentIterationCounter))
                    CurrentRelativeError = DetermineRelativeError(CurrentIterationCounter)
                    print('%-47s:     Current relative error: %s' % (multiprocessing.current_process().name+' Iterator',CurrentRelativeError))
                    IteratedValueForNElectrons = NextNElectrons(Currentga,Usedn0,ValuesForxSyncSuperIt,ValuesFornSyncPs,UsedDotni,UsedDotNi,1,LastValuesForgaIt,CurrentRelativeError,Include1stIntRangeOf4thTerm=NElectronsIterated[Currentga]['Take1stIntRangeOf4thTermIntoAccount'],ValuesForNextNElectrons=CurrentValuesForNElectrons) # Now perform the computation of the next iterated value of the electrons' spectral number-density. IterationCounter=1 is used as NElectrons is an interpolated object in every cases and never an analytical function.
                    NElectronsIterated[Currentga]['StackOfIteratedValues'] = np.append(NElectronsIterated[Currentga]['StackOfIteratedValues'],IteratedValueForNElectrons) # Save the computed value at the end of the stack.
                    NElectronsIterated[Currentga]['IterationCounter'] = CurrentIterationCounter # Save the CurrentIterationCounter.
                    CurrentIterationCounter += 1 # Increment the counter.
                NElectronsIterated[Currentga]['IterationFinished'] = True # This signalises that convergence at this gamma has been achieved.
                ExportDataOfIteration(Currentga, NElectronsIterated, "Temporary", SuperItCounter) # Save an intermediate belay point for the case of a sudden abort of the iteration.
        else: # Use multiprocessing:
            MyManager=multiprocessing.Manager() # Initialise the manager.
            MyManagersNElectronsIterated = MyManager.dict(NElectronsIterated) # Initialise the manager's dictionary with the argument NElectronsIterated.
            #print('MyManagersNElectronsIterated after initialisation =',MyManagersNElectronsIterated)
            MyInputQueue = multiprocessing.Queue() # This is the queue which will contain the values of gamma, at which the iteration still has to be performed.
            for i in np.flipud(ValuesForgaItToBeWalkedThrough): # Fill the queue with the reversed list. 
                MyInputQueue.put(i)
            NumberOfActualSamplingPoints = MyInputQueue.qsize() # This is the number of gamma-values, at which the iteration still has to be performed, i.e. the number of elements in the queue.
            ListOfProcesses = [multiprocessing.Process(target=WorkerToIteratePointwise, args=(MyInputQueue,MyManagersNElectronsIterated,ValuesForxSyncSuperIt,ValuesFornSyncPs,SuperItCounter,ValuesForgaIt)) for i in range(NumberOfIteratingProcesses)] # Create the processes for the iteration. Use InputQueue as input of ga and MyManagersNElectronsIterated as output and WorkerToIteratePointwise as the targeted function. The arguments ValuesForxSyncSuperIt and ValuesFornSyncPs define the target-SyncPs. The processes work through the function WorkerToIteratePointwise for the queue of input values and store the result in MyManagersNElectronsIterated.
            for p in ListOfProcesses: # Start the processes.
                p.start()
            for p in ListOfProcesses: # Finish the processes.
                p.join()
            NElectronsIterated=MyManagersNElectronsIterated # Retransform the manager's dict to the original dict.
        DeleteFile('data, temporary after gamma =') # Delete the temporary belay point savings...
        if SuperItCounter:
            DictOfNElectrons[SuperItCounter] = NElectronsIterated.copy() # Save the results of this common iteration in the super-iteration-dictionary.
        else:
            DictOfNElectrons = NElectronsIterated
        ExportDataOfIteration('Walked through all elements of ValuesForgaIt', DictOfNElectrons, "DateAndTime", SuperItCounter) # ...and save all data again with a final filename.        
        PrintActuallyAchievedAccuracy(NElectronsIterated,ValuesForgaIt)
        print('\nList of number of necessary iteration steps:', [NElectronsIterated[ga]['IterationCounter'] for ga in ValuesForgaIt])
        CurrentValuesForNElectrons = np.asarray([NElectronsIterated[ga]['StackOfIteratedValues'][-1] for ga in ValuesForgaIt]) # Update the list of latest values of NElectrons.
        CurrentNElectrons = interp1d(ValuesForgaIt, CurrentValuesForNElectrons, kind='linear', bounds_error=False, fill_value=0.0) # Update the interpolation.
        EndingTime = time.time() # Now, determine the point of time after the end of the iteration.
        TimeInterval = (EndingTime-StartingTime)/3600.0 # Determine the time in hours, that was spent by the iteration. 
        print("\nTime spent on the iteration:", TimeInterval, "hours")
        return CurrentNElectrons, CurrentValuesForNElectrons, DictOfNElectrons

def ExportDataOfIteration(CurrentIterationCounterOrCurrentga, DictToSave, FilenameExtension, SuperItCounter):
    '''Remember again, that the use of super-iterations is only allowed together with 'PointsFromRToLIterationPointwise'.
    FilenameExtension has to be "DateAndTime" or "Temporary".
    In the case 'IterationStepByStepPointsFromLToR', CurrentIterationCounterOrCurrentga has to be given CurrentIterationCounter in a call of this function. Then both the value of CurrentIterationCounter and the complete dictionary DictToSave (containing both the list of evaluation-points in gamma and the computed data points and the interpolated object for each iteration) are written as a string into a .dat-file.
    In the case 'PointsFromRToLIterationPointwise', CurrentIterationCounterOrCurrentga has to be given Currentga in a call of this function. Then the complete dictionary DictToSave are written as a string into a .dat-file. DictToSave can either be NElectronsIterated or NElectronsSuperIt. The name of the created file depends on what is given to FilenameExtension.'''
    np.set_printoptions(threshold=np.inf, precision=20)
    if SuperItCounter==None:
        if FilenameExtension=="DateAndTime":
            Outputfile = open("Run %s, N(gamma) versus gamma from %s (data).dat" % (RunIdentifier,GetCurrentDate()), 'w')
        elif FilenameExtension=="Temporary":
            if UsedIterationScheme == 'IterationStepByStepPointsFromLToR':
                Outputfile = open("Run %s, N(gamma) versus gamma (data, temporary after iteration %s).dat" % (RunIdentifier,CurrentIterationCounterOrCurrentga-1), 'w')
            elif UsedIterationScheme == 'PointsFromRToLIterationPointwise':
                Outputfile = open("Run %s, N(gamma) versus gamma (data, temporary after gamma = %s).dat" % (RunIdentifier,CurrentIterationCounterOrCurrentga), 'w')
        if UsedIterationScheme == 'IterationStepByStepPointsFromLToR':
            Outputfile.write('CurrentIterationCounter=%s\n' % CurrentIterationCounterOrCurrentga)
    else:
        if FilenameExtension=="DateAndTime":
            Outputfile = open("Run %s, N(gamma) versus gamma from %s (SuperIt=%s, data).dat" % (RunIdentifier,GetCurrentDate(),SuperItCounter), 'w')
        elif FilenameExtension=="Temporary":
            Outputfile = open("Run %s, N(gamma) versus gamma (SuperIt=%s, data, temporary after gamma = %s).dat" % (RunIdentifier,SuperItCounter,CurrentIterationCounterOrCurrentga), 'w')
    Outputfile.write('%s' % DictToSave)
    Outputfile.close()

def EvaluatePreprocessOldIterationForNewOne(ValuesForgaIt,ValuesForNextNElectrons, NElectronsIterated):
    '''To reduce the number of necessary iterations, it is advisable to choose a course of NElectrons, that lies as near as possible at the final course of NElectrons. For a given set of Usedn0, UsedDotni and UsedDotNi, the final NElectrons might all look similar, up to a stretching/compression, which might be dependent on the input values. If, for a given set of Usedn0, UsedDotni and UsedDotNi, any final NElectrons is known and imported via EvaluateImportDataOfIteration, this NElectrons (To be precise ValuesForgaIt and ValuesForNextNElectrons) is now stretched/compressed. The resulting ValuesForgaIt and ValuesForNextNElectrons is returned and all initialisations for a new iteration are performed. Thereby, it is pretended as if the stretched/compressed NElectrons was the result of a 0. iteration. The parameters with the prefix Old... have to be specified from the input values of the iteration-run, whose NElectrons was imported. Attention: This function should be used cautiously and the values of LogOffsetOfHighestga and LowestStretchingFactor are dependent on the set of Usedn0, UsedDotni and UsedDotNi. HighestStretchingFactor was intuitively found and might change for another set of Usedn0, UsedDotni and UsedDotNi.
    Comment of v2_4: This definition was not maintained since the introduction of "PointsFromRToLIterationPointwise" as it seems outdated now.'''
    Oldx0 = 10**(-4) # The value of x0 of the iteration-run, that is supposed to be used for the stretching/compression. Has to be dimensionless.
    OldExPsEnergyDensity = 1936 # The value of the energy-density of the soft photon-field of the iteration-run, that is supposed to be used for the new run. Units have to be J/m^3.
    OldHEPhotonsxga0Delta = 1.0*10**15/(me*(c**2)/e) # The value of HEPhotonsxga0Delta of the iteration-run, that is supposed to be used for the new run. Has to be dimensionless.
    OldDotniDelta = 10**1 # The value of DotniDelta of the iteration-run, that is supposed to be used for the new run. Has to be in units of 1/(s*m^3).
    OldValuesForgaIt = ValuesForgaIt # Rename the old list.
    OldValuesForNextNElectrons = ValuesForNextNElectrons # Rename the old list.
    NElectronsIterated = {} # Initialise a new, empty dictionary.
    pl.figure(figsize=(12, 9), num="Cascade-equation: N(gamma) versus particle energy (stretched distribution)")    
    pl.loglog(OldValuesForgaIt,OldValuesForNextNElectrons, label='Old NElectrons') # Plot the old values, which are supposed to be stretched.
    # Now, stretching along ga is done:
    LogOffsetOfLowestga = np.log10(Oldx0/x0) # This value is the translation in logarithmic space of Lowestgaga. It is natural to use the ratio of the old x0 to the new x0, here.
    LogOffsetOfHighestga = np.log10(HEPhotonsxga0Delta/OldHEPhotonsxga0Delta) # This value is the translation of the highest value of ga in logarithmic space. It is advisable to choose the ratio of new injected energy to old injected energy. Tricky value!
    def StretchingValuesOfgaLog(gaLog): 
        '''This function does the stretching in logarithmic space.'''
        Slope = 1.0+(LogOffsetOfHighestga-LogOffsetOfLowestga)/(np.log10(OldValuesForgaIt[-1])-np.log10(OldValuesForgaIt[0]))
        yOffset = np.log10(OldValuesForgaIt[-1])+LogOffsetOfHighestga-Slope*np.log10(OldValuesForgaIt[-1])
        return Slope*gaLog+yOffset
    ValuesForgaIt = 10**(StretchingValuesOfgaLog(np.log10(OldValuesForgaIt))) # The logarithmic values are stretched, then converted to linear space and then assigned 
    LastValuesForgaIt = ValuesForgaIt
    # Now, stretching along NElectrons is done:    
    LowestStretchingFactor = (HEPhotonsxga0Delta*DotniDelta)/(OldHEPhotonsxga0Delta*OldDotniDelta)*np.sqrt(OldExPsEnergyDensity/ExPsEnergyDensity) # This is the value which is multiplied with the first value of ValuesForgaIt (that is with that value of NElectrons that corresponds to the lowest gamma). It might roughly depend on the ratio of new injected energy-density to old injected energy-density and on the square root of the ratio of old background energy-density to new background energy-density. Tricky value!
    HighestStretchingFactor = LowestStretchingFactor**(1.0/2) # This is the value which is multiplied with the last value of ValuesForgaIt (that is with that value of NElectrons that corresponds to the highest gamma). It was found by trying out. Tricky value!
    def DetermineLogOffsetForNElectrons(gaLog): 
        '''This function determines the offset in logarithmic space (= the stretching factor in linear space) for each value of OldValuesForNextNElectrons that corresponds to a value of gaLog.'''
        Slope = (np.log10(HighestStretchingFactor)-np.log10(LowestStretchingFactor))/(np.log10(ValuesForgaIt[-1])-np.log10(ValuesForgaIt[0]))
        yOffset = np.log10(LowestStretchingFactor)-np.log10(ValuesForgaIt[0])*Slope
        return Slope*gaLog+yOffset
    LogOffsetForNElectrons = DetermineLogOffsetForNElectrons(np.log10(ValuesForgaIt)) # Determine the offset-values 
    ValuesForNextNElectrons = OldValuesForNextNElectrons*10**LogOffsetForNElectrons # Stretch the old values with the determined factor.
    InterpolatedObjectForValuesForNextNElectrons = interp1d(ValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0) # An interpolation is done...
    CurrentNElectrons = InterpolatedObjectForValuesForNextNElectrons # ...and assigned.
    CurrentIterationCounter = 1 # Pretend that only the 0th iteration was performed.
    NElectronsIterated['%s. Iteration' % (CurrentIterationCounter-1)] = {} # Save...
    NElectronsIterated['%s. Iteration' % (CurrentIterationCounter-1)]['Values for gamma'] = ValuesForgaIt # ...the...
    NElectronsIterated['%s. Iteration' % (CurrentIterationCounter-1)]['Interpolated object'] = InterpolatedObjectForValuesForNextNElectrons # ...stretched...
    NElectronsIterated['%s. Iteration' % (CurrentIterationCounter-1)]['Values for next NElectrons'] = ValuesForNextNElectrons # ...values...
    pl.loglog(ValuesForgaIt,ValuesForNextNElectrons, label='Stretched NElectrons') # ...and plot them.
    pl.legend(loc="best", fontsize=12)
    pl.xlabel("Energy $\gamma$", fontsize=16)
    pl.ylabel("$N(\gamma)$", fontsize=16)
    return CurrentIterationCounter, CurrentNElectrons, LastValuesForgaIt, NElectronsIterated, ValuesForgaIt, ValuesForNextNElectrons

if "MainProcess" in multiprocessing.current_process().name and UsedIterationScheme == 'IterationStepByStepPointsFromLToR' and not(IncludeSynchrotron): # The first condition prevents evaluation in child-processes.
    def EvaluatePlotForIteration(CurrentIterationCounter, NElectronsIterated, CurrentNElectrons): # Now, plotting of the iterated curves. This plotting function does not support plotting of super-iteration steps.
        pl.figure(figsize=(14, 12), num="Cascade-equation: N(gamma) versus particle energy in the case n_0(x)=%s" % Usedn0.__name__)
        # In figure 1 of Zdziarski, gamma*N(gamma) is plottet versus gamma*x0:
        Lowestga = DetermineLowestga(x0)
        if StartingNElectrons != NElectronsZero:
            if StartingNElectrons==NElectronsRecommendedZdzExact or StartingNElectrons==NElectronsRecommendedZdzApprox or StartingNElectrons==NElectronsRecommendedAllProcesses:
                ValuesForgamma = np.logspace(np.log10(Injectedga1),np.log10(Injectedga0),1000)
            elif StartingNElectrons == NElectronsUnity or StartingNElectrons == NElectronsPremonition:
                ValuesForgamma = np.logspace(np.log10(Lowestga),np.log10(CurrentNElectronsga0),1000)
            pl.loglog(ValuesForgamma*x0,ValuesForgamma*StartingNElectrons(ValuesForgamma), label='Initialisation', color='#FFBBBB', linewidth=4, alpha=0.6)
        for IterationCounter in range(0,CurrentIterationCounter):
            pl.loglog(NElectronsIterated['%s. Iteration' % IterationCounter]['Values for gamma']*x0,NElectronsIterated['%s. Iteration' % IterationCounter]['Values for gamma']*NElectronsIterated['%s. Iteration' % IterationCounter]['Values for next NElectrons'], label='After %s. iteration' % IterationCounter)
        ValuesForgamma = np.logspace(np.log10(Lowestga),np.log10(CurrentNElectronsga0),1000)
        #pl.loglog(ValuesForgamma*x0,ValuesForgamma*CurrentNElectrons(ValuesForgamma), label='CurrentNElectrons')
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=26, direction='inout', width=3, length=10, pad=7)
        ax.yaxis.set_tick_params(labelsize=26, direction='inout', width=3, length=10)
        ax.xaxis.set_tick_params(which='minor', labelsize=20, direction='inout', width=2, length=6, pad=7)
        ax.yaxis.set_tick_params(which='minor', labelsize=20, direction='inout', width=2, length=6)
        pl.ylim(10.0**(6),2*10.0**(8))
        pl.xlim(0.01,CurrentNElectronsga0*x0*1.2)
        #pl.legend(loc="best", fontsize=12)
        # pl.xlabel("Lorentz-factor times dimensionless LEP-energy $\gamma \cdot x_0$", family='serif', fontsize=30)
        # pl.ylabel("Lorentz-factor times spectral\nnumber-density $\gamma \cdot N(\gamma)$ in $\mathrm{m}^{-3}$", family='serif', fontsize=32)
        pl.xlabel("$\gamma \, x_1$", family='serif', fontsize=40)
        pl.ylabel("$\gamma \, N(\gamma) / \mathrm{m}^{-3}$", fontsize=40)
        pl.subplots_adjust(top=0.97, bottom=0.12, left=0.13, right=0.97)
        # Save the figure:
        pl.savefig("Run %s, N(gamma) versus particle energy from %s (iterations, plot).svg" % (RunIdentifier,GetCurrentDate()))
        pl.savefig("Run %s, N(gamma) versus particle energy from %s (iterations, plot).pdf" % (RunIdentifier,GetCurrentDate()))
elif "MainProcess" in multiprocessing.current_process().name and UsedIterationScheme == 'PointsFromRToLIterationPointwise': # The first condition prevents evaluation in child-processes.
    def EvaluatePlotForIteration(ValuesForgaIt, NElectronsIteratedOrNElectronsSuperIt, CurrentNElectrons, SuperItCounter=False, xSync0Used=None): # Now, plotting of the iterated curves.
        # In figure 1 of Zdziarski, gamma*N(gamma) is plottet versus gamma*x0, this is done now. If a SuperIt is plotted, xSync0Used has to be given.
        # First, the plot is initialised:
        if SuperItCounter:
            Figure = pl.figure(figsize=(14, 12), num="Cascade-equation: gamma*N(gamma) versus particle energy gamma*x_0 (Super-iteration counter = %s)" % SuperItCounter)
            NElectronsIterated=NElectronsIteratedOrNElectronsSuperIt[SuperItCounter] # Makes the following lines less complicated.
        else:
            Figure = pl.figure(figsize=(13, 11), num="Cascade-equation: gamma*N(gamma) versus particle energy gamma*x_0 (Iterated)")
            NElectronsIterated=NElectronsIteratedOrNElectronsSuperIt # Makes the following lines less complicated.
        pl.rc('font', family='serif')
        LeftyAxis = Figure.add_subplot(111)
        LeftyAxis.xaxis.set_tick_params(labelsize=35, direction='inout', width=3, length=18, pad=-2)
        LeftyAxis.yaxis.set_tick_params(labelsize=35, direction='inout', width=3, length=18)
        LeftyAxis.xaxis.set_tick_params(which='minor', labelsize=22, direction='inout', width=2, length=10)
        LeftyAxis.yaxis.set_tick_params(which='minor', labelsize=22, direction='inout', width=2, length=10)
        LeftyAxis.set_xlim(ValuesForgaIt[0]*x0,ValuesForgaIt[-1]*x0*1.2)
        pl.xlabel("$\gamma \, x_{1}$", fontsize=45)
        pl.ylabel("$\gamma \, N(\gamma) / \mathrm{m}^{-3}$", fontsize=45)
        #pl.xlabel("Lorentz factor times dimensionless LEP energy $\gamma \cdot x_{0,1}$", fontsize=34)
        #pl.ylabel("Lorentz factor times spectral\nnumber density $\gamma \cdot N(\gamma)$ in $\mathrm{m}^{-3}$", fontsize=34)
        
        # Plotting of the converged courses of NElectrons of previous SuperIt-steps:
        if SuperItCounter>1: # For all SuperIt-steps except for the first one:
            ColourList=['green','cyan','magenta','yellow']
            for i in range(1,SuperItCounter): # For all SuperIt-steps except the latest.
                ValuesForgammaOfPreviousConvergence = np.sort(np.asarray([j for j in NElectronsIteratedOrNElectronsSuperIt[i].keys()])) # Extract the list of values of gamma.
                ValuesOfPreviousConvergence = np.asarray([NElectronsIteratedOrNElectronsSuperIt[i][ga]['StackOfIteratedValues'][-1] for ga in ValuesForgammaOfPreviousConvergence]) # And extract the corresponding converged values of NElectrons.
                LeftyAxis.loglog(ValuesForgammaOfPreviousConvergence*x0,ValuesForgammaOfPreviousConvergence*ValuesOfPreviousConvergence, label='Converged $N$ after SuperIt %s' % i, color=ColourList[(i-1)%4])
        
        # Plotting of the initialisation curves:
        if SuperItCounter>1:
            ValuesForInitialisation = np.asarray([NElectronsIterated[ga]['StackOfIteratedValues'][1] for ga in ValuesForgaIt])
            LeftyAxis.loglog(ValuesForgaIt*x0,ValuesForgaIt*ValuesForInitialisation, label='Initialisation $N_{\mathrm{init}}(\gamma)$', color='red', linewidth=3)
        else:
            if StartingNElectrons != NElectronsZero:
                if StartingNElectrons==NElectronsRecommendedZdzExact or StartingNElectrons==NElectronsRecommendedZdzApprox or StartingNElectrons==NElectronsRecommendedAllProcesses or StartingNElectrons==OldCurrentNElectrons: # These are the cases where StartingNElectrons is based on UsedDotNi.
                    ValuesForgaInit = np.logspace(np.log10(Injectedga1),np.log10(Injectedga0),1000)
                elif StartingNElectrons == NElectronsUnity or StartingNElectrons == NElectronsPremonition:
                    ValuesForgaInit = np.logspace(np.log10(ValuesForgaIt[0]),np.log10(ValuesForgaIt[-1]),1000)
                ValuesForInitialisation = [StartingNElectrons(ga) for ga in ValuesForgaInit]
                LeftyAxis.loglog(ValuesForgaInit*x0,ValuesForgaInit*ValuesForInitialisation, label='Initialisation $N_{\mathrm{init}}(\gamma)$', color='#FFBBBB', linewidth=4, alpha=0.6)
        
        # Plotting of the injected electrons:
        if UsedDotNi != DotNiZero:
            ValuesForgaInje = np.logspace(np.log10(Injectedga1),np.log10(Injectedga0),301,endpoint=False)
            ValuesForInjection = np.asarray([NElectronsRecommendedAllProcesses(ga) for ga in ValuesForgaInje])
            LeftyAxis.loglog(ValuesForgaInje*x0,ValuesForgaInje*ValuesForInjection, label='Spect. inject. rate / spect. loss rate', color='blue', linewidth=2)
        
        # Plotting of the iteration:
        ValuesForgaItToBeWalkedThrough = ValuesForgaIt[np.asarray([NElectronsIterated[ga]['IterationFinished'] for ga in ValuesForgaIt])] # This yields a list of those values of ValuesForgaIt, for which the iteration has been finished yet and can be plotted.
        for ga in ValuesForgaItToBeWalkedThrough:
            NumberOfPointsToDraw = len(NElectronsIterated[ga]['StackOfIteratedValues'])-1 # The first point is always =0.
            for i in range(1,NumberOfPointsToDraw+1):
                IteratedValueForNElectrons = NElectronsIterated[ga]['StackOfIteratedValues'][i]
                if IteratedValueForNElectrons!=0: # In the case =0, it cannot be drawn into a loglog-plot.
                    LeftyAxis.scatter(ga*x0,ga*IteratedValueForNElectrons, 3, color=ColourRedToBlack(i,NumberOfPointsToDraw+1)) # Draw the point.
        ValuesForgaItDensely = np.logspace(np.log10(ValuesForgaIt[0]),np.log10(ValuesForgaIt[-1]),2000)
        if SuperItCounter:
            LeftyAxis.loglog(ValuesForgaItDensely*x0,ValuesForgaItDensely*CurrentNElectrons(ValuesForgaItDensely), label='Converged $N(\gamma)$', color='black', linewidth=2)
        else:
            LeftyAxis.loglog(ValuesForgaItDensely*x0,ValuesForgaItDensely*CurrentNElectrons(ValuesForgaItDensely), label='Converged $N(\gamma)$', color='black', linewidth=2)
        
        # Finalisation of the plot:
        LeftyAxis.set_ylim(10.0**5,2*10.0**8)
        pl.legend(loc="lower left", fontsize=26)
        LeftyAxis.spines['left'].set_linewidth(3)
        LeftyAxis.spines['right'].set_linewidth(3)
        LeftyAxis.spines['top'].set_linewidth(3)
        LeftyAxis.spines['bottom'].set_linewidth(3)
        
        # The following is to create a different top x-axis:
        xAxis2=LeftyAxis.twiny()
        xAxis2.set_xscale('log')
        xAxis2.set_xlim(LeftyAxis.get_xlim())
        xAxis2MajorTicks=np.asarray([10**(-4),10**(-3),10**(-2),10**(-1),10**(0)]) # The tick-labels in TeV.
        xAxis2MinorTicks=np.asarray([2*10**(-4),3*10**(-4),4*10**(-4),5*10**(-4),6*10**(-4),7*10**(-4),8*10**(-4),9*10**(-4),2*10**(-3),3*10**(-3),4*10**(-3),5*10**(-3),6*10**(-3),7*10**(-3),8*10**(-3),9*10**(-3),2*10**(-2),3*10**(-2),4*10**(-2),5*10**(-2),6*10**(-2),7*10**(-2),8*10**(-2),9*10**(-2),2*10**(-1),3*10**(-1),4*10**(-1),5*10**(-1),6*10**(-1),7*10**(-1),8*10**(-1),9*10**(-1),2*10**(0),3*10**(0),4*10**(0),5*10**(0),6*10**(0)]) # The minor tick-labels in TeV.
        xAxis2.set_xticks(xAxis2MajorTicks*x0*e*10**12/(me*c**2)) # Conversion from values of the top axis to the positions on the bottom axis.
        xAxis2.set_xticks(xAxis2MinorTicks*x0*e*10**12/(me*c**2), minor=True) # Conversion from values of the top axis to the positions on the bottom axis.
        xAxis2.set_xticklabels(xAxis2MajorTicks)
        xAxis2.set_xticklabels([], minor=True)
        xAxis2.xaxis.set_tick_params(labelsize=35, direction='inout', width=3, length=18, pad=8)
        xAxis2.xaxis.set_tick_params(which='minor', labelsize=22, direction='inout', width=2, length=10, pad=8)
        xAxis2.set_xlabel("$\gamma \, m_{\mathrm{e}} c^2 / \, {\mathrm{TeV}}$", fontsize=45, labelpad=20)
        
        if SuperItCounter:
            # The following is to create a different top x-axis:
            xAxis3=LeftyAxis.twiny()
            xAxis3.set_xscale('log')
            xAxis3.set_xlim(LeftyAxis.get_xlim())
            xAxis3Ticks=np.asarray([10**(1),10**(2),10**(3),10**(4),10**(5)]) # The tick-labels.
            xAxis3.set_xticks(xAxis3Ticks*x0/xSync0Used) # Conversion from values of the top axis to the positions on the bottom axis.
            xAxis3.set_xticklabels(xAxis3Ticks)
            xAxis3.xaxis.set_tick_params(labelsize=28, direction='inout', width=3, length=10, pad=3)
            xAxis3.xaxis.set_tick_params(which='minor', labelsize=22, direction='inout', width=2, length=6, pad=3)
            xAxis3.set_xlabel("Lorentz factor times dim.less SyncP energy $\gamma \cdot x_{\mathrm{sync}\,0,\,\mathrm{used}}$", fontsize=34)
            xAxis3.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
            xAxis3.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
            xAxis3.spines['bottom'].set_position(('outward', 120))
        
        pl.subplots_adjust(top=0.85, bottom=0.15, left=0.16, right=0.96)
        if not(SuperItCounter):
            pl.savefig("Run %s, N(gamma) versus gamma from %s (plot).svg" % (RunIdentifier,GetCurrentDate()))
            #pl.savefig("Run %s, N(gamma) versus gamma from %s (plot).pdf" % (RunIdentifier,GetCurrentDate()))
        elif SuperItCounter:
            pl.savefig("Run %s, N(gamma) versus gamma from %s (SuperIt=%s, plot).svg" % (RunIdentifier,GetCurrentDate(),SuperItCounter))

        # # Now, gamma*N(gamma) is plottet versus the non-normalised particle energy gamma:
        # # First, the plot is initialised:
        # if SuperItCounter:
        #     pl.figure(figsize=(18, 14), num="Cascade-equation: gamma*N(gamma) versus particle energy gamma (Super-iteration counter = %s)" % SuperItCounter)
        #     pl.rc('font', family='serif')
        #     pl.title('Super-iteration counter = %s' % SuperItCounter, fontsize=34)
        # else:
        #     pl.figure(figsize=(18, 14), num="Cascade-equation: gamma*N(gamma) versus particle energy gamma (Iterated)")
        #     pl.rc('font', family='serif')
        # ax = pl.gca()
        # ax.xaxis.set_tick_params(labelsize=28, direction='inout', width=3, length=10, pad=7)
        # ax.yaxis.set_tick_params(labelsize=28, direction='inout', width=3, length=10)
        # ax.xaxis.set_tick_params(which='minor', labelsize=22, direction='inout', width=2, length=6, pad=7)
        # ax.yaxis.set_tick_params(which='minor', labelsize=22, direction='inout', width=2, length=6)
        # pl.xlim(ValuesForgaIt[0],ValuesForgaIt[-1]*1.2)
        # pl.xlabel("Lorentz factor $\gamma$", fontsize=34)
        # pl.ylabel("Lorentz factor times spectral\nnumber density $\gamma \cdot N(\gamma)$ in $\mathrm{m}^{-3}$", fontsize=34)
        # 
        # # Plotting of the converged courses of NElectrons of previous SuperIt-steps:
        # if SuperItCounter>1: # For all SuperIt-steps except for the first one:
        #     for i in range(1,SuperItCounter): # For all SuperIt-steps except the latest.
        #         pl.loglog(ValuesForgammaOfPreviousConvergence,ValuesForgammaOfPreviousConvergence*ValuesOfPreviousConvergence, label='Converged $N$ after SuperIt %s' % i, color=ColourList[(i-1)%4])

        # Plotting of the initialisation curves:
        # if SuperItCounter > 1:
        #     pl.loglog(ValuesForgaIt,ValuesForgaIt*ValuesForInitialisation, label='Initialisation $\gamma \cdot N_{\mathrm{init}}(\gamma)$', color='red')
        # else:
        #     if StartingNElectrons != NElectronsZero:
        #         pl.loglog(ValuesForgaInit,ValuesForgaInit*ValuesForInitialisation, label='Initialisation $\gamma \cdot N_{\mathrm{init}}(\gamma)$', color='red')
        # 
        # # Plotting of the injected electrons:
        # if UsedDotNi != NElectronsZero:
        #     pl.loglog(ValuesForgaInje,ValuesForgaInje*ValuesForInjection, label='Sp. injection rate $\dot N(\gamma)$ / Sp. loss rate', color='blue')
        # 
        # # Plotting of the iteration:
        # for ga in ValuesForgaItToBeWalkedThrough:
        #     NumberOfPointsToDraw = len(NElectronsIterated[ga]['StackOfIteratedValues'])-1 # The first point is always =0.
        #     for i in range(1,NumberOfPointsToDraw+1):
        #         IteratedValueForNElectrons = NElectronsIterated[ga]['StackOfIteratedValues'][i]
        #         if IteratedValueForNElectrons!=0: # In the case =0, it cannot be drawn into a loglog-plot.
        #             pl.scatter(ga,ga*IteratedValueForNElectrons, 3, color=ColourRedToBlack(i,NumberOfPointsToDraw+1)) # Draw the point.
        # pl.loglog(ValuesForgaItDensely,ValuesForgaItDensely*CurrentNElectrons(ValuesForgaItDensely), label='$\gamma \cdot N(\gamma)$ after convergence', color='black')
        # 
        # # Finalisation of the plot:
        # pl.ylim(10.0**3,10.0**7)
        # pl.legend(loc="best", fontsize=26)
        # pl.subplots_adjust(top=0.95, bottom=0.11, left=0.13, right=0.96)
        # if not(SuperItCounter):
        #     pl.savefig("Run %s, N(gamma) versus gamma from %s (plot, non-normalised).svg" % (RunIdentifier,GetCurrentDate()))
        # elif SuperItCounter:
        #     pl.savefig("Run %s, N(gamma) versus gamma from %s (SuperIt=%s, plot, non-normalised).svg" % (RunIdentifier,GetCurrentDate(),SuperItCounter))

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    def EvaluateExportInputValues():
        '''Save all actually used global input values in a file.'''
        Outputfile = open("Run %s, N(gamma) versus gamma from %s (input values).dat" % (RunIdentifier,GetCurrentDate()), 'w')
        Outputfile.write('alphaViscosity: %s\n' % alphaViscosity)
        Outputfile.write('betaPressure: %s\n' % betaPressure)
        Outputfile.write('ADAFrInner: %s\n' % ADAFrInner)
        Outputfile.write('ADAFrOuter: %s\n' % ADAFrOuter)
        Outputfile.write('ADAFTemperatureElectron: %s\n' % ADAFTemperatureElectron)
        Outputfile.write('M9: %s\n' % M9)
        Outputfile.write('Dotm: %s\n' % Dotm)
        Outputfile.write('UseApproxxM: %s\n' % UseApproxxM)
        Outputfile.write('etaff: %s\n' % etaff)
        Outputfile.write('Toy1Distance: %s\n' % Toy1Distance)
        Outputfile.write('B4: %s\n' % B4)
        Outputfile.write('hGap: %s\n' % hGap)
        Outputfile.write('rCurvature: %s\n' % rCurvature)
        Outputfile.write('Toy2tauPPOnADAF: %s\n' % Toy2tauPPOnADAF)
        Outputfile.write('Toy2TestParticleInjectionRate: %s\n' % Toy2TestParticleInjectionRate)
        Outputfile.write('UsedAlternativeFor4thTermOfEq1: %s\n' % UsedAlternativeFor4thTermOfEq1)
        Outputfile.write('UsedSampleIntBordersIn4thTerm: %s\n' % UsedSampleIntBordersIn4thTerm.__name__)
        Outputfile.write('UsedSampleIntBordersOfxSync: %s\n' % UsedSampleIntBordersOfxSync.__name__)
        Outputfile.write('IncludeHEPhotonEscape: %s\n' % IncludeHEPhotonEscape)
        Outputfile.write('IncludeElectronEscape: %s\n' % IncludeElectronEscape)
        Outputfile.write('IncludeSynchrotron: %s\n' % IncludeSynchrotron)
        Outputfile.write('PerformanceMode: %s\n' % PerformanceMode)
        Outputfile.write('InjectionType: %s\n' % InjectionType)
        Outputfile.write('FactorNumberDensityToInjectionRate: %s\n' % FactorNumberDensityToInjectionRate)
        Outputfile.write('ExPsOuterRadius: %s\n' % ExPsOuterRadius)
        Outputfile.write('ExPsVeryOuterRadius: %s\n' % ExPsVeryOuterRadius)
        Outputfile.write('InteractionVolume: %s\n' % InteractionVolume)
        Outputfile.write('InjectedalphaPL: %s\n' % InjectedalphaPL)
        Outputfile.write('InjectedalphaExp: %s\n' % InjectedalphaExp)
        Outputfile.write('InjectedTemperature: %s\n' % InjectedTemperature)
        Outputfile.write('Injectedga0: %s\n' % Injectedga0)
        Outputfile.write('Injectedga1: %s\n' % Injectedga1)
        Outputfile.write('InjectedGaussianLocation: %s\n' % InjectedGaussianLocation)
        Outputfile.write('InjectedGaussianWidth: %s\n' % InjectedGaussianWidth)
        Outputfile.write('UseInjectedGaussianCoefficientManualInput: %s\n' % UseInjectedGaussianCoefficientManualInput)
        Outputfile.write('InjectedGaussianNormManualInput: %s\n' % InjectedGaussianNormManualInput)
        Outputfile.write('InjectedGaussianCoefficientManualInput: %s\n' % InjectedGaussianCoefficientManualInput)
        Outputfile.write('InjectedExpCoefficient: %s\n' % InjectedExpCoefficient)
        Outputfile.write('InjectedPLCoefficient: %s\n' % InjectedPLCoefficient)
        Outputfile.write('InjectedHeaviCoefficient: %s\n' % InjectedHeaviCoefficient)
        Outputfile.write('InjectedRelMaxCoefficient: %s\n' % InjectedRelMaxCoefficient)
        Outputfile.write('UsedDotNi: %s\n' % UsedDotNi.__name__)
        Outputfile.write('HEPhotonsalphaPL: %s\n' % HEPhotonsalphaPL)
        Outputfile.write('HEPhotonsbetaPL: %s\n' % HEPhotonsbetaPL)
        Outputfile.write('HEPhotonsPLSharpnessIndex: %s\n' % HEPhotonsPLSharpnessIndex)
        Outputfile.write('HEPhotonsExpCutOff: %s\n' % HEPhotonsExpCutOff)
        Outputfile.write('HEPhotonsPlanckTemperature: %s\n' % HEPhotonsPlanckTemperature)
        Outputfile.write('HEPhotonsxga0NonDelta: %s\n' % HEPhotonsxga0NonDelta)
        Outputfile.write('HEPhotonsxga0Delta: %s\n' % HEPhotonsxga0Delta)
        Outputfile.write('HEPhotonsxga1: %s\n' % HEPhotonsxga1)
        Outputfile.write('HEPhotonsxgaBreak: %s\n' % HEPhotonsxgaBreak)
        Outputfile.write('HEPhotonsGaussianLocation: %s\n' % HEPhotonsGaussianLocation)
        Outputfile.write('HEPhotonsGaussianWidth: %s\n' % HEPhotonsGaussianWidth)
        Outputfile.write('UseHEPhotonsGaussianCoefficientManualInput: %s\n' % UseHEPhotonsGaussianCoefficientManualInput)
        Outputfile.write('HEPhotonsGaussianNormManualInput: %s\n' % HEPhotonsGaussianNormManualInput)
        Outputfile.write('HEPhotonsGaussianCoefficientManualInput: %s\n' % HEPhotonsGaussianCoefficientManualInput)
        Outputfile.write('HEPhotonsExpCoefficient: %s\n' % HEPhotonsExpCoefficient)
        Outputfile.write('HEPhotonsPLCoefficient: %s\n' % HEPhotonsPLCoefficient)
        Outputfile.write('DotniDelta: %s\n' % DotniDelta)
        Outputfile.write('UsedDotni: %s\n' % NameOfUsedDotni)
        Outputfile.write('ImportedLineDataFile: %s\n' % ImportedLineDataFile)
        Outputfile.write('alphaPL: %s\n' % alphaPL)
        Outputfile.write('alphaExp: %s\n' % alphaExp)
        Outputfile.write('PlanckTemperature: %s\n' % PlanckTemperature)
        Outputfile.write('x0NonDelta: %s\n' % x0NonDelta)
        Outputfile.write('x0Delta: %s\n' % x0Delta)
        Outputfile.write('x0MultiDelta: %s\n' % x0MultiDelta)
        Outputfile.write('x1: %s\n' % x1)
        Outputfile.write('ExpCoefficient: %s\n' % ExpCoefficient)
        Outputfile.write('PLCoefficient: %s\n' % PLCoefficient)
        Outputfile.write('UseDeltaCoefficientManualInput: %s\n' % UseDeltaCoefficientManualInput)
        Outputfile.write('DeltaCoefficient: %s\n' % DeltaCoefficient)
        Outputfile.write('CoefficientKLines: %s\n' % CoefficientKLines)
        Outputfile.write('Usedn0: %s\n' % Usedn0.__name__)
        Outputfile.write('InteractionRadius: %s\n' % InteractionRadius)
        Outputfile.write('JetHalfOpeningAngle: %s\n' % JetHalfOpeningAngle)
        Outputfile.write('DistanceSourceEarthMpc: %s\n' % DistanceSourceEarthMpc)
        Outputfile.write('gamma: %s\n' % gamma)
        Outputfile.write('xgamma: %s\n' % xgamma)
        Outputfile.write('LowestgaParameter1: %s\n' % LowestgaParameter1)
        Outputfile.write('LowestgaParameter2: %s\n' % LowestgaParameter2)
        if UsedIterationScheme == 'PointsFromRToLIterationPointwise': 
            Outputfile.write('Range5Parameter: %s\n' % Range5Parameter)
            Outputfile.write('Range6Parameter: %s\n' % Range6Parameter)
            Outputfile.write('Range7Parameter: %s\n' % Range7Parameter)
            Outputfile.write('NumberOfSamplingPointsForRange1: %s\n' % NumberOfSamplingPointsForRange1)
            Outputfile.write('NumberOfSamplingPointsForRange2: %s\n' % NumberOfSamplingPointsForRange2)
            Outputfile.write('NumberOfSamplingPointsForRange3: %s\n' % NumberOfSamplingPointsForRange3)
            Outputfile.write('NumberOfSamplingPointsForRange4: %s\n' % NumberOfSamplingPointsForRange4)
            Outputfile.write('NumberOfSamplingPointsForRange5: %s\n' % NumberOfSamplingPointsForRange5)
            Outputfile.write('NumberOfSamplingPointsForRange6: %s\n' % NumberOfSamplingPointsForRange6)
            Outputfile.write('NumberOfSamplingPointsForRange7: %s\n' % NumberOfSamplingPointsForRange7)
        Outputfile.write('StartingNElectrons: %s\n' % StartingNElectrons.__name__)
        Outputfile.write('ArtificialPushUpOfInitialisation: %s\n' % ArtificialPushUpOfInitialisation)
        Outputfile.write('DefaultTake1stIntRangeOf4thTermIntoAccount: %s\n' % DefaultTake1stIntRangeOf4thTermIntoAccount)
        Outputfile.write('FactornHEPsToFlux: %s\n' % FactornHEPsToFlux)
        Outputfile.close()

# Documentation of version 1: The evaluation of ThirdTermOfEq1 seems to work correctly for n0Delta. For the three other cases, the output is not reasonable. The evaluation seems to have problems if the integration range is too wide. Thus, in version 2, the integration range is again divided internally into parts.

# Documentation of version 2: Integration of the third term is now completed and seems reasonable. Third term was evaluated and considered. Possible forms of NElectrons were included. Minor changes were performed. The integral in the brackets of term 4 were performed. 

# Documentation of version 3: The order of structuring was reorganised. The fourth term was completed and evaluated.

# Documentation of version 4: The part "Injection of primary particles" was revised. The iteration was done. Minor changes in the whole file were done. Most test-blocks have been deleted. Minor fixes.

# Documentation of version 5: The findings of "1988ApJ...335..786Z - Determining an upper cut-off.png" were realised. 

# Documentation of version 6: The subsections "Injection of primary particles" and "Injection of primary HE-photons" were moved to part 1. UpperCutOffOf4thTermOfEq1, _ExportData, _ImportData and ImportFileToString were inserted. Minor changes.

# Documentation of version 7: Bug fix in FourthTermOfEq1.

# Documentation of version 8: Revision of NElectronsRecommended. Auxiliary functions were interchanged with part 1 and 3. The import-section was changed. _EvaluateIntegrandOfBracketsOfEq1ForCertainxga was added.

# Documentation of version 9: It was realised that global variables, that change their value during the program flow (like e. g. CurrentIterationCounter) and that appear inside an if-test inside of a function, have to be added to a function's arguments. Otherwise, the if-test is inevitably decided already at the definition of the function. For this reason, IterationCounter was added to the arguments of ThirdTermOfEq1, IntegralOfBracketsOfEq1, IntegrandOfFirstSummandOf4thTerm, CompleteIntegrandOf4thTerm, FourthTermOfEq1Alternative1, FourthTermOfEq1Alternative2, UsedAlternativeFor4thTermOfEq1, FourthTermOfEq1 and NextNElectrons.
# Additionally, in ThirdTermOfEq1 and in IntegralOfBracketsOfEq1 the if-test-conditions, that decide which borders are additionally added to the list of integration borders, have been altered.

# Documentation of version 10: Minor changes in _EvaluatePlotForIteration. ValuesForgaLogComputation was moved into the iteration-loop and the number of iteration points was made dependent on CurrentIterationCounter. Based on this, everything connected with the dictionary NElectronsIterated and with ValuesForgaLogComputation was adapted. Especially, ValuesForgaLogComputation was added to the arguments of those functions that use ValuesForgaLogComputation as their internal integration borders.

# Documentation of version 11: Multiprocessing was added in the iteration-loop. For all functions, beginning with underscore, the underscore was deleted and the function-name was prefixed with "Evaluate", if not done yet. This is necessary for making it possible to use the Evaluate-functions of part2 from part3. Additionally, the functions EvaluateIteration, EvaluateExportDataOfIteration, EvaluateImportDataOfIteration and EvaluatePlotForIteration have been adapted in such a way as to be usable via function-call and to make them able to access and alter the global variables. Minor changes.

# Documentation of version 12: The multiprocessing implementation didn't work in version 11. To be precise, the InterpolatedObject could not be pickled. This means that it could not be imported into the child-processes. Therefore the program had to be changed such that not the InterpolatedObject but only the values, that are interpolated, are imported into the child-processes. Then, the interpolation has to be done within the child-processes.

# Documentation of version 13: IntelligentStorage and ApplyIntelligentStorage have been added. Instead, the complicated and intricate procedure within SecondTermOfNumeratorOfEq8 to use AuxiliaryQuantityOfSecondTermOfNumeratorOfEq8 could be erased. The if-test in EvaluateExportInputValueshas been deleted as it is not needed any more.

# Documentation of version 14: I realised that the fourth term FourthTermOfEq1 cannot be evaluated for ga below 1/4x0. However, the third term can be evaluated there and is non-vanishing above 1/(8*x0). The first term as well as COfA21 can be evaluated below 1/4x0, too. Hence, to make FourthTermOfEq1 evaluable below 1/4x0, too, the if-test discriminating between > / <= 1/(4*x0) was introduced in FourthTermOfEq1. By this, NextNElectrons can and should be evaluated below 1/4x0, too. Furthermore, Lowestga was lowered to 1/(8*x0) in the iteration.

# Documentation of version 15: Up to now, there occurred differences in the results between python2- and python3-evaluations. These were due to python2 treating 1/2 as integer-division. To find remedy / make the code also valid for python2, all numbers that are definitively no integer are rendered with a decimal point. Furthermore, quantities that can be simplified, like e. g. 1/2, are simplified, e. g. to 0.5.
# In the case of n0Planck, convergence is quite slow below ga=1/(2*x0). Therefore a subroutine is introduced, which pushes NElectrons up in the range Lowestga<ga<10*Lowestga during the first couple of iteration steps. This subroutine is called AuxiliaryPushUpVersion1 / 2.
# EvaluateExportDataOfIteration was modified and renamed to ExportDataOfIteration and included into the iteration-loop.
# EvaluatePlotForIteration has been slightly modified.
# ValuesForgaLogLowerRange was added to ValuesForgaLogComputation.
# An additional integration border was introduced in IntegralOfBracketsOfEq1.
# EvaluatePreprocessOldIterationForNewOne was created.
# The epsrel-keyword-argument has been introduced in some integrals.
# Minor changes.

# Documentation from v1_0 on will proceed in part 3...