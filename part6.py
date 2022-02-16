## ------------------------------------------------------------------------------------     o   o
## Electromagnetic cascades                                                                   I
##                                                                                          \.../
## Numerical solution of the kinetic equation for synchrotron-photons
## ------------------------------------------------------------------------------------

import os
import multiprocessing

from part5 import * # Import of the file.


## Considering the synchrotron kinetic equation

# Consider synchrotron-photons (SyncPs). The kinetic equation of the synchrotron-photons (KESP) is yielded by balancing the production-rate with the loss-rate, cf. "Idea of iteration.png".

def SyncPsSpectralProductionRateIntegrand(ga,xSync,NElectrons,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the integrand of the integral in the numerator. It is NElectrons*SynchrotronEmissivity/xSync. SynchrotronEmissivity(ga)/xSync gives the number of synchrotron-photons of energy xSync, that are emitted by one electron of energy ga per unit time and per unit dimensionless photon-energy-interval. Then, NElectrons*SynchrotronEmissivity/xSync gives the number of synchrotron-photons of energy xSync, that are emitted by all electrons in the energy-interval around ga per unit space-volume, per unit time, per unit dimensionless electron-energy-interval and per unit dimensionless photon-energy-interval. It is in units of 1/(m^3*s).
    In the following integral it is to be integrated over ga, so this has to be the first argument.'''
    return NElectrons(ga)*SynchrotronEmissivity(ga,xSync,RelativeError)/xSync

def SyncPsSpectralProductionRate(xSync,NElectrons,ValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is in units of 1/(s*m^3). It describes the spectral density-rate of SyncP-production-events from the electron population of spectral number-density NElectrons, in other words the number of SyncPs, that are produced per unit space-volume, per unit time-interval and per unit dimensionless photon-energy. (This quantity is called emissivity by Rieger.)
    NElectrons can either be a defined function multiplied with a Heaviside-function (with jumps at Injectedga1 and Injectedga0) or an interpolated object with kinks at all elements of ValuesForgaIt. Integration should be performed along the entire non-vanishing range of NElectrons.'''
    if isinstance(NElectrons, interp1d):
        ListOfIntegrationBorders=ValuesForgaIt # In this case, NElectrons is an interpolated object, which has kinks at each item of the ValuesForgaIt. So, take every kink as an element of ListOfIntegrationBorders.
    elif isinstance(NElectrons, types.FunctionType):
        LowestIntegrationBorder=Injectedga1
        BiggestIntegrationBorder=Injectedga0
        ListOfIntegrationBorders = SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
    #print('    Integration borders of SyncPsSpectralProductionRate: ', ListOfIntegrationBorders)
    SyncPsSpectralProductionRateResult = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        SyncPsSpectralProductionRateResult += integrate.quad(SyncPsSpectralProductionRateIntegrand, LeftBorder, RightBorder, args=(xSync,NElectrons,RelativeError), epsrel=RelativeError)[0]
    return SyncPsSpectralProductionRateResult

def SyncPsSpectralICScatteringRate(xga,xSync,ga):
    '''This is a measure for the spectral rate of IC-events, that have an incident electron with energy ga and an original photon with energy xSync and result in a final photon of energy xga, per unit time interval and per unit final energy xga. It is in units of m^3/s. It is the same as the PartOfIntegrandOfA1 behind n_0 in Zdz's eq. A1, just with the use of gaP=ga-xga. This quantity mustn't be evaluated at xga=0. But as this is not physically interesting, anyway, it will not be excluded. It does behave well for all reasonable xga. However, for big ga, it has a very sharp peak at values of xga, that are very slightly below xmaxOfA20(ga,xSync)'''
    EAstOverEAlt = EAstOfA3Alt(ga,xga)/EOfA2(ga,xSync) # Abbreviation.
    return IntegrandOfA1Coefficient/(EOfA2(ga,xSync)*ga) * ( rOfA4(ga,(ga-xga)) + (2.0-rOfA4(ga,(ga-xga)))*EAstOverEAlt - 2.0*EAstOverEAlt**2.0 - 2.0*EAstOverEAlt*np.log(1.0/EAstOverEAlt) )

def SyncPsICScatteringRate(xSync,ga,RelativeError=integrate.quad.__defaults__[3]):
    '''This is a measure for the rate of IC-events, that have an incident electron with energy ga and an original photon with energy xSync and result in any final photon of kinematically allowed energy per unit time interval. It is in units of m^3/s. '''
    BiggestIntegrationBorder = xmaxOfA20(ga,xSync)
    LowestIntegrationBorder = 0.001*BiggestIntegrationBorder # Physically, this should be =0. This might be a delicate value. By using EvaluateSyncPsICScatteringRate, one can see that very tiny deviations in the result of SyncPsICScatteringRate appear at small values of ga, if it is e.g. LowestIntegrationBorder=0.01.
    ListOfIntegrationBorders = SampleIntegrationBorders8(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
    #print('    Integration borders of SyncPsICScatteringRate: ', ListOfIntegrationBorders)
    SyncPsICScatteringRateResult = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        SyncPsICScatteringRateResult += integrate.quad(SyncPsSpectralICScatteringRate, LeftBorder, RightBorder, args=(xSync,ga), epsrel=RelativeError)[0]
    return SyncPsICScatteringRateResult

def SyncPsSpectralICScatteringRateAlt(gaP,xSync,ga):
    '''This is a measure for the spectral rate of IC-events, that have an incident electron with energy ga and an original photon with energy xSync and result in a final electron of energy gaP, per unit time interval and per unit final energy gaP. It is in units of m^3/s. It is the same as SyncPsSpectralICScatteringRate, just with the use of gaP=ga-xga.'''
    if gaP < gaPminOf6(xSync,ga):
        return 0.0
    elif gaP < ga:
        return PartOfIntegrandOfA1(xSync,ga,gaP)
    elif gaP == ga:
        return PartOfIntegrandOfA1InFrontOfBrackets(xSync,ga)
    else:
        return 0.0

def SyncPsICScatteringRateAlt(xSync,ga,RelativeError=integrate.quad.__defaults__[3]):
    '''This is a measure for the rate of IC-events, that have an incident electron with energy ga and an original photon with energy xSync and result in any final electron of kinematically allowed energy per unit time interval. It is in units of m^3/s. '''
    BiggestIntegrationBorder = ga
    LowestIntegrationBorder = gaPminOf6(xSync,ga)
    ListOfIntegrationBorders = SampleIntegrationBorders7(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
    #print('    Integration borders of SyncPsICScatteringRate: ', ListOfIntegrationBorders)
    SyncPsICScatteringRateResultAlt = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        SyncPsICScatteringRateResultAlt += integrate.quad(SyncPsSpectralICScatteringRateAlt, LeftBorder, RightBorder, args=(xSync,ga), epsrel=RelativeError)[0]
    return SyncPsICScatteringRateResultAlt

def SyncPsTotalICScatteringRateIntegrand(ga,xSync,NElectrons,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the integrand of the outer integral in the denominator. It is NElectrons*SyncPsICScatteringRateAlt. It gives the probability of scattering events with original photon energy xSync and original electron energy ga per unit time interval and per original electron energy interval. It is in units of 1/s.
    In the following integral it is to be integrated over ga, so this has to be the first argument. SyncPsICScatteringRate could also be used, however SyncPsICScatteringRateAlt is slightly faster.'''
    return NElectrons(ga)*SyncPsICScatteringRateAlt(xSync,ga,RelativeError)

def SyncPsTotalICScatteringRate(xSync,NElectrons,ValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is in units of 1/s. It describes the probability-rate of IC-upscattering-events of SyncPs of energy xSync off electrons of spectral number-density NElectrons, in other words the probability that a SyncP is upscattered per unit time-interval. (This quantity might be called spectral IC up-scattering rate, in analogy to the 2019 A&A paper.)
    NElectrons can either be a defined function multiplied with a Heaviside-function (with jumps at Injectedga1 and Injectedga0) or an interpolated object with kinks at all elements of ValuesForgaIt. Integration should be performed along the entire non-vanishing range of NElectrons.'''
    if isinstance(NElectrons, interp1d):
        ListOfIntegrationBorders=ValuesForgaIt # In this case, NElectrons is an interpolated object, which has kinks at each item of the ValuesForgaIt. So, take every kink as an element of ListOfIntegrationBorders.
    elif isinstance(NElectrons, types.FunctionType):
        LowestIntegrationBorder=Injectedga1
        BiggestIntegrationBorder=Injectedga0
        ListOfIntegrationBorders = SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
    #print('    Integration borders of SyncPsTotalICScatteringRate: ', ListOfIntegrationBorders)
    SyncPsTotalICScatteringRateResult = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        SyncPsTotalICScatteringRateResult += integrate.quad(SyncPsTotalICScatteringRateIntegrand, LeftBorder, RightBorder, args=(xSync,NElectrons,RelativeError), epsrel=RelativeError)[0]
    return SyncPsTotalICScatteringRateResult

def SyncPsLossRate(xSync,NElectrons,ValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is in units of 1/s and gives the spectral probability-rate of events that effectively make SyncPs disappear from the energy xSync. It includes IC-up-scattering as well as escape from the interaction region.'''
    return SynchrotronPhotonEscapeRate(xSync)+SyncPsTotalICScatteringRate(xSync,NElectrons,ValuesForgaIt,RelativeError)

def SyncPsOpticalDepth(Process,*PositionalArgumentsOfProcess):
    '''This is the optical depth of the synchrotron-photons with respect to a loss-process Process. It is dimensionless.
    It was yielded based on the relation loss-probability-rate-of-process = c*dtau/dlength, where tau is the optical depth and length the path-length. Solving for dtau and integrating yields tau = loss-probability-rate-of-process*MeanEscapeLength/c, where MeanEscapeLength is the mean distance a SyncP has to pass until it leaves the interaction region. Using c = MeanEscapeLength/MeanEscapeTime gives the following expression.
    Process can be any contribution to SyncPsLossRate or SyncPsLossRate itself, hence IC-scattering off the electrons (Process=SyncPsTotalICScatteringRate, arguments have to be xSync, NElectrons and ValuesForgaIt), escape from the region (Process=SynchrotronPhotonEscapeRate, arguments have to be xSync) or all these (Process=SyncPsLossRate, arguments have to be xSync, NElectrons and ValuesForgaIt).'''
    return Process(*PositionalArgumentsOfProcess)*MeanEscapeTime # This essentially is the equation tau = loss-probability-rate-of-process * MeanEscapeTime.

def SyncPsSpectralNumberDensity(xSync,NElectrons,ValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
    '''The spectral number-density of the SyncPs in units of 1/m^3.'''
    return SyncPsSpectralProductionRate(xSync,NElectrons,ValuesForgaIt,RelativeError)/SyncPsLossRate(xSync,NElectrons,ValuesForgaIt,RelativeError)

def SyncPsNumberDensity(MultiplicatorOfLowestBorder,MultiplicatorOfBiggestBorder,ValuesForgaIt,ValuesForNextNElectrons,RelativeError=integrate.quad.__defaults__[3]):
    '''The total number-density of the SyncPs in units of 1/m^3. Actually, one should integrate from 0 to infinity. However, the critical frequencies of the upper and lower cut-off of the electron distribution help in finding appropriate integration borders. For MultiplicatorOfLowestBorder and MultiplicatorOfBiggestBorder a good choice is 10^(-6) and 5, respectively, to keep the relative error smaller than 0.01. RelativeError can be set to 0.1'''
    NElectrons = interp1d(ValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
    BiggestIntegrationBorder = SynchrotronCriticalxSync(ValuesForgaIt[-1])*MultiplicatorOfBiggestBorder # Notice that ValuesForgaIt[-1] = NElectronsga0
    LowestIntegrationBorder = SynchrotronCriticalxSync(ValuesForgaIt[0])*MultiplicatorOfLowestBorder # Notice that ValuesForgaIt[0] = Lowestga
    ListOfIntegrationBorders = SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
    #print('    Integration borders of SyncPsNumberDensity: ', ListOfIntegrationBorders)
    SyncPsNumberDensityResult = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        SyncPsNumberDensityResult += integrate.quad(SyncPsSpectralNumberDensity, LeftBorder, RightBorder, args=(NElectrons,ValuesForgaIt,RelativeError), epsrel=RelativeError)[0]
    return SyncPsNumberDensityResult

def SyncPsSpectralEnergyDensity(xSync,NElectrons,ValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
    '''The spectral energy-density of the SyncPs in units of 1/m^3.'''
    return xSync*SyncPsSpectralNumberDensity(xSync,NElectrons,ValuesForgaIt,RelativeError)

def SyncPsEnergyDensity(MultiplicatorOfLowestBorder,MultiplicatorOfBiggestBorder,ValuesForgaIt,ValuesForNextNElectrons,RelativeError=integrate.quad.__defaults__[3]):
    '''The total energy-density of the SyncPs in units of 1/m^3. Actually, one should integrate from 0 to infinity. However, the critical frequencies of the upper and lower cut-off of the electron distribution help in finding appropriate integration borders. For MultiplicatorOfLowestBorder and MultiplicatorOfBiggestBorder a good choice is 1 and 5, respectively, to keep the relative error smaller than 0.001. RelativeError can be set to 0.1'''
    NElectrons = interp1d(ValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
    BiggestIntegrationBorder = SynchrotronCriticalxSync(ValuesForgaIt[-1])*MultiplicatorOfBiggestBorder
    LowestIntegrationBorder = SynchrotronCriticalxSync(ValuesForgaIt[0])*MultiplicatorOfLowestBorder
    ListOfIntegrationBorders = SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
    #print('    Integration borders of SyncPsEnergyDensity: ', ListOfIntegrationBorders)
    SyncPsEnergyDensityResult = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        SyncPsEnergyDensityResult += integrate.quad(SyncPsSpectralEnergyDensity, LeftBorder, RightBorder, args=(NElectrons,ValuesForgaIt,RelativeError), epsrel=RelativeError)[0]
    return SyncPsEnergyDensityResult


## Soft photons energy injection

# The amount of energy from soft photons (LEPs) that is going into the cascade is computed now. It is called a power-density, but could also be called an energy-injection-rate. In truth, it gives the energy that is injected per unit space volume and per unit time. It makes use of the function SyncPsTotalICScatteringRate, which gives the scattering probability-rate of soft photons on the electron-field NElectrons. Although SyncPsTotalICScatteringRate is defined in dependence on xSync it can be used for any soft photon energy as well.

def LEPPowerDensityIntegrand(x,n0,NElectrons,ValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is in units of 1/(s*m^3) and gives the spectral energy-density, that is injected via LEPs per unit time, per unit space volume and per unit dimensionless energy x. SyncPsTotalICScatteringRate can be equally used for ExPs as for SyncPs.'''
    return x*n0(x)*SyncPsTotalICScatteringRate(x,NElectrons,ValuesForgaIt,RelativeError)

# Now, the case-like definition of the power-density, measured in 1/(s*m^3).
if Usedn0 == n0Delta:
    def LEPPowerDensity(n0,ValuesForgaIt,ValuesForNextNElectrons,RelativeError=integrate.quad.__defaults__[3]):
        NElectrons = interp1d(ValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
        return LEPPowerDensityIntegrand(x0,n0,NElectrons,ValuesForgaIt,RelativeError)#*me*c**2
elif Usedn0 == n0MultiDelta:
    def LEPPowerDensity(n0,ValuesForgaIt,ValuesForNextNElectrons,RelativeError=integrate.quad.__defaults__[3]):
        NElectrons = interp1d(ValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
        LEPPowerDensitySum = 0
        for i in x0MultiDelta:
            LEPPowerDensitySum += LEPPowerDensityIntegrand(i,n0,NElectrons,ValuesForgaIt,RelativeError)
        return LEPPowerDensitySum#*me*c**2
elif Usedn0==n0Exp or Usedn0==n0PL or Usedn0==n0Planck:  # In the following cases, the integration has to be performed numerically.
    def LEPPowerDensity(n0,ValuesForgaIt,ValuesForNextNElectrons,RelativeError=integrate.quad.__defaults__[3]):
        NElectrons = interp1d(ValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
        LowestIntegrationBorder = x1
        BiggestIntegrationBorder = x0
        if LowestIntegrationBorder >= BiggestIntegrationBorder:
            return 0.0
        else:
            NumberOfBorders = max(6,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
            ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
            LEPPowerDensityResult = 0.0
            for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                LEPPowerDensityResult += integrate.quad(LEPPowerDensityIntegrand, LeftBorder, RightBorder, args=(n0,NElectrons,ValuesForgaIt), epsrel=RelativeError)[0]
            return LEPPowerDensityResult#*me*c**2


## Evaluation and tests of the synchrotron-related functions:

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n\n----------- Consider the terms of the kinetic equation of synchrotron-photons -----------\n')

    print('Test for SyncPsSpectralProductionRate:', SyncPsSpectralProductionRate(0.001,CurrentNElectrons,CurrentValuesForgaIt,RelativeError=0.01),'1/(m^3*s)')
    print('Test for SyncPsICScatteringRate:      ', SyncPsICScatteringRate(0.001,10**6),'m^3/s')
    print('Test for SyncPsICScatteringRateAlt:   ', SyncPsICScatteringRateAlt(0.001,10**6),'m^3/s')
    print('Test for SyncPsTotalICScatteringRate: ', SyncPsTotalICScatteringRate(0.001,CurrentNElectrons,CurrentValuesForgaIt),'1/s')
    
    print('\nTest for SyncPsOpticalDepth with respect to IC-scattering:', SyncPsOpticalDepth(SyncPsTotalICScatteringRate,0.001,CurrentNElectrons,CurrentValuesForgaIt))
    print('Test for SyncPsOpticalDepth with respect to escape:       ', SyncPsOpticalDepth(SynchrotronPhotonEscapeRate,0.001))
    
    print('\nSynchrotron critical energy xSync of Lowestga:     ', SynchrotronCriticalxSync(CurrentLowestga))
    print('Synchrotron critical energy xSync of NElectronsga0:', SynchrotronCriticalxSync(CurrentNElectronsga0))

    def EvaluateSyncPsSpectralProductionRateAndEnergyDensity():
        if isinstance(CurrentNElectrons, interp1d):
            ValuesForxSync = np.logspace(np.log10(SynchrotronCriticalxSync(Lowestga))-3,np.log10(SynchrotronCriticalxSync(Highestga))+1,21)  # A sampling-range for the dimensionless SyncP-energy. The extent of these distributions roughly is from the critical xSync at the lowest ga, where NElectrons is non-vanishing, (minus three orders of magnitude) to the highest ga, where NElectrons is non-vanishing (, plus one order of magnitude).
        elif isinstance(CurrentNElectrons, types.FunctionType):
            ValuesForxSync = np.logspace(np.log10(SynchrotronCriticalxSync(Injectedga1))-3,np.log10(SynchrotronCriticalxSync(Injectedga0))+1,21)
        print('Sample values for xSync:\n', ValuesForxSync)
        ValuesForSyncPsSpectralProductionRate = [SyncPsSpectralProductionRate(i,CurrentNElectrons,CurrentValuesForgaIt) for i in ValuesForxSync]
        print('Corresponding values for SyncPsSpectralProductionRate:\n', ValuesForSyncPsSpectralProductionRate)
        ValuesForSyncPsSpectralEnergyDensity = [SyncPsSpectralEnergyDensity(i,CurrentNElectrons,CurrentValuesForgaIt) for i in ValuesForxSync]
        print('Corresponding values for SyncPsSpectralEnergyDensity:\n', ValuesForSyncPsSpectralEnergyDensity)
        pl.rc('font', family='serif')
        Fig, LeftAxis = pl.subplots(figsize=(19, 11), num="KESP: SyncPs' spectral production-rate and energy-density versus SyncPs' dimensionless energy")
        LeftAxis.set_xlabel("Synchrotron-photons dimensionless energy $x_{\mathrm{sync}}$", fontsize=35)        
        #LeftAxis.set_xlim((,))
        LeftAxis.set_ylabel("Spectral production-rate times\nSyncPs dim.less energy $x_{\mathrm{sync}}$ in s$^{-1} \cdot $m$^{-3}$", fontsize=35, color='red')
        LeftAxis.loglog(ValuesForxSync,ValuesForSyncPsSpectralProductionRate*ValuesForxSync, label="Spectral production rate for $B = %s \mathrm{T}$" % B4, color='red', linewidth=4)
        LeftAxis.tick_params(axis='y', labelsize=28, direction='inout', width=3, length=10, labelcolor='red')
        LeftAxis.tick_params(axis='y', which='minor', labelsize=22, direction='inout', width=2, length=6, labelcolor='red')
        LeftAxis.tick_params(axis='x', labelsize=28, direction='inout', width=3, length=10, pad=7)
        LeftAxis.tick_params(axis='x', which='minor', labelsize=22, direction='inout', width=2, length=6, pad=7)
        RightAxis = LeftAxis.twinx()  # instantiate a second axes that shares the same x-axis.
        RightAxis.set_ylabel("Spectral energy-density in m$^{-3}$", fontsize=35, color='green')
        #RightAxis.set_xlim((,))
        RightAxis.loglog(ValuesForxSync,ValuesForSyncPsSpectralEnergyDensity, label='Spectral number-density', color='green', linewidth=4, linestyle=':')
        RightAxis.tick_params(axis='y', labelsize=28, direction='inout', width=3, length=10, labelcolor='green')
        RightAxis.tick_params(axis='y', which='minor', labelsize=22, direction='inout', width=2, length=6, labelcolor='green')
        #LeftAxis.legend(loc="upper left", fontsize=22)
        #RightAxis.legend(loc="lower right", fontsize=22)
        #Fig.tight_layout()  # otherwise the right y-label is slightly clipped
        Fig.subplots_adjust(top=0.96, bottom=0.13, left=0.14, right=0.88)
        pl.savefig("Run %s, KESP - SyncPs Dotn times xSync and spectral energy-density versus photon energy from %s.svg" % (RunIdentifier,GetCurrentDate()))

    def EvaluateSyncPsSpectralICScatteringRate(Testx):
        # Look at SyncPsSpectralICScatteringRate. It will be integrated along xga, so plot it along xga:
        for Testga in [1.0/Testx,0.5/Testx,0.2/Testx,0.1/Testx,0.01/Testx]:#[1.0/Testx,2.0/Testx,4.0/Testx,10.0/Testx,100.0/Testx,1000.0/Testx]:
            ValuesForxga=np.logspace(-3,np.log10(xmaxOfA20(Testga,Testx)),1000)
            ValuesForSyncPsSpectralICScatteringRate=SyncPsSpectralICScatteringRate(ValuesForxga,Testx,Testga)
            #print(ValuesForSyncPsSpectralICScatteringRate)
            pl.figure(num="IC-scattering spectral rate")
            pl.loglog(ValuesForxga,ValuesForSyncPsSpectralICScatteringRate, label="$\gamma = %.2f$, $x = %.7f$" % (Testga, Testx))
            pl.xscale('log')
            pl.xlabel("High-energetic photon energy $x_{\gamma}$")
            pl.ylabel("SyncPsSpectralICScatteringRate")
            pl.legend(loc="best")
            #pl.scatter(xmaxOfA20(Testga,Testx),0.0) # This is the biggest integration border along xga.
        
    def EvaluateSyncPsICScatteringRate():
        # Look at SyncPsICScatteringRate. It will be integrated along ga, so plot it along ga:
        # StartingTime = time.time()
        for Testx in [0.01,0.001,0.0001]:
            ValuesForga=np.logspace(np.log10(0.002/Testx),np.log10(10000000.0/Testx),500)
            ValuesForSyncPsICScatteringRate=[SyncPsICScatteringRate(Testx,i) for i in ValuesForga]
            #print(ValuesForSyncPsICScatteringRate)
            pl.figure(num="IC-scattering rate")
            pl.loglog(ValuesForga,ValuesForSyncPsICScatteringRate, label="$x = %.7f$,\nLowerIntBorder=0.001*BiggestIntBorder,\nSampleIntegrationBorders8" % Testx)
            pl.xlabel("Electron energy $\gamma$")
            pl.ylabel("SyncPsICScatteringRate")
            pl.legend(loc="best")
        # pl.savefig("SyncP IC-scattering - SyncPsICScatteringRate from %s.svg" % GetCurrentDate())
        # EndingTime = time.time()
        # TimeInterval = (EndingTime-StartingTime)
        # print("\nTime spent:", TimeInterval, "s\n") # Evaluation of SyncPsICScatteringRate is about 10 % slower than evaluation of SyncPsICScatteringRateAlt.

    def EvaluateSyncPsSpectralICScatteringRateAlt(Testx):
        # Look at SyncPsSpectralICScatteringRateAlt. It will be integrated along gaP, so plot it along gaP:
        for Testga in [1.0/Testx,0.5/Testx,0.2/Testx,0.1/Testx,0.05/Testx,0.02/Testx]:#[1.0/Testx,2.0/Testx,4.0/Testx,10.0/Testx,100.0/Testx,1000.0/Testx]:
            ValuesForgaP=np.logspace(np.log10(gaPminOf6(Testx,Testga)),np.log10(Testga),10000)
            ValuesForSyncPsSpectralICScatteringRateAlt=[SyncPsSpectralICScatteringRateAlt(i,Testx,Testga) for i in ValuesForgaP]
            #print(ValuesForSyncPsSpectralICScatteringRateAlt)
            pl.figure(num="IC-scattering spectral rate (alternative function)")
            pl.loglog(ValuesForgaP,ValuesForSyncPsSpectralICScatteringRateAlt, label="$\gamma = %.2f$, $x = %.7f$" % (Testga, Testx))
            pl.xlabel("Final electron energy $\gamma'$")
            pl.ylabel("SyncPsSpectralICScatteringRateAlt")
            pl.legend(loc="best")
            #pl.scatter(xmaxOfA20(Testga,Testx),0.0) # This is the biggest integration border along xga.

    def EvaluateSyncPsICScatteringRateAlt():
        # Look at SyncPsICScatteringRateAlt. It will be integrated along ga, so plot it along ga:
        # StartingTime = time.time()
        for Testx in [0.01,0.001,0.0001]:
            ValuesForga=np.logspace(np.log10(0.002/Testx),np.log10(10000000.0/Testx),500)
            ValuesForSyncPsICScatteringRateAlt=[SyncPsICScatteringRateAlt(Testx,i) for i in ValuesForga]
            #print(ValuesForSyncPsICScatteringRateAlt)
            pl.figure(num="IC-scattering rate (alternative integration)")
            pl.loglog(ValuesForga,ValuesForSyncPsICScatteringRateAlt, label="$x = %.7f$,\nSampleIntegrationBorders7" % Testx)
            pl.xlabel("Electron energy $\gamma$")
            pl.ylabel("SyncPsICScatteringRateAlt")
            pl.legend(loc="best")
        # pl.savefig("SyncP IC-scattering - SyncPsICScatteringRateAlt from %s.svg" % GetCurrentDate())
        # EndingTime = time.time()
        # TimeInterval = (EndingTime-StartingTime)
        # print("\nTime spent:", TimeInterval, "s\n") # Evaluation of SyncPsICScatteringRateAlt is about 10 % faster than evaluation of SyncPsICScatteringRate.

def EvaluateSyncPsLossRate(ValuesForgaIt,ValuesForNextNElectrons):
    print('\nCalling EvaluateSyncPsLossRate:')
    InternalInterpolatedObjectForValuesForNextNElectrons = interp1d(ValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
    CurrentNElectrons = InternalInterpolatedObjectForValuesForNextNElectrons
    ValuesForxSync = np.logspace(-7,-1,40)  # A sampling-range for xSync, which ranges over reasonable values. Logarithmic division of the interval is used.
    ValuesForSyncPsTotalICScatteringRate = [SyncPsTotalICScatteringRate(i,CurrentNElectrons,ValuesForgaIt) for i in ValuesForxSync]
    # print('Corresponding values for SyncPsTotalICScatteringRate:\n', ValuesForSyncPsTotalICScatteringRate)
    ValuesForSynchrotronPhotonEscapeRate = [SynchrotronPhotonEscapeRate(i) for i in ValuesForxSync] # The escape-probability-rate.
    ValuesForSyncPsLossRate = [SyncPsLossRate(i,CurrentNElectrons,ValuesForgaIt) for i in ValuesForxSync]
    pl.rc('font', family='serif')
    Figure = pl.figure(figsize=(18, 10), num="KESP: Synchrotron-photons' spectral loss-rate versus synchrotron-photons' dimensionless energy")
    # On the left y-axis, the probability-rate of scattering-events is drawn:
    LeftyAxis = Figure.add_subplot(111)
    LeftyAxis.loglog(ValuesForxSync,ValuesForSyncPsLossRate, label='Spectral loss rate (sum)', linewidth=2)
    LeftyAxis.loglog(ValuesForxSync,ValuesForSynchrotronPhotonEscapeRate, label='Spectral escape rate', linewidth=2)
    LeftyAxis.loglog(ValuesForxSync,ValuesForSyncPsTotalICScatteringRate, label='Spectral IC up-scattering rate', linewidth=2)
    LeftyAxis.xaxis.set_tick_params(labelsize=30, direction='inout', width=3, length=10, pad=7)
    LeftyAxis.yaxis.set_tick_params(labelsize=30, direction='inout', width=3, length=10, right=False)
    pl.legend(loc="best", fontsize=18)
    pl.xlabel("Synchrotron-photons' dimensionless energy $x_{\mathrm{sync}}$", fontsize=35)
    pl.ylabel("Probability rate in s$^{-1}$", fontsize=35)
    Ticks=LeftyAxis.get_yticks() # Store the locations of the ticks.
    LeftyAxis.set_ylim(Ticks[1],Ticks[-2]) # This setting of the limits might change from plot to plot. Usage for pion-decay: Ticks[1],Ticks[7]
    # The right y-axis shows the optical depth, corresponding to the probability-rate.
    RightyAxis = Figure.add_subplot(111, sharex=LeftyAxis, frameon=False)
    RightyAxis.loglog()
    RightyAxis.yaxis.tick_right()
    RightyAxis.yaxis.set_label_position("right")
    RightyAxisTicks=np.asarray([10**(-3),10**(-2),10**(-1),10**(0),10**(1),10**(2),10**(3)]) # The values of the optical depth, that are supposed to be marked with ticks.
    RightyAxis.yaxis.set_ticks(RightyAxisTicks/MeanEscapeTime) # The values of the spectral loss-rate, that correspond to the ticked optical depth. This essentially is loss-probability-rate-of-process = tau / MeanEscapeTime, where tau represents RightyAxisTicks.
    RightyAxis.yaxis.set_ticklabels(RightyAxisTicks)
    RightyAxis.set_ylim(Ticks[1],Ticks[-2]) # Usage for pion-decay: Ticks[1],Ticks[7]
    pl.ylabel("Corresponding optical depth", fontsize=35)
    RightyAxis.xaxis.set_tick_params(which='both', bottom=False, top=False, labelbottom=False)
    RightyAxis.yaxis.set_tick_params(labelsize=30, direction='inout', width=3, length=10)
    pl.subplots_adjust(top=0.96, bottom=0.13, left=0.12, right=0.88)
    pl.savefig("Run %s, KESP - SyncPs spectral loss-rate versus photon energy from %s.svg" % (RunIdentifier,GetCurrentDate()))


## Super-iteration

# The super-iteration is an iteration of iterations. It determines the steady-state electron distribution NElectrons and synchrotron-photon distribution nSyncPs.

ValuesForxSyncSuperIt = UsedSampleIntBordersOfxSync(CurrentxSync1,CurrentxSync0Used)[0] # The values along which nSyncPs will be sampled.

SuperItCounter = 0 # A number that counts the walk-throughs of the super-iteration.

nSyncPsSuperIt = {} # A data-structure to save the findings of the SyncP-distribution of the super-iteration. It will be a dictionary of dictionaries. The outer dictionary has the number of the super-iteration-step (super-iteration counter) as keys. The value of each key is a dictionary, which will contain the keys 'ValuesForxSyncSuperIt' and 'ValuesFornSyncPs'. 'ValuesForxSyncSuperIt' contains an array with the values of xSync that were used for the respective super-iteration-step. 'ValuesFornSyncPs' contains an array with the values of the spectral SyncP number-density that were obtained in the respective super-iteration-step.
# Initialise nSyncPsSuperIt:
nSyncPsSuperIt[SuperItCounter] = {}
nSyncPsSuperIt[SuperItCounter]['ValuesForxSyncSuperIt'] = ValuesForxSyncSuperIt
ValuesFornSyncPs = np.asarray([SyncPsSpectralNumberDensity(xSync,StartingNElectrons,CurrentValuesForgaIt,RelativeError=0.01) for xSync in ValuesForxSyncSuperIt]) # The values of the zeroth entry are based on the starting function of NElectrons.
nSyncPsSuperIt[SuperItCounter]['ValuesFornSyncPs'] = ValuesFornSyncPs

nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0) # Create an interpolation that can be called like a function.

NElectronsSuperIt = {} # A data-structure to save the findings of the electron-distribution of the super-iteration. It will be a dictionary of dictionaries. The outer dictionary has the number of super-iteration-step as keys. The value of each key will be a NElectronsIterated-dictionary, which will contain all the results of one common iteration.

if "MainProcess" in multiprocessing.current_process().name and UsedIterationScheme == 'PointsFromRToLIterationPointwise': # The first condition prevents evaluation in child-processes.
    def EvaluateSuperIteration(CurrentNElectrons, ValuesForgaIt, LastValuesForgaIt, CurrentValuesForNElectrons, NElectronsSuperIt, NElectronsga0, Lowestga, SuperItCounter, nSyncPs, ValuesForxSyncSuperIt, ValuesFornSyncPs, nSyncPsSuperIt, xSync0, xSync1, xSync0Used, xTotal0, ResumeSuperIt=False):
        '''This function takes turns at determining the spectral number-density nSyncPs of SyncPs (via evaluating SyncPsSpectralNumberDensity) and determining the spectral number-density of electrons (by evaluating EvaluateIteration).'''
        while (SuperItCounter<2 and ResumeSuperIt==False) or (SuperItCounter<3 and ResumeSuperIt==True): # The first test is for the case of a SuperIt, that is running from its beginning. The second test is for the case of a resumed SuperIt. Alternatively, one can invent a function like MoreSuperItsNecessary(), which can decide whether convergence is achieved or not:
            if ResumeSuperIt==False:
                SuperItCounter += 1 # Increment the counter.
            # First, do an iteration of the electron-distribution:
            if SuperItCounter==1:
                CurrentNElectrons, CurrentValuesForNElectrons, NElectronsSuperIt = EvaluateIteration(CurrentNElectrons, ValuesForgaIt, CurrentValuesForNElectrons, NElectronsSuperIt, [1000, 0.001], [0.0, 0.0], SuperItCounter) # Perform a common iteration with all the initial values and without inclusion of SyncPs as soft background photons. They are at first not included because the SyncP-field of StartingNElectrons is so weak in comparison to Usedn0, that it is not worth the effort. Usage of ValuesForxSyncSuperIt=[1000, 0.001], and ValuesFornSyncPs=[0.0, 0.0] mimics that the SyncP-field vanishes.
            else:
                NElectronsga0 = DetermineNElectronsga0(xSync0Used) # Determine NElectronsga0 anew.
                Lowestga = DetermineLowestga(xSync0Used) # Determine Lowestga anew.
                ValuesForgaIt = CompileSamplingRange(Lowestga,NElectronsga0,xSync0Used) # Sample ValuesForgaIt anew.
                LastValuesForgaIt = ValuesForgaIt # This is necessary for some functions to be callable also in this 'PointsFromRToLIterationPointwise' case.
                if ResumeSuperIt==False:
                    CurrentNElectrons, CurrentValuesForNElectrons, NElectronsSuperIt = EvaluateIteration(CurrentNElectrons, ValuesForgaIt, CurrentValuesForNElectrons, NElectronsSuperIt, ValuesForxSyncSuperIt, ValuesFornSyncPs, SuperItCounter) # Perform a common iteration with inclusion of SyncPs as soft background photons.
                elif ResumeSuperIt==True:
                    CurrentNElectrons, CurrentValuesForNElectrons, NElectronsSuperIt = EvaluateIteration(CurrentNElectrons, ValuesForgaIt, CurrentValuesForNElectrons, NElectronsSuperIt, ValuesForxSyncSuperIt, ValuesFornSyncPs, SuperItCounter, ResumeSuperIt) # Resume a common iteration that has partly been performed earlier.
                    ResumeSuperIt=False # Make sure that in subsequent SuperIt-steps, the usual common iteration is performed.
            EvaluatePlotForIteration(ValuesForgaIt, NElectronsSuperIt, CurrentNElectrons, SuperItCounter, xSync0Used) # Plot the results of this common iteration.
            # Second, consider the SyncPs:
            xSync0 = DeterminexSync0(NElectronsga0) # Determine xSync0 anew.
            xSync1 = DeterminexSync1(Lowestga) # Determine xSync1 anew.
            xSync0Used = DeterminexSync0Used(xSync0) # Determine xSync0Used anew.
            xTotal0 = DeterminexTotal0(xSync0Used) # Determine xTotal0 anew.
            ValuesForxSyncSuperIt = UsedSampleIntBordersOfxSync(xSync1,xSync0Used)[0] # Sample the values of nSyncPs anew.
            nSyncPsSuperIt[SuperItCounter] = {} # Save...
            nSyncPsSuperIt[SuperItCounter]['ValuesForxSyncSuperIt'] = ValuesForxSyncSuperIt # ...the values of xSync, ...
            ValuesFornSyncPs = np.asarray([SyncPsSpectralNumberDensity(xSync,CurrentNElectrons,ValuesForgaIt) for xSync in ValuesForxSyncSuperIt]) # ...compute the synchrotron-spectrum...
            nSyncPsSuperIt[SuperItCounter]['ValuesFornSyncPs'] = ValuesFornSyncPs # ...and save its values, too.
            nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0) # Create an interpolation that can be called like a function.
            ExportDataOfnSyncPs(nSyncPsSuperIt, SuperItCounter) # Store the SyncP-field in a file...
            EvaluatePlotFornSyncPs(nSyncPsSuperIt, SuperItCounter) # ... and plot them.
            if ResumeSuperIt==True:
                ResumeSuperIt=False # This is necessary that the usual SuperIt-path (with incrementation of the SuperItCounter) is taken after the SuperIt was resumed.
        return CurrentNElectrons, ValuesForgaIt, LastValuesForgaIt, CurrentValuesForNElectrons, NElectronsSuperIt, NElectronsga0, Lowestga, SuperItCounter, nSyncPs, ValuesForxSyncSuperIt, ValuesFornSyncPs, nSyncPsSuperIt, xSync0, xSync1, xSync0Used, xTotal0

def ExportDataOfnSyncPs(nSyncPsSuperIt, SuperItCounter):
    '''Both the value of SuperItCounter and the complete dictionary nSyncPsSuperIt are written as a string into a .dat-file.'''
    np.set_printoptions(threshold=np.inf, precision=20)
    Outputfile = open("Run %s, nSyncPs versus xSync from %s (SuperIt=%s, data).dat" % (RunIdentifier,GetCurrentDate(),SuperItCounter), 'w')
    Outputfile.write('%s' % nSyncPsSuperIt)
    Outputfile.close()

def EvaluatePlotFornSyncPs(nSyncPsSuperIt, SuperItCounter): # Now, plotting of the super-iterated curves.
    pl.figure(figsize=(18, 14), num="KESP: nSyncPs(xSync) versus photon-energy xSync (Super-Iteration counter = %s)" % SuperItCounter)
    pl.title('Super-iteration counter = %s' % SuperItCounter, fontsize=34)
    pl.rc('font', family='serif')
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=28, direction='inout', width=3, length=10, pad=7)
    ax.yaxis.set_tick_params(labelsize=28, direction='inout', width=3, length=10)
    ax.xaxis.set_tick_params(which='minor', labelsize=22, direction='inout', width=2, length=6, pad=7)
    ax.yaxis.set_tick_params(which='minor', labelsize=22, direction='inout', width=2, length=6)
    #pl.xlim(Lowestga*x0,Highestga*x0*1.2)
    pl.xlabel("Dimensionless SyncP-energy $x_{sync}$", fontsize=34)
    pl.ylabel("Spectral number density $n(x_{sync})$ times\ndim.less SyncP-energy $x_{sync}$ in $\mathrm{m}^{-3}$", fontsize=34)
    for i in range(SuperItCounter+1):
        if i!=0 or StartingNElectrons != NElectronsZero:
            pl.loglog(nSyncPsSuperIt[i]['ValuesForxSyncSuperIt'],nSyncPsSuperIt[i]['ValuesForxSyncSuperIt']*nSyncPsSuperIt[i]['ValuesFornSyncPs'], label='Super-iteration counter = %s' % i)
    #pl.ylim(10.0**(5),2*10.0**(8.0))
    pl.legend(loc="best", fontsize=26)
    pl.subplots_adjust(top=0.95, bottom=0.11, left=0.13, right=0.96)
    pl.savefig("Run %s, nSyncPs versus xSync from %s (SuperIt=%s, plot).svg" % (RunIdentifier,GetCurrentDate(),SuperItCounter))


## Tests of the IC-scattering- and pair-production-related functions:

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.

    print('\n\n-------------------------- Consider inverse-Compton-scattering --------------------------\n')
    
    print('Test for COfA1:                  ', COfA1(100000.0,10000.0,Usedn0,nSyncPs,ValuesForxSyncSuperIt), '1/s')
    if IncludeSynchrotron:
        print('Test for COfA1WithSynchrotron:   ', COfA1WithSynchrotron(100000.0,10000.0,nSyncPs,ValuesForxSyncSuperIt), '1/s')
    print('Test for COfA1WithoutSynchrotron:', COfA1WithoutSynchrotron(100000.0,10000.0,Usedn0), '1/s')
    
    if not(IncludeSynchrotron):
        print('\nTest for COfA1DependingOnxga:', COfA1DependingOnxga(10000.0,1000.0,Usedn0,None,None), '1/s')
    elif IncludeSynchrotron:
        print('\nTest for COfA1DependingOnxga:', COfA1DependingOnxga(10000.0,1000.0,Usedn0,nSyncPs,ValuesForxSyncSuperIt), '1/s')
    
    print('\nTest for COfA21Alternative:       ', COfA21Alternative(gamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt), '1/s')
    if not(IncludeSynchrotron):
        print('Test for COfA21:                  ', COfA21(gamma,Usedn0,None,None), '1/s')
    elif IncludeSynchrotron:
        print('Test for COfA21:                  ', COfA21(gamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt), '1/s')
        print('Test for COfA21WithSynchrotron:   ', COfA21WithSynchrotron(gamma,nSyncPs,ValuesForxSyncSuperIt), '1/s')
    print('Test for COfA21WithoutSynchrotron:', COfA21WithoutSynchrotron(gamma,Usedn0), '1/s')
    
    print('\nTest for OpticalDepthOfElectrons with respect to escape:                       ', OpticalDepthOfElectrons(ElectronEscapeRate,gamma))
    print('Test for OpticalDepthOfElectrons with respect to IC-scattering on soft photons:', OpticalDepthOfElectrons(COfA21WithoutSynchrotron,gamma,Usedn0))
    if IncludeSynchrotron:
        print('Test for OpticalDepthOfElectrons with respect to synchrotron-emission:         ', OpticalDepthOfElectrons(SynchrotronSpectralLossRate,gamma))
        print('Test for OpticalDepthOfElectrons with respect to IC-scattering (total):        ', OpticalDepthOfElectrons(COfA21,gamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt))
        print('Test for OpticalDepthOfElectrons with respect to IC-scattering (SyncPs):       ', OpticalDepthOfElectrons(COfA21WithSynchrotron,gamma,nSyncPs,ValuesForxSyncSuperIt))
        print('Test for OpticalDepthOfElectrons with respect to all processes:                ', OpticalDepthOfElectrons(ElectronLossRate,gamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt))
    if not(IncludeSynchrotron):
        print('Test for OpticalDepthOfElectrons with respect to IC-scattering (total):        ', OpticalDepthOfElectrons(COfA21,gamma,Usedn0,None,None))
        print('Test for OpticalDepthOfElectrons with respect to all processes:                ', OpticalDepthOfElectrons(ElectronLossRate,gamma,Usedn0,None,None))

    print('\n\n------------------------------ Consider pair-production ---------------------------------\n')
    
    print('Test for POfB1:                  ', POfB1(500000.0,1000000.0,Usedn0,nSyncPs,ValuesForxSyncSuperIt), '1/s')
    if IncludeSynchrotron:
        print('Test for POfB1WithSynchrotron:   ', POfB1WithSynchrotron(500000.0,1000000.0,nSyncPs,ValuesForxSyncSuperIt), '1/s')
    print('Test for POfB1WithoutSynchrotron:', POfB1WithoutSynchrotron(500000.0,1000000.0,Usedn0), '1/s')
    
    print('\nTest for POfB11Alternative:       ', POfB11Alternative(10000000000000.0,Usedn0,nSyncPs,ValuesForxSyncSuperIt), '1/s')
    if not(IncludeSynchrotron):
        print('Test for POfB11:                  ', POfB11(10000000000000.0,Usedn0,None,None), '1/s')
    elif IncludeSynchrotron:
        print('Test for POfB11:                  ', POfB11(10000000000000.0,Usedn0,nSyncPs,ValuesForxSyncSuperIt), '1/s')
        print('Test for POfB11WithSynchrotron:   ', POfB11WithSynchrotron(10000000000000.0,nSyncPs,ValuesForxSyncSuperIt), '1/s')
    print('Test for POfB11WithoutSynchrotron:', POfB11WithoutSynchrotron(10000000000000.0,Usedn0), '1/s')
    
    print('\nTest for OpticalDepthOfHEPhotons with respect to escape:         ', OpticalDepthOfHEPhotons(HEPhotonEscapeRate,xgamma))
    if IncludeSynchrotron:
        print('Test for OpticalDepthOfHEPhotons with respect to pair-production:', OpticalDepthOfHEPhotons(POfB11,xgamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt))
        print('Test for OpticalDepthOfHEPhotons with respect to all processes:  ', OpticalDepthOfHEPhotons(HEPhotonLossRate,xgamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt))
    if not(IncludeSynchrotron):
        print('Test for OpticalDepthOfHEPhotons with respect to pair-production:', OpticalDepthOfHEPhotons(POfB11,xgamma,Usedn0,None,None))
        print('Test for OpticalDepthOfHEPhotons with respect to all processes:  ', OpticalDepthOfHEPhotons(HEPhotonLossRate,xgamma,Usedn0,None,None))
