## ------------------------------------------------------------------------------------------  
## Toy-model connecting vacuum-gap-models, an AGN-jet-model and an SSC-model (toy-model 1)
##                                                                                           
## Implementation of various estimating formulae (originally intended for the case of IC310)            
## ------------------------------------------------------------------------------------------

# This is the updated version of "Estimates_5_after_MiniWorkshop.py" and now it is incorporated into the complete framework of code. 

import os
import multiprocessing

from part1 import * # Import of the file, whose path was added above.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n')
    print('\n_________________________________________________________________________________________\n')
    print('---------------------    Vacuum-gaps, AGN-jets and SSC-model    -------------------------')
    print('_________________________________________________________________________________________\n\n')


## Blackbody laws

def BBEnergyOfMaximumPower(T):
    '''This yields the energy of maximum spectral power of an emitting blackbody according to Wien's displacement law in units of the electron's rest-energy.'''
    return (5.87893*10**10*T*h)/(me*c**2)

def BBTotalNumberDensity(T):
    '''This gives the total number-density of photons of a Planck-distribution in units of 1/m^3.'''
    return 1.202*16*np.pi*(kB*T/h/c)**3

def BBSpectralNumberDensity(epsilon,T):
    '''This gives the spectral number-density of photons of a Planck-distribution in units of 1/J/m^3. epsilon is the photon-energy in J.'''
    return (8*np.pi*epsilon**2)/(c**3*h**3*(np.exp(epsilon/(kB*T))-1))

def BBTotalEnergyDensity(T):
    '''This gives the total energy-density of photons of a Planck-distribution in units of J/m^3.'''
    return aRadiation*T**4

def BBTotalEnergyDensityDiluted(T,DistanceGap,Toy1Distance,M9):
    '''This gives the total energy-density in units of J/m^3 of blackbody photons at distance Toy1Distance (in Schwarzschild-radii) where T is the effective temperature of the emitting gas in units of of K. DistanceGap is the distance in units of Schwarzschild-radii at which the energy-density has the value BBTotalEnergyDensity(T). It is for my estimating toy-model for Dorit.'''
    return BBTotalEnergyDensity(T)*(rSchwarzschildTom(DistanceGap,M9)**2)/(rSchwarzschildTom(Toy1Distance,M9)**2)

def BBMeanEnergy(T):
    '''This gives the mean energy of photons of a Planck-distribution in units of eV.'''
    return 2.7*kB*T/e


## The effective temperature of an optically thick accretion disk

def TeffOfEquation5ByKadowaki(M9,etaff,Dotm,ADAFrInner):
    '''This is equation 5 of 2015ApJ...802..113K solved for T_eff. So, this gives the effective temperature of the accretion disk at its surface in units of K.'''
    return ((4*3*G*ToThe9SolarMassesTokg(M9)*DotMEddingtonTokgPers(M9,etaff,Dotm)*(1-np.sqrt(3*rSchwarzschild(M9)/rSchwarzschildTom(ADAFrInner,M9))))/(c*8*np.pi*rSchwarzschildTom(ADAFrInner,M9)**3*aRadiation))**(1/4)

def TeffOfEquation3dot50ByKato(M9,etaff,Dotm,ADAFrInner,r):
    '''This is equation 3.50 of 2008bhad.book.....K solved for T_eff. So, this gives the effective temperature of the accretion disk at its surface in units of K.'''
    return ((4*3*G*ToThe9SolarMassesTokg(M9)*DotMEddingtonTokgPers(M9,etaff,Dotm)*(1-np.sqrt(rSchwarzschildTom(ADAFrInner,M9)/rSchwarzschildTom(r,M9))))/(c*8*np.pi*rSchwarzschildTom(r,M9)**3*aRadiation))**(1/4)

def TvirOfEquation3dot95ByKato(M9,r):
    '''This is equation 3.95 of 2008bhad.book.....K. So, this gives the virial temperature of the accretion disk in units of K.'''
    return (G*ToThe9SolarMassesTokg(M9)*mProton)/(kB*rSchwarzschildTom(r,M9))


## Equations according to 2011ApJ...730..123L (Levinson and Rieger, 2011)

def LRMagneticFluxDensityOfEquation3(M9,Dotm,r):
    '''This is an approximation for the magnetic flux-density inside the accretion disk. r is the radial distance to the black hole in units of Schwarzschild-radii. It is in units of Gauss.'''
    return 4.0*10**4*(Dotm/M9)**0.5*r**(-1.25)

def LRxgamma(TElectron):
    '''This is an approximation for the energy of the HE-photons. It is in units of the electron's rest-energy (in contrast to L&R). Therefore, x is used instead of epsilon.'''
    return 3.0*Theta(TElectron)

def LRTotalPhotonDensityRIAF(Dotm,M9):
    '''Equation 5 of 2011ApJ...730..123L (Levinson and Rieger, 2011). It gives the total number-density of bremsstrahlung-photons that are present in a black hole magnetosphere and stemming from a RIAF. Probably, L&R give the equation in units of cm^-3. Hence, here it is in units of m^-3 (due to the 10^6 conversion)'''
    return 1.4*10**11*Dotm**2*M9**(-1)*10**6

def LRTotalInjectionRatePPRIAF(Dotm,M9):
    '''Between equation 5 and eq. 6, L&R give a total number of particles that are pair-produced per unit time interval in a certain space-volume by the RIAF-photons. Similarly, here we determine the corresponding total number-density of particles that are pair-produced per unit time interval by the RIAF-photons. Based on LRTotalPhotonDensityRIAF, it might be in units of 1/(m^3*s). The approximation for the pair-production cross-section is from page 4 of L&R and is also used elsewhere.'''
    return 0.2*sigmaT*c*LRTotalPhotonDensityRIAF(Dotm,M9)**2

def LRCoefficientOfEquation13():
    '''This is the numerical value of equation 13 according to my computation. It is in accordance with L&R. Cf. "2011ApJ...730..123L - Inferring equation 13.png".'''
    print("SI-units: ",((3*1.7*10**21*G*10**9*2*10**30*4*np.pi*epsilon0)/(e*c*c))**0.25)
    print("Gaussian units: ",((3*1.7*10**21*GGaussian*10**9*2*10**33/299.792458)/(eGaussian*cGaussian*cGaussian))**0.25)

def LRCoefficientOfEquation14():
    '''This is the numerical value of equation 14 according to my computation. It is in accordance with L&R. Cf. "2011ApJ...730..123L - Inferring equation 14.png"'''
    print("SI-units: ",((c**2*e*1.7*10**21*10**(-6))/(2*G*sigmaT*3*10**9*2*10**30*10**(-7)))**0.5)
    print("Gaussian units: ",((cGaussian**2*eGaussian*1.7*10**21/299.792458)/(2*GGaussian*sigmaTGaussian*3*10**9*2*10**33))**0.5)

def LRgammaCurOfEquation13(B4,M9,hGap,rCurvature):
    '''This determines the equilibrium gamma due to curvature-radiation. B4 might be the magnetic flux-density in 10**4 Gauss, M9 is the mass of the black hole in 10**9 solar masses, hGap is the gap height in Schwarzschild-radii, rCurvature is the curvature-radius of the magnetic field-lines in units of Schwarzschild-radii.'''
    return 5*10**10*B4**0.25*M9**0.5*hGap**0.25*rCurvature**0.5

def LRgammaICOfEquation14(B4,hGap,uSoft):
    '''This determines the equilibrium gamma due to inverse-Compton-scattering. B4 might be the magnetic flux-density in 10**4 Gauss, hGap is the gap height in Schwarzschild-radii, uSoft is the energy-density of the soft background photons in units of erg/cm**3. This is not exactly equation 14 but the equivalent equation using uSoft.'''
    return np.sqrt(3)*2*10**9*B4**0.5*hGap**0.5*uSoft**(-0.5)

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('Estimate for total number-density of RIAF-photons (bremsstrahlung) of L&R:', LRTotalPhotonDensityRIAF(Dotm,M9), '1/m^3')
    print('Estimate for total injection-rate of pair-produced particles of L&R:', LRTotalInjectionRatePPRIAF(Dotm,M9), '1/(m^3/s)')

## Equations concerning the size of the blob

def Toy1SizeOfBlob(Toy1Distance,DistanceGap,SizeGap,thetaOpening):
    '''This gives the size (that is the diameter of the blob/jet in units of Schwarzschild-radii. SizeGap is the size of the gap in units of Schwarzschild-radii. thetaOpening is the half opening angle of the jet. Toy1Distance and DistanceGap is the distance of the blob and the distance of the gap to the centre in units of Schwarzschild-radii.'''
    return (SizeGap+2.0*(Toy1Distance-DistanceGap)*np.tan(thetaOpening))

def Toy1RadiusOfBlob(SizeGap,Toy1Distance,DistanceGap,thetaOpening,M9):
    '''This gives the approximated (that is for big distances) radius of the blob in units of cm. SizeGap is the size of the gap in units of Schwarzschild-radii. thetaOpening is the half opening angle of the jet. Toy1Distance is the distance of the blob to the centre in units of Schwarzschild-radii.'''
    return Toy1SizeOfBlob(Toy1Distance,DistanceGap,SizeGap,thetaOpening)/2*rSchwarzschild(M9)*100


## Equations concerning the magnetic field

def Toy1MagneticFluxDensityDipolarDilution(M9,Dotm,DistanceCoinc,Toy1Distance,theta):
    '''This is an approximation for the magnetic flux-density according to my estimation. It is in units of Gauss. theta is the polar angle, Toy1Distance is the distance to the black hole in units of Schwarzschild-radii and DistanceCoinc is the radial distance from the centre (black hole) to the point, at which both expressions for the flux-density coincide, in units of Schwarzschild-radii.'''
    return (LRMagneticFluxDensityOfEquation3(M9,Dotm,DistanceCoinc*np.sin(np.pi/8))*DistanceCoinc**3)/(np.sqrt(1.0+3*(np.cos(np.pi/8))**2)*Toy1Distance**3)*np.sqrt(1.0+3*(np.cos(theta))**2)

def Toy1MagneticFluxDensityDipolarDilutionOnAxis(M9,Dotm,DistanceCoinc,Toy1Distance):
    '''This is an approximation for the magnetic flux-density in the jet / on the axis of symmetry (i. e. for theta = 0) according to my estimation for Dorit. It is in units of Gauss. Toy1Distance is the distance to the black hole in units of Schwarzschild-radii and DistanceCoinc is the radial distance from the centre (black hole) to the point, at which both expressions for the flux-density coincide, in units of Schwarzschild-radii.'''
    return Toy1MagneticFluxDensityDipolarDilution(M9,Dotm,DistanceCoinc,Toy1Distance,0)

def Toy1MagneticFluxDensityExpansionDilution(M9,Dotm,DistanceCoinc,DistanceGap,thetaGap,SizeGap,thetaOpening,Toy1Distance,pB):
    '''This is an approximation for the magnetic flux-density according to my estimation, which assumes a dipolar field near the gap and a decline which is direct proportional to Toy1SizeOfBlob^-pB. It is in units of Gauss. thetaGap is the polar angle of the location of the gap, DistanceGap is the distance between the gap and the black hole in units of Schwarzschild-radii and DistanceCoinc is the radial distance from the centre (black hole) to the point, at which both expressions for the flux-density coincide, in units of Schwarzschild-radii. SizeGap is a typical length-scale of the gap in units of Schwarzschild-radii. thetaOpening is the half opening angle of the jet. Toy1Distance is the distance of the blob to the centre in units of Schwarzschild-radii and pB is the exponent parametrising the decline of the magnetic flux-density.'''
    return Toy1MagneticFluxDensityDipolarDilution(M9,Dotm,DistanceCoinc,DistanceGap,thetaGap)*(SizeGap/Toy1SizeOfBlob(Toy1Distance,DistanceGap,SizeGap,thetaOpening))**pB

def MagneticFluxDensityBHGap(betaMagnetic,M9):
    '''This is equation 22 by 2017ApJ...841...61A Aharonian. It gives the magnetic flux-density due to an ADAF near a black hole with the assumption of an unscreened gap in units of T. betaMagnetic is the disk magnetisation.'''
    return 0.19*(betaMagnetic/M9)**(4/7)

## Electric field-strength and voltage

HirotaniElectricFieldStrengthMax = 1.35*10**2 # The maximum electric field-strength in units of statV/cm. According to the vacuum case with a gap width of about 3/8 . Cf. page 5 of 2016ApJ...818...50H Hirotani, Pu.
HirotaniElectricFieldStrengthMaxSI = HirotaniElectricFieldStrengthMax*(c/(10**6))*100 # The maximum electric field-strength in units of V/m.

def LRVoltage(B4,M9,hGap):
    '''The voltage along the gap in units of V according to equation 10 of 2011ApJ...730..123L Levinson, Rieger. hGap is the gap height in Schwarzschild-radii.'''
    return 1.7*10**21*B4*M9*(hGap)**2
    
def LRElectricFieldStrength(B4,M9,hGap):
    '''This gives the electric field-strength in the gap of length hGap (in units of Schwarzschild-radii), according to 2011ApJ...730..123L Levinson, Rieger in units of V/m'''
    return LRVoltage(B4,M9,hGap)/(hGap*rSchwarzschild(M9))

def CompareElectricFieldStrengths():
    print('Maximum electric field-strength of Hirotani:', HirotaniElectricFieldStrengthMaxSI, 'V/m')
    print('Electric field-strength of L&R:', LRElectricFieldStrength(1,0.3,3/8), 'V/m for B4=1, M9=0.3 and hGap=3/8')

def Toy1ElectricFieldStrength(Toy1ElectricFieldStrengthMax,Toy1Distance,DistanceGap,hGap):
    '''This gives the field-strength of an electric field in a gap. Toy1Distance and DistanceGap are the distance to the centre and the distance of the gap to the centre, respectively, in units of Schwarzschild-radii. hGap is the gap height in Schwarzschild-radii. Toy1ElectricFieldStrengthMax is the maximally achieved electric field-strength, that is the field-strength at Toy1Distance=DistanceGap. Both Toy1ElectricFieldStrength and Toy1ElectricFieldStrengthMax are in units of V/m. This function is motivated by 2016ApJ...818...50H Hirotani, Pu'''
    Result = np.where(abs(Toy1Distance-DistanceGap)>=hGap/2.0, 0.0, Toy1ElectricFieldStrengthMax*(1.0-(2.0*(Toy1Distance-DistanceGap)/hGap)**2))
    if np.size(Result)==1:
        return float(Result)
    else:
        return Result


## Equations concerning the BZ-power

def BZPowerHirotani(aStar,B4,M9):
    '''This is equation 1 by 2016ApJ...833..142H Hirotani et al. It is in units of W.'''
    return 10**38*aStar**2*M9*B4**2

def BZPowerHirotaniEquipartition(aStar,Dotm,M9):
    '''This is equation 3 by 2016ApJ...833..142H Hirotani et al, hence assuming an equipartition magnetic field. It is in units of W.'''
    return 1.7*10**39*aStar**2*Dotm*M9

def MaximumGapPowerAharonian(betaMagnetic,kappa,hGap,M9,PolarAngle):
    '''This is similar to equation 23 by 2017ApJ...841...61A Aharonian. It yields an estimate of the absolute maximum gap power after some assumptions, especially after assuming an open gap and a magnetic field of an ADAF. In comparison to eq. 23 the variability time scale was replaced by the gap height, i.e. t_{var,5}=1.7*M8*h/r_g. It is in units of W. betaMagnetic is the disk magnetisation, kappa is the multiplicity as defined by Aharonian, hGap is the height of the gap in units of Schwarzschild-radii. PolarAngle is about pi/2 according to Rieger 2011.'''
    return 4.8*10**37*betaMagnetic**(8/7)*kappa*hGap*M9**(6/7)*(np.sin(PolarAngle))**2

## Equations concerning the curvature-radius

def Toy1CurvatureRadiusOfDipolarFieldNearAxis(Toy1Distance,theta):
    '''This gives the curvature-radius of a dipole field according to my analysis (my diploma thesis equation 2.23 with rho=0 in the fraction and with the parametrisation 2.14 without alpha) in units of Schwarzschild-radii. It is only valid for theta approx 0. Toy1Distance is the distance to the black hole in units of Schwarzschild-radii and theta is the polar angle.'''
    return (4*Toy1Distance)/(3.0*np.sin(theta))

def Toy1CurvatureRadiusExpansionDilution(Toy1Distance,thetaOpening,SizeGap,DistanceGap,thetaGap,M9):
    '''This gives the curvature-radius that is stretched due to expansion according to my analysis in units m. It is only valid for thetaGap approx 0. Toy1Distance and DistanceGap is the distance and the distance of the gap to the black hole in units of Schwarzschild-radii and thetaGap is the typical polar angle of the gap. thetaOpening is the half opening angle of the jet. SizeGap is the size of the gap in units of Schwarzschild-radii.'''
    return Toy1CurvatureRadiusOfDipolarFieldNearAxis(rSchwarzschildTom(DistanceGap,M9),thetaGap)*(Toy1SizeOfBlob(Toy1Distance,DistanceGap,SizeGap,thetaOpening)/SizeGap)


## Dynamics of the Lorentz-factor (simplified and obsolete)

def Toy1gammaICKN(M9,T,B4,hGap):
    '''This determines the equilibrium gamma due to inverse-Compton-scattering in the Klein-Nishina-regime. B4 might be the magnetic flux-density in 10**4 Gauss, hGap is the gap height in Schwarzschild-radii and T is the effective temperature of the radiation field, that IC-scatters the electrons.'''
    return me*c**2/(4*kB*T)*np.exp(((h**3*e*2*1.7*10**21*B4*M9*hGap)/(np.pi**3*kB**2*T**2*rSchwarzschild(M9)*sigmaT*me**2*c))+5.0/6+0.58+0.57)

def Toy1gammaTPP(M9,T,B4,hGap):
    '''This determines the equilibrium gamma due to triplet pair-production. B4 might be the magnetic flux-density in 10**4 Gauss, hGap is the gap height in Schwarzschild-radii and T is the effective temperature of the radiation field.'''
    return ((h**3*e*3*1.7*10**21*B4*M9*hGap)/(24*rSchwarzschild(M9)*(1/137)*sigmaT*me**(3/2)*np.pi**(3/2)*1.34*kB**(5/2)*T**(5/2)))**2

def Toy1gammaSynchrotron(M9,B4,hGap):
    '''This determines the equilibrium gamma due to synchrotron-radiation. B4 might be the magnetic flux-density in 10**4 Gauss, hGap is the gap height in Schwarzschild-radii.'''
    return np.sqrt((9*np.pi*me**2*c**2*epsilon0*1.7*10**21*M9*hGap)/(rSchwarzschild(M9)*e**3*B4))

def Toy1gammaICCooled(M9,B4,hGap,uSoft,Toy1Distance,DistanceGap):
    '''This determines the gamma cooled via inverse-Compton-scattering in the Thomson-regime. B4 might be the magnetic flux-density in 10**4 Gauss, hGap is the gap height in Schwarzschild-radii, uSoft is the energy-density of the soft background photons in the gap in units of erg/cm**3. ublackbody is the energy-density of the soft background photons as they are emitted of the disk in units of erg/cm**3. Toy1Distance and DistanceGap are the distance to the centre and the distance of the gap to the centre, respectively, in units of Schwarzschild-radii.'''
    return 1/(1/LRgammaICOfEquation14(B4,hGap,uSoft)+(sigmaT*ergPercm3ToJPerm3(uSoft)*rSchwarzschildTom(DistanceGap,M9)/(me*c**2))*(1-rSchwarzschildTom(DistanceGap,M9)/rSchwarzschildTom(Toy1Distance,M9)))

def Toy1gammaCurvatureCooled(M9,B4,SizeGap,thetaOpening,hGap,Toy1Distance,DistanceGap,thetaGap):
    '''This determines the gamma cooled via curvature-radiation. B4 might be the magnetic flux-density in 10**4 Gauss, hGap is the gap height in Schwarzschild-radii. SizeGap is the size of the gap in units of Schwarzschild-radii. thetaOpening is the half opening angle of the jet. Toy1Distance and DistanceGap are the distance to the centre and the distance of the gap to the centre, respectively, in units of Schwarzschild-radii. thetaGap is the polar angle coordinate of the gap.'''
    return (1/(LRgammaCurOfEquation13(B4,M9,hGap,Toy1CurvatureRadiusOfDipolarFieldNearAxis(DistanceGap,thetaGap)))**3+(9*e**2*(np.sin(thetaGap))**2*rSchwarzschildTom(SizeGap,M9)/(64*np.pi*epsilon0*me*c**2*(rSchwarzschildTom(DistanceGap,M9))**2*np.tan(thetaOpening)))*(1-1/(1+2*(rSchwarzschildTom(Toy1Distance,M9)-rSchwarzschildTom(DistanceGap,M9))*np.tan(thetaOpening)/rSchwarzschildTom(SizeGap,M9))))**(-1/3)


## Curvature-radiation spectrum

def ECritOfCurvatureRadiation(M9,gamma,rCurvature):
    '''This determines the critical energy of curvature-radiation in units of J, cf. e.g. Longair's chapter 8.4.4. gamma is the Lorentz-factor of the radiating particle. rCurvature is the curvature-radius of the magnetic field-lines in units of Schwarzschild-radii.'''
    return (3*h*c*gamma**3)/(4*np.pi*rSchwarzschildTom(rCurvature,M9))


## Energy-loss-rates

Toy1IntegrandOfEnergyLossRateTPPAbbreviation1 = 218.0/84.0 # This makes evaluation of the following functions faster.
def Toy1IntegrandOfEnergyLossRateTPP(epsilon,gamma,SpectralNumberDensity,*PositionalArgumentsOfSpectralNumberDensityExceptepsilon):
    '''This gives the integrand for the triplet pair-production energy-loss-rate (especially for the function EnergyLossRateTPP). It is in units of 1/(m^3*J^2). epsilon is the photons' energy in units of J, gamma is the Lorentz-factor.
    SpectralNumberDensity can be any spectral number-density (spectral in energy), it just has to take epsilon as first argument. Then the following arguments can be delivered via *PositionalArgumentsOfSpectralNumberDensityExceptepsilon.
    In case SpectralNumberDensity is given ADAFNumberDensitySpectralInEnergy, then *PositionalArgumentsOfSpectralNumberDensityExceptepsilon is a place holder for the arguments SpectralLuminosity, alphaViscosity, betaPressure, ADAFrInner, ADAFrOuter, ADAFTemperatureElectron, M9, Dotm.'''
    Toy1IntegrandOfEnergyLossRateTPPAbbreviation2 = np.log(2.0*epsilon*gamma/(me*c**2))
    return SpectralNumberDensity(epsilon,*PositionalArgumentsOfSpectralNumberDensityExceptepsilon)*(Toy1IntegrandOfEnergyLossRateTPPAbbreviation2)**2*(Toy1IntegrandOfEnergyLossRateTPPAbbreviation2-Toy1IntegrandOfEnergyLossRateTPPAbbreviation1)/epsilon

def Toy1ADAFEnergyLossRateTPP(Toy1Distance,gamma,DistanceGap,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,PrintIntermediateResults=False):
    '''This determines the absolute value of the energy-loss-rate due to triplet pair-production for an isotropically diluted distribution of ADAF-photons according to 1997ApJ...477..585M. It is in units of J/s. gamma is the Lorentz-factor of the particle, it mustn't be an array. Toy1Distance and DistanceGap are the distance to the centre and the distance of the gap to the centre, respectively, in units of Schwarzschild-radii.'''
    PartInFrontOfIntegral = 3*0.386*(1.0/137.0)*sigmaT*me**2*c**5*(rSchwarzschildTom(DistanceGap,M9))**2/(8*(rSchwarzschildTom(Toy1Distance,M9))**2)
    EnergyTPPThreshold = 7*me*c**2/gamma # The triplet pair-production threshold.
    EnergyBremsstrahlungCutOff = 300*kB*ADAFTemperatureElectron # This may be a delicate value. It is an approximate cut-off of the bremsstrahlung-distribution.
    EnergyComptonisedCutOff = 3*kB*ADAFTemperatureElectron # This is the upper border of the reign of the comptonised contribution.
    EnergynuPeak = ApplyADAFIntelligentStorage('ADAFnuPeak')*h # This is the border between the reign of the cyclosynchrotron- and the comptonised part.
    EnergynuMin = ApplyADAFIntelligentStorage('ADAFnuMin')*h # This is the kink in the cyclosynchrotron-part.
    # Bremsstrahlung-contribution:
    if EnergyBremsstrahlungCutOff<=EnergyTPPThreshold:
        BremsstrahlungContribution = 0
    else: # The case EnergyTPPThreshold<EnergyBremsstrahlungCutOff:
        BremsstrahlungContributionLowestBorder = EnergyTPPThreshold
        BremsstrahlungContributionHighestBorder = EnergyBremsstrahlungCutOff
        BremsstrahlungContributionNumberOfBorders = max(4,np.round(np.log10(BremsstrahlungContributionHighestBorder)-np.log10(BremsstrahlungContributionLowestBorder))) # However, this range has to be subdivided into smaller ranges, for the integration not to fail. Make at least 4 borders, that is 3 intervals.  
        BremsstrahlungContributionListOfBorders = np.logspace(np.log10(BremsstrahlungContributionLowestBorder),np.log10(BremsstrahlungContributionHighestBorder),BremsstrahlungContributionNumberOfBorders)
        BremsstrahlungContribution = 0
        for LeftBorder, RightBorder in zip(BremsstrahlungContributionListOfBorders[:-1], BremsstrahlungContributionListOfBorders[1:]):
            BremsstrahlungContribution += integrate.quad(Toy1IntegrandOfEnergyLossRateTPP, LeftBorder, RightBorder, args=(gamma,ADAFNumberDensitySpectralInEnergy,ADAFSpectralLuminosityBremsstrahlungOf30,alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0] # This is the numerical integration of the TPP-integrand for the bremsstrahlung-case. The upper integration border should actually be infinity. However, when it is too high a number, then the integration can fail. Therefore the upper border was reduced to the value EnergyBremsstrahlungCutOff.
    # Comptonised and cyclosynchrotron-contribution:
    if EnergyComptonisedCutOff<=EnergyTPPThreshold:
        CyclosynchrotronContributionBelownuMin = 0
        CyclosynchrotronContributionAbovenuMin = 0
        ComptonisedContribution = 0
    else: # The case EnergyTPPThreshold<EnergyComptonisedCutOff:
        if EnergynuPeak<=EnergyTPPThreshold:
            CyclosynchrotronContributionBelownuMin = 0
            CyclosynchrotronContributionAbovenuMin = 0
            ComptonisedLowestBorder = EnergyTPPThreshold # For the numerical integration one has to integrate from EnergyTPPThreshold up to EnergyComptonisedCutOff.
        else: # The case EnergyTPPThreshold<EnergynuPeak:
            ComptonisedLowestBorder = EnergynuPeak
            if EnergynuMin<=EnergyTPPThreshold:
                CyclosynchrotronContributionBelownuMin = 0
                CyclosynchrotronContributionAbovenuMinLowestBorder = EnergyTPPThreshold
            else: # The case EnergyTPPThreshold<EnergynuMin:
                CyclosynchrotronContributionAbovenuMinLowestBorder = EnergynuMin
                CyclosynchrotronContributionBelownuMinLowestBorder = EnergyTPPThreshold
                CyclosynchrotronContributionBelownuMinHighestBorder = EnergynuMin
                CyclosynchrotronContributionBelownuMinNumberOfBorders = max(4,np.round(np.log10(CyclosynchrotronContributionBelownuMinHighestBorder)-np.log10(CyclosynchrotronContributionBelownuMinLowestBorder))) # However, this range has to be subdivided into smaller ranges, for the integration not to fail. Make at least 4 borders, that is 3 intervals.  
                CyclosynchrotronContributionBelownuMinListOfBorders = np.logspace(np.log10(CyclosynchrotronContributionBelownuMinLowestBorder),np.log10(CyclosynchrotronContributionBelownuMinHighestBorder),CyclosynchrotronContributionBelownuMinNumberOfBorders)
                CyclosynchrotronContributionBelownuMin = 0
                for LeftBorder, RightBorder in zip(CyclosynchrotronContributionBelownuMinListOfBorders[:-1], CyclosynchrotronContributionBelownuMinListOfBorders[1:]):
                    CyclosynchrotronContributionBelownuMin += integrate.quad(Toy1IntegrandOfEnergyLossRateTPP, LeftBorder, RightBorder, args=(gamma,ADAFNumberDensitySpectralInEnergy,ADAFSpectralLuminosityCyclosynchrotronBelownuMin,alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0]
            CyclosynchrotronContributionAbovenuMinHighestBorder = EnergynuPeak
            CyclosynchrotronContributionAbovenuMinNumberOfBorders = max(4,np.round(np.log10(CyclosynchrotronContributionAbovenuMinHighestBorder)-np.log10(CyclosynchrotronContributionAbovenuMinLowestBorder))) # However, this range has to be subdivided into smaller ranges, for the integration not to fail. Make at least 4 borders, that is 3 intervals.  
            CyclosynchrotronContributionAbovenuMinListOfBorders = np.logspace(np.log10(CyclosynchrotronContributionAbovenuMinLowestBorder),np.log10(CyclosynchrotronContributionAbovenuMinHighestBorder),CyclosynchrotronContributionAbovenuMinNumberOfBorders)
            CyclosynchrotronContributionAbovenuMin = 0
            for LeftBorder, RightBorder in zip(CyclosynchrotronContributionAbovenuMinListOfBorders[:-1], CyclosynchrotronContributionAbovenuMinListOfBorders[1:]):
                CyclosynchrotronContributionAbovenuMin += integrate.quad(Toy1IntegrandOfEnergyLossRateTPP, LeftBorder, RightBorder, args=(gamma,ADAFNumberDensitySpectralInEnergy,ADAFSpectralLuminosityCyclosynchrotronOf25,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm))[0] # Numerical integration of ADAFSpectralLuminosityCyclosynchrotron above nu_min.
        ComptonisedHighestBorder = EnergyComptonisedCutOff
        ComptonisedNumberOfBorders = max(8,np.round(np.log10(ComptonisedHighestBorder)-np.log10(ComptonisedLowestBorder))) # However, this range has to be subdivided into smaller ranges, for the integration not to fail. Make at least 8 borders, that is 3 intervals.  
        ComptonisedListOfBorders = np.logspace(np.log10(ComptonisedLowestBorder),np.log10(ComptonisedHighestBorder),ComptonisedNumberOfBorders)
        ComptonisedContribution = 0
        for LeftBorder, RightBorder in zip(ComptonisedListOfBorders[:-1], ComptonisedListOfBorders[1:]):
            ComptonisedContribution += integrate.quad(Toy1IntegrandOfEnergyLossRateTPP, LeftBorder, RightBorder, args=(gamma,ADAFNumberDensitySpectralInEnergy,ADAFSpectralLuminosityComptonised,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0] # Numerical integration of ADAFSpectralLuminosityComptonised for each part-range.
    if PrintIntermediateResults==True:
        print('BremsstrahlungContribution',BremsstrahlungContribution)
        print('ComptonisedContribution',ComptonisedContribution)
        print('CyclosynchrotronContributionAbovenuMin',CyclosynchrotronContributionAbovenuMin)
        print('CyclosynchrotronContributionBelownuMin',CyclosynchrotronContributionBelownuMin)
    IntegralValue = BremsstrahlungContribution+ComptonisedContribution+CyclosynchrotronContributionAbovenuMin+CyclosynchrotronContributionBelownuMin
    Value = PartInFrontOfIntegral*IntegralValue
    return Value

def Toy1EnergyLossRateCurvature(Toy1Distance,gamma,DistanceGap,SizeGap,thetaGap,thetaOpening,M9):
    '''This determines the absolute value of the energy-loss-rate due to curvature-radiation in units of J/s. gamma is the Lorentz-factor of the particle. SizeGap is the size of the gap in units of Schwarzschild-radii. thetaOpening is the half opening angle of the jet. Toy1Distance and DistanceGap are the distance to the centre and the distance of the gap to the centre, respectively, in units of Schwarzschild-radii. thetaGap is the polar angle coordinate of the gap.'''
    return (2.0*e**2*c*gamma**4)/(12.0*np.pi*epsilon0*(Toy1CurvatureRadiusExpansionDilution(Toy1Distance,thetaOpening,SizeGap,DistanceGap,thetaGap,M9))**2)

def EnergyLossRateICThomson(gamma,EnergyDensity):
    '''This is the absolute value of the energy-loss-rate of IC-scattering in the Thomson-limit, i.e. equation 2.18 by B&G1970. It is in units of J/s. EnergyDensity is in units of J/m^3.'''
    return 4*gamma**2*c*sigmaT*EnergyDensity/3

Toy1IntegrandOfEnergyLossRateICKNAbbreviation = 11.0/6.0 # This makes evaluation of the following functions faster.
def Toy1IntegrandOfEnergyLossRateICKN(epsilon,gamma,SpectralNumberDensity,*PositionalArgumentsOfSpectralNumberDensityExceptepsilon):
    '''This gives the integrand for the inverse-Compton-Klein-Nishina-energy-loss-rate (especially for the function EnergyLossRateIC), i.e. this is the integrand of equation 2.57 by B&G1970. It is in units of 1/(m^3*J^2). epsilon is the photons' energy in units of J, gamma is the Lorentz-factor.
    SpectralNumberDensity can be any spectral number-density (spectral in energy), it just has to take epsilon as first argument. Then the following arguments can be delivered via *PositionalArgumentsOfSpectralNumberDensityExceptepsilon.
    In case SpectralNumberDensity is given ADAFNumberDensitySpectralInEnergy, then *PositionalArgumentsOfSpectralNumberDensityExceptepsilon is a place holder for the arguments SpectralLuminosity, alphaViscosity, betaPressure, ADAFrInner, ADAFrOuter, ADAFTemperatureElectron, M9, Dotm.'''
    return SpectralNumberDensity(epsilon,*PositionalArgumentsOfSpectralNumberDensityExceptepsilon)*(np.log(4*epsilon*gamma/(me*c**2))-Toy1IntegrandOfEnergyLossRateICKNAbbreviation)/epsilon

def IntegrandOfEnergyLossRateICKN(x,gamma,n0):
    '''This gives the integrand for the inverse-Compton-Klein-Nishina-energy-loss-rate (especially for the function EnergyLossRateIC), i.e. this is the integrand of equation 2.57 by B&G1970. It is in units of 1/(m^3). It is similar to Toy1IntegrandOfEnergyLossRateICKN, just for dimensionless energies and n0 being spectral in dimensionless energy.'''
    return n0(x)*(np.log(4*x*gamma)-Toy1IntegrandOfEnergyLossRateICKNAbbreviation)/x

# The energy-loss-rate in the extreme Klein-Nishina-limit is computed in part4 after the soft photon-field was defined.

def Toy1ADAFEnergyLossRateIC(Toy1Distance,gamma,DistanceGap,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,PrintIntermediateResults=False):
    '''This determines the absolute value of the energy-loss-rate due to inverse-Compton-scattering for an isotropically diluted ADAF-spectrum of background photons. It is in units of J/s. gamma is the Lorentz-factor of the particle. Toy1Distance and DistanceGap are the distance to the centre and the distance of the gap to the centre, respectively, in units of Schwarzschild-radii.'''
    ThPartInFrontOfIntegral = 4.0*sigmaT*c*gamma**2*(rSchwarzschildTom(DistanceGap,M9))**2/(3.0*(rSchwarzschildTom(Toy1Distance,M9))**2) # Term in front of the integral in the Thomson-energy-loss-rate. It is in units of m^3/s.
    KNPartInFrontOfIntegral = 3.0*sigmaT*me**2*c**5*(rSchwarzschildTom(DistanceGap,M9))**2/(8*(rSchwarzschildTom(Toy1Distance,M9))**2) # Term in front of the integral in the Klein-Nishina-energy-loss-rate. It is in units of J^2*m^3/s.
    EnergyLimit = 2*me*c**2/gamma # The limit between Thomson- and Klein-Nishina-regime.
    EnergyBremsstrahlungCutOff = 300*kB*ADAFTemperatureElectron # This may be a delicate value. It is an approximate cut-off of the bremsstrahlung-distribution.
    EnergyComptonisedCutOff = 3*kB*ADAFTemperatureElectron # This is the upper border of the reign of the comptonised contribution.
    EnergynuPeak = ApplyADAFIntelligentStorage('ADAFnuPeak')*h # This is the border between the reign of the cyclosynchrotron- and the comptonised part.
    EnergynuMin = ApplyADAFIntelligentStorage('ADAFnuMin')*h # This is the kink in the cyclosynchrotron-part.
    # Bremsstrahlung-contribution:
    ThBremsstrahlungLowestBorder = 0
    if EnergyBremsstrahlungCutOff<=EnergyLimit:
        KNBremsstrahlungContribution = 0
        ThBremsstrahlungHighestBorder = EnergyBremsstrahlungCutOff
    else:
        KNBremsstrahlungHighestBorder = EnergyBremsstrahlungCutOff
        KNBremsstrahlungLowestBorder = EnergyLimit
        ThBremsstrahlungHighestBorder = EnergyLimit
        # The numerical integration of the ICKN-integrand (which is in units of 1/(m^3*J^2)) for the bremsstrahlung-case:
        KNBremsstrahlungNumberOfBorders = max(4,np.round(np.log10(KNBremsstrahlungHighestBorder)-np.log10(KNBremsstrahlungLowestBorder)))
        KNBremsstrahlungListOfBorders = np.logspace(np.log10(KNBremsstrahlungLowestBorder),np.log10(KNBremsstrahlungHighestBorder),KNBremsstrahlungNumberOfBorders)
        KNBremsstrahlungContribution = 0 # In units of 1/(m^3*J).
        for LeftBorder, RightBorder in zip(KNBremsstrahlungListOfBorders[:-1], KNBremsstrahlungListOfBorders[1:]):
            KNBremsstrahlungContribution += integrate.quad(Toy1IntegrandOfEnergyLossRateICKN, LeftBorder, RightBorder, args=(gamma,ADAFNumberDensitySpectralInEnergy,ADAFSpectralLuminosityBremsstrahlungOf30,alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0]
    # The numerical integration of the IC-Thomson-integrand (which is in units of 1/m^3) for the bremsstrahlung-case:
    ThBremsstrahlungNumberOfBorders = 8
    ThBremsstrahlungListOfBorders = np.linspace(ThBremsstrahlungLowestBorder,ThBremsstrahlungHighestBorder,ThBremsstrahlungNumberOfBorders)
    ThBremsstrahlungContribution = 0 # In units of J/m^3.
    for LeftBorder, RightBorder in zip(ThBremsstrahlungListOfBorders[:-1], ThBremsstrahlungListOfBorders[1:]):
        ThBremsstrahlungContribution += integrate.quad(ADAFEnergyDensitySpectralInEnergy, LeftBorder, RightBorder, args=(ADAFSpectralLuminosityBremsstrahlungOf30,alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0]
    # The complete bremsstrahlung-contribution:
    BremsstrahlungContribution = ThPartInFrontOfIntegral*ThBremsstrahlungContribution+KNPartInFrontOfIntegral*KNBremsstrahlungContribution # In units of J/s.
    # Comptonised and cyclosynchrotron-contribution, Klein-Nishina-regime:
    if EnergyComptonisedCutOff<=EnergyLimit:
        KNCyclosynchrotronContributionBelownuMin = 0
        KNCyclosynchrotronContributionAbovenuMin = 0
        KNComptonisedContribution = 0
    else: # The case EnergyLimit<EnergyComptonisedCutOff:
        KNComptonisedHighestBorder = EnergyComptonisedCutOff
        if EnergynuPeak<=EnergyLimit:
            KNCyclosynchrotronContributionBelownuMin = 0
            KNCyclosynchrotronContributionAbovenuMin = 0
            KNComptonisedLowestBorder = EnergyLimit
        else: # The case EnergyLimit<EnergynuPeak:
            KNComptonisedLowestBorder = EnergynuPeak
            KNCyclosynchrotronAbovenuMinHighestBorder = EnergynuPeak
            if EnergynuMin<=EnergyLimit:
                KNCyclosynchrotronContributionBelownuMin = 0
                KNCyclosynchrotronAbovenuMinLowestBorder = EnergyLimit
            else: # The case EnergyLimit<EnergynuMin:
                KNCyclosynchrotronAbovenuMinLowestBorder = EnergynuMin
                KNCyclosynchrotronBelownuMinHighestBorder = EnergynuMin
                KNCyclosynchrotronBelownuMinLowestBorder = EnergyLimit
                # The numerical integration of the ICKN-integrand (which is in units of 1/(m^3*J^2)) for the cyclosynchrotron-case below nu_min:
                KNCyclosynchrotronBelownuMinNumberOfBorders = max(4,np.round(np.log10(KNCyclosynchrotronBelownuMinHighestBorder)-np.log10(KNCyclosynchrotronBelownuMinLowestBorder)))
                KNCyclosynchrotronBelownuMinListOfBorders = np.logspace(np.log10(KNCyclosynchrotronBelownuMinLowestBorder),np.log10(KNCyclosynchrotronBelownuMinHighestBorder),KNCyclosynchrotronBelownuMinNumberOfBorders)
                KNCyclosynchrotronContributionBelownuMin = 0 # In units of 1/(m^3*J).
                for LeftBorder, RightBorder in zip(KNCyclosynchrotronBelownuMinListOfBorders[:-1], KNCyclosynchrotronBelownuMinListOfBorders[1:]):
                    KNCyclosynchrotronContributionBelownuMin += integrate.quad(Toy1IntegrandOfEnergyLossRateICKN, LeftBorder, RightBorder, args=(gamma,ADAFNumberDensitySpectralInEnergy,ADAFSpectralLuminosityCyclosynchrotronBelownuMin,alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0]
            # The numerical integration of the ICKN-integrand (which is in units of 1/(m^3*J^2)) for the cyclosynchrotron-case above nu_min:
            KNCyclosynchrotronAbovenuMinNumberOfBorders = max(4,np.round(np.log10(KNCyclosynchrotronAbovenuMinHighestBorder)-np.log10(KNCyclosynchrotronAbovenuMinLowestBorder)))
            KNCyclosynchrotronAbovenuMinListOfBorders = np.logspace(np.log10(KNCyclosynchrotronAbovenuMinLowestBorder),np.log10(KNCyclosynchrotronAbovenuMinHighestBorder),KNCyclosynchrotronAbovenuMinNumberOfBorders)
            KNCyclosynchrotronContributionAbovenuMin = 0 # In units of 1/(m^3*J).
            for LeftBorder, RightBorder in zip(KNCyclosynchrotronAbovenuMinListOfBorders[:-1], KNCyclosynchrotronAbovenuMinListOfBorders[1:]):
                KNCyclosynchrotronContributionAbovenuMin += integrate.quad(Toy1IntegrandOfEnergyLossRateICKN, LeftBorder, RightBorder, args=(gamma,ADAFNumberDensitySpectralInEnergy,ADAFSpectralLuminosityCyclosynchrotronOf25,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm))[0]
        # The numerical integration of the ICKN-integrand (which is in units of 1/(m^3*J^2)) for the comptonised case:
        KNComptonisedNumberOfBorders = max(4,np.round(np.log10(KNComptonisedHighestBorder)-np.log10(KNComptonisedLowestBorder)))
        KNComptonisedListOfBorders = np.logspace(np.log10(KNComptonisedLowestBorder),np.log10(KNComptonisedHighestBorder),KNComptonisedNumberOfBorders)
        KNComptonisedContribution = 0 # In units of 1/(m^3*J).
        for LeftBorder, RightBorder in zip(KNComptonisedListOfBorders[:-1], KNComptonisedListOfBorders[1:]):
            KNComptonisedContribution += integrate.quad(Toy1IntegrandOfEnergyLossRateICKN, LeftBorder, RightBorder, args=(gamma,ADAFNumberDensitySpectralInEnergy,ADAFSpectralLuminosityComptonised,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0]
    # Comptonised and cyclosynchrotron-contribution, Thomson-regime:
    ThCyclosynchrotronBelownuMinLowestBorder = 0
    if EnergyLimit<=EnergynuMin:
        ThComptonisedContribution = 0
        ThCyclosynchrotronContributionAbovenuMin = 0
        ThCyclosynchrotronBelownuMinHighestBorder = EnergyLimit
    else: # The case EnergynuMin<EnergyLimit:
        ThCyclosynchrotronBelownuMinHighestBorder = EnergynuMin
        ThCyclosynchrotronAbovenuMinLowestBorder = EnergynuMin
        if EnergyLimit<=EnergynuPeak:
            ThComptonisedContribution = 0
            ThCyclosynchrotronAbovenuMinHighestBorder = EnergyLimit
        else: # The case EnergynuPeak<EnergyLimit:
            ThCyclosynchrotronAbovenuMinHighestBorder = EnergynuPeak
            ThComptonisedLowestBorder = EnergynuPeak
            if EnergyLimit<=EnergyComptonisedCutOff:
                ThComptonisedHighestBorder = EnergyLimit
            else: # The case EnergyComptonisedCutOff<EnergyLimit:
                ThComptonisedHighestBorder = EnergyComptonisedCutOff
            # The numerical integration of the IC-Thomson-integrand (which is in units of 1/m^3) for the comptonised case:
            ThComptonisedNumberOfBorders = max(4,np.round(np.log10(ThComptonisedHighestBorder)-np.log10(ThComptonisedLowestBorder)))
            ThComptonisedListOfBorders = np.logspace(np.log10(ThComptonisedLowestBorder),np.log10(ThComptonisedHighestBorder),ThComptonisedNumberOfBorders)
            ThComptonisedContribution = 0 # In units of J/m^3.
            for LeftBorder, RightBorder in zip(ThComptonisedListOfBorders[:-1], ThComptonisedListOfBorders[1:]):
                ThComptonisedContribution += integrate.quad(ADAFEnergyDensitySpectralInEnergy, LeftBorder, RightBorder, args=(ADAFSpectralLuminosityComptonised,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0]
        # The numerical integration of the IC-Thomson-integrand (which is in units of 1/m^3) for the cyclosynchrotron-case above nu_min:
        ThCyclosynchrotronAbovenuMinNumberOfBorders = max(4,np.round(np.log10(ThCyclosynchrotronAbovenuMinHighestBorder)-np.log10(ThCyclosynchrotronAbovenuMinLowestBorder)))
        ThCyclosynchrotronAbovenuMinListOfBorders = np.logspace(np.log10(ThCyclosynchrotronAbovenuMinLowestBorder),np.log10(ThCyclosynchrotronAbovenuMinHighestBorder),ThCyclosynchrotronAbovenuMinNumberOfBorders)
        ThCyclosynchrotronContributionAbovenuMin = 0 # In units of J/m^3.
        for LeftBorder, RightBorder in zip(ThCyclosynchrotronAbovenuMinListOfBorders[:-1], ThCyclosynchrotronAbovenuMinListOfBorders[1:]):
            ThCyclosynchrotronContributionAbovenuMin += integrate.quad(ADAFEnergyDensitySpectralInEnergy, LeftBorder, RightBorder, args=(ADAFSpectralLuminosityCyclosynchrotronOf25,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm))[0]
    # The numerical integration of the IC-Thomson-integrand (which is in units of 1/m^3) for the cyclosynchrotron-case below nu_min:
    ThCyclosynchrotronBelownuMinNumberOfBorders = 6
    ThCyclosynchrotronBelownuMinListOfBorders = np.linspace(ThCyclosynchrotronBelownuMinLowestBorder,ThCyclosynchrotronBelownuMinHighestBorder,ThCyclosynchrotronBelownuMinNumberOfBorders)
    ThCyclosynchrotronContributionBelownuMin = 0 # In units of J/m^3.
    for LeftBorder, RightBorder in zip(ThCyclosynchrotronBelownuMinListOfBorders[:-1], ThCyclosynchrotronBelownuMinListOfBorders[1:]):
        ThCyclosynchrotronContributionBelownuMin += integrate.quad(ADAFEnergyDensitySpectralInEnergy, LeftBorder, RightBorder, args=(ADAFSpectralLuminosityCyclosynchrotronBelownuMin,alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0]
    # The complete cyclosynchrotron-contribution below nu_min:
    CyclosynchrotronContributionBelownuMin = ThPartInFrontOfIntegral*ThCyclosynchrotronContributionBelownuMin+KNPartInFrontOfIntegral*KNCyclosynchrotronContributionBelownuMin # In units of J/s.
    # The complete cyclosynchrotron-contribution above nu_min:
    CyclosynchrotronContributionAbovenuMin = ThPartInFrontOfIntegral*ThCyclosynchrotronContributionAbovenuMin+KNPartInFrontOfIntegral*KNCyclosynchrotronContributionAbovenuMin # In units of J/s.
    # The complete comptonised contribution:
    ComptonisedContribution = ThPartInFrontOfIntegral*ThComptonisedContribution+KNPartInFrontOfIntegral*KNComptonisedContribution # In units of J/s.
    # The complete value:
    Value = CyclosynchrotronContributionBelownuMin+CyclosynchrotronContributionAbovenuMin+ComptonisedContribution+BremsstrahlungContribution # In units of J/s.
    if PrintIntermediateResults==True:
        print('Bremsstrahlung-part:')
        print('    KN integral value:',KNBremsstrahlungContribution,'/(m^3*J)')
        print('    Thomson integral value:',ThBremsstrahlungContribution,'J/m^3')
        print('    Complete energy-loss-rate:',BremsstrahlungContribution,'J/s')
        print('Comptonised part:')
        print('    KN integral value:',KNComptonisedContribution,'/(m^3*J)')
        print('    Thomson integral value:',ThComptonisedContribution,'J/m^3')
        print('    Complete energy-loss-rate:',ComptonisedContribution,'J/s')
        print('Cyclosynchrotron-part above nu_min:')
        print('    KN integral value:',KNCyclosynchrotronContributionAbovenuMin,'/(m^3*J)')
        print('    Thomson integral value:',ThCyclosynchrotronContributionAbovenuMin,'J/m^3')
        print('    Complete energy-loss-rate:',CyclosynchrotronContributionAbovenuMin,'J/s')
        print('Cyclosynchrotron-part below nu_min:')
        print('    KN integral value:',KNCyclosynchrotronContributionBelownuMin,'/(m^3*J)')
        print('    Thomson integral value:',ThCyclosynchrotronContributionBelownuMin,'J/m^3')
        print('    Complete energy-loss-rate:',CyclosynchrotronContributionBelownuMin,'J/s')
    return Value


## Dynamics of the Lorentz-factor

def Toy1dgammaOverddTPPADAF(Toy1Distance,gamma,DistanceGap,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm):
    '''This gives the derivative of gamma along the Toy1Distance due to triplet pair-production-energy-losses. It is in units of 1/m. Toy1Distance and DistanceGap are the distance to the centre and the distance of the gap to the centre, respectively, in units of Schwarzschild-radii.'''
    return  -Toy1ADAFEnergyLossRateTPP(Toy1Distance,gamma,DistanceGap,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)/(me*c**3)

def Toy1dgammaOverddCurvature(Toy1Distance,gamma,DistanceGap,SizeGap,thetaGap,thetaOpening,M9):
    '''This gives the derivative of gamma along the Toy1Distance due to curvature-radiation-energy-losses. It is in units of 1/m. SizeGap is the size of the gap in units of Schwarzschild-radii. thetaOpening is the half opening angle of the jet. Toy1Distance and DistanceGap are the distance to the centre and the distance of the gap to the centre, respectively, in units of Schwarzschild-radii. thetaGap is the polar angle coordinate of the gap.'''
    return -Toy1EnergyLossRateCurvature(Toy1Distance,gamma,DistanceGap,SizeGap,thetaGap,thetaOpening,M9)/(me*c**3)

def Toy1dgammaOverddICADAF(Toy1Distance,gamma,DistanceGap,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm):
    '''This gives the derivative of gamma along the Toy1Distance due to IC-energy-losses. It is in units of 1/m. Toy1Distance and DistanceGap are the distance to the centre and the distance of the gap to the centre, respectively, in units of Schwarzschild-radii.'''
    return -Toy1ADAFEnergyLossRateIC(Toy1Distance,gamma,DistanceGap,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)/(me*c**3)

def Toy1dgammaOverddAcceleration(Toy1ElectricFieldStrengthMax,Toy1Distance,DistanceGap,hGap):
    '''This gives the derivative of gamma along the Toy1Distance due to energy-gains due to electrostatic acceleration of an electron. It is in units of 1/m. Toy1Distance and DistanceGap are the distance to the centre and the distance of the gap to the centre, respectively, in units of Schwarzschild-radii. hGap is the gap height in Schwarzschild-radii. ElectricFieldmax is the maximally achieved electric field, that is the electric field at Toy1Distance=DistanceGap in units of V/m.'''
    return -e/(me*c**2)*Toy1ElectricFieldStrength(Toy1ElectricFieldStrengthMax,Toy1Distance,DistanceGap,hGap) 

def Toy1RightHandSideOfDEQForGamma(Toy1Distance,gamma,Toy1ElectricFieldStrengthMax,DistanceGap,hGap,SizeGap,thetaGap,thetaOpening,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm):
    '''This gives the derivative of gamma along the Toy1Distance in 1/m. The multiplication with rSchwarzschild is necessary because the left-hand side originally was in SI-units and is now in Schwarzschild-radii, too.'''
    return rSchwarzschild(M9)*(Toy1dgammaOverddAcceleration(Toy1ElectricFieldStrengthMax,Toy1Distance,DistanceGap,hGap)+Toy1dgammaOverddCurvature(Toy1Distance,gamma,DistanceGap,SizeGap,thetaGap,thetaOpening,M9)+Toy1dgammaOverddTPPADAF(Toy1Distance,gamma,DistanceGap,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)+Toy1dgammaOverddICADAF(Toy1Distance,gamma,DistanceGap,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))
    # Comment of v1_5: It was realised, that the summand Toy1dgammaOverddICADAF was mistakenly drawn out of the second factor of the product. this occurred sometimes after Estimates_3_(new_approach_for_determining_gamma).py. It is corrected now.

def Toy1SolvedgammaOverdd(Toy1ElectricFieldStrengthMax,DistanceGap,hGap,SizeGap,thetaGap,thetaOpening,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm):
    '''This procedure solves the DEQ and thus determines gamma(Toy1Distance). t corresponds to Toy1Distance and y corresponds to gamma. t0 is the starting Toy1Distance and t1 is the Toy1Distance up to which the solution is determined. y0 is the initial value of gamma, that is gamma(t0). dt is the Toy1Distance-step. t0, t1 and dt is measured in units of Schwarzschild-radii.'''
    y0, t0, t1, dt = 10, DistanceGap-hGap/2, 10000, 0.02*hGap
    InitialStepSize = dt
    tSolution, ySolution = np.asarray([]), np.asarray([])
    MyDEQ = integrate.ode(Toy1RightHandSideOfDEQForGamma).set_integrator('dopri5')
    MyDEQ.set_initial_value(y0, t0)
    MyDEQ.set_f_params(Toy1ElectricFieldStrengthMax,DistanceGap,hGap,SizeGap,thetaGap,thetaOpening,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)
    StepCounter = 1
    while MyDEQ.successful() and MyDEQ.t <= t1:
        print(StepCounter)
        if np.mod(StepCounter,np.round(1.5*hGap/InitialStepSize))==0:
            dt *= 3
            print('t is', MyDEQ.t)
            print('dt is now', dt)
        tSolution = np.append(tSolution, MyDEQ.t)
        ySolution = np.append(ySolution, MyDEQ.y[0])
        MyDEQ.integrate(MyDEQ.t+dt)
        StepCounter +=1
    return tSolution, ySolution, DistanceGap, hGap


## Evaluations

def _EvaluateLRMagneticFluxDensityOfEquation3():
    rRange=np.linspace(3,300,20000)
    ValuesForLRMagneticFluxDensityOfEquation31=LRMagneticFluxDensityOfEquation3(M9,0.0002,rRange)
    ValuesForLRMagneticFluxDensityOfEquation32=LRMagneticFluxDensityOfEquation3(M9,0.0001,rRange)
    pl.figure(figsize=(12, 9), num="Magnetic flux-density versus radius")
    pl.title('Figure 2: Magnetic flux-density versus radius', fontsize=24)
    pl.loglog(rRange,ValuesForLRMagneticFluxDensityOfEquation31, label='$\dot m = 0.0002$, $M_9 = 0.3$')
    pl.loglog(rRange,ValuesForLRMagneticFluxDensityOfEquation32, label='$\dot m = 0.0001$, $M_9 = 0.3$')
    pl.legend(loc="best", fontsize=20)
    pl.xlim(rRange[0],rRange[-1])
    pl.xlabel('radius $r / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('magnetic flux-density $B_{\mathrm{L&R}} / \mathrm{G}$', fontsize=24)
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def _EvaluateToy1MagneticFluxDensityDipolarDilution():
    Toy1DistanceRange=np.linspace(1,20,20000)
    ValuesForToy1MagneticFluxDensityDipolarDilution1=Toy1MagneticFluxDensityDipolarDilution(M9,0.0001,3,Toy1DistanceRange,0.1)
    ValuesForToy1MagneticFluxDensityDipolarDilution2=Toy1MagneticFluxDensityDipolarDilution(M9,0.0001,6,Toy1DistanceRange,0.1)
    ValuesForToy1MagneticFluxDensityDipolarDilution3=Toy1MagneticFluxDensityDipolarDilution(M9,0.0001,3,Toy1DistanceRange,0.7)
    ValuesForToy1MagneticFluxDensityDipolarDilution4=Toy1MagneticFluxDensityDipolarDilution(M9,0.0001,6,Toy1DistanceRange,0.7)
    ValuesForToy1MagneticFluxDensityDipolarDilution5=Toy1MagneticFluxDensityDipolarDilution(M9,0.0002,3,Toy1DistanceRange,0.1)
    ValuesForToy1MagneticFluxDensityDipolarDilution6=Toy1MagneticFluxDensityDipolarDilution(M9,0.0002,6,Toy1DistanceRange,0.1)
    ValuesForToy1MagneticFluxDensityDipolarDilution7=Toy1MagneticFluxDensityDipolarDilution(M9,0.0002,3,Toy1DistanceRange,0.7)
    ValuesForToy1MagneticFluxDensityDipolarDilution8=Toy1MagneticFluxDensityDipolarDilution(M9,0.0002,6,Toy1DistanceRange,0.7)
    pl.figure(figsize=(12, 9), num="Magnetic flux-density versus distance near vacuum gap")
    pl.title('Figure 3: Magnetic flux-density versus\ndistance near vacuum gap', fontsize=24)
    pl.plot(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilution1, label='$\dot m = 0.0001$, $d_{\mathrm{coinc}} = 3 r_{\mathrm{S}}$, $M_9 = 0.3$, $\\vartheta = 0.1$')
    pl.plot(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilution2, label='$\dot m = 0.0001$, $d_{\mathrm{coinc}} = 6 r_{\mathrm{S}}$, $M_9 = 0.3$, $\\vartheta = 0.1$')
    pl.plot(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilution3, label='$\dot m = 0.0001$, $d_{\mathrm{coinc}} = 3 r_{\mathrm{S}}$, $M_9 = 0.3$, $\\vartheta = 0.7$')
    pl.plot(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilution4, label='$\dot m = 0.0001$, $d_{\mathrm{coinc}} = 6 r_{\mathrm{S}}$, $M_9 = 0.3$, $\\vartheta = 0.7$')
    pl.plot(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilution5, label='$\dot m = 0.0002$, $d_{\mathrm{coinc}} = 3 r_{\mathrm{S}}$, $M_9 = 0.3$, $\\vartheta = 0.1$')
    pl.plot(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilution6, label='$\dot m = 0.0002$, $d_{\mathrm{coinc}} = 6 r_{\mathrm{S}}$, $M_9 = 0.3$, $\\vartheta = 0.1$')
    pl.plot(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilution7, label='$\dot m = 0.0002$, $d_{\mathrm{coinc}} = 3 r_{\mathrm{S}}$, $M_9 = 0.3$, $\\vartheta = 0.7$', linestyle='--', linewidth=1.5)
    pl.plot(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilution8, label='$\dot m = 0.0002$, $d_{\mathrm{coinc}} = 6 r_{\mathrm{S}}$, $M_9 = 0.3$, $\\vartheta = 0.7$', linestyle=':', linewidth=2)
    pl.yscale('log')
    pl.legend(loc="best", fontsize=18)
    pl.xlim(Toy1DistanceRange[0],Toy1DistanceRange[-1])
    pl.xlabel('distance $d / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('magnetic flux-density $B_{\mathrm{dipole}} / \mathrm{G}$', fontsize=24)
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def _EvaluateToy1MagneticFluxDensityDipolarDilutionOnAxis():
    Toy1DistanceRange=np.logspace(2,4,50)
    ValuesForToy1MagneticFluxDensityDipolarDilutionOnAxis1=Toy1MagneticFluxDensityDipolarDilutionOnAxis(M9,0.0002,3,Toy1DistanceRange)
    ValuesForToy1MagneticFluxDensityDipolarDilutionOnAxis2=Toy1MagneticFluxDensityDipolarDilutionOnAxis(M9,0.0001,3,Toy1DistanceRange)
    ValuesForToy1MagneticFluxDensityDipolarDilutionOnAxis3=Toy1MagneticFluxDensityDipolarDilutionOnAxis(M9,0.0002,6,Toy1DistanceRange)
    ValuesForToy1MagneticFluxDensityDipolarDilutionOnAxis4=Toy1MagneticFluxDensityDipolarDilutionOnAxis(M9,0.0001,6,Toy1DistanceRange)
    pl.figure(figsize=(12, 9), num="Magnetic flux-density versus distance on axis")
    pl.title('Figure 4: Magnetic flux-density\nversus distance on axis', fontsize=24)
    pl.loglog(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilutionOnAxis1, label='$\dot m = 0.0002$, $d_{\mathrm{coinc}} = 3 r_{\mathrm{S}}$, $M_9 = 0.3$')
    pl.loglog(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilutionOnAxis2, label='$\dot m = 0.0001$, $d_{\mathrm{coinc}} = 3 r_{\mathrm{S}}$, $M_9 = 0.3$')
    pl.loglog(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilutionOnAxis3, label='$\dot m = 0.0002$, $d_{\mathrm{coinc}} = 6 r_{\mathrm{S}}$, $M_9 = 0.3$')
    pl.loglog(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityDipolarDilutionOnAxis4, label='$\dot m = 0.0001$, $d_{\mathrm{coinc}} = 6 r_{\mathrm{S}}$, $M_9 = 0.3$')
    pl.legend(loc="best", fontsize=20)
    pl.xlim(Toy1DistanceRange[0],Toy1DistanceRange[-1])
    pl.xlabel('distance $d / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('magnetic flux-density $B_{\mathrm{dipole}} / \mathrm{G}$', fontsize=24)
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def _EvaluateToy1MagneticFluxDensityExpansionDilution():
    Toy1DistanceRange=np.logspace(2,4,20)
    ValuesForToy1MagneticFluxDensityExpansionDilution1=Toy1MagneticFluxDensityExpansionDilution(M9,0.0001,3,15,0.04,0.1,0.03,Toy1DistanceRange,1)
    ValuesForToy1MagneticFluxDensityExpansionDilution2=Toy1MagneticFluxDensityExpansionDilution(M9,0.0002,6,3,0.01,10,0.01,Toy1DistanceRange,1)
    ValuesForToy1MagneticFluxDensityExpansionDilution4=Toy1MagneticFluxDensityExpansionDilution(M9,0.0001,3,5,0.02,0.1,0.02,Toy1DistanceRange,1)
    ValuesForToy1MagneticFluxDensityExpansionDilution5=Toy1MagneticFluxDensityExpansionDilution(M9,0.0002,6,5,0.02,0.1,0.02,Toy1DistanceRange,1)
    pl.figure(figsize=(12, 9), num="Magnetic flux-density in blob versus distance")
    pl.title('Figure 6b: Magnetic flux-density in blob versus distance', fontsize=24)
    pl.loglog(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityExpansionDilution4, label='$M_9 = 0.3$, \n $\dot m = 0.0001$, \n $d_{\mathrm{coinc}} = 3 r_{\mathrm{S}}$, \n $d_{\mathrm{gap}} = 5 r_{\mathrm{S}}$, \n $\\vartheta_{\mathrm{gap}} = 0.02$, \n $s_1 = 0.1 r_{\mathrm{S}}$, \n $\\theta = 0.02$, \n $p_{\mathrm{B}} = 1$ \n')
    pl.loglog(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityExpansionDilution5, label='$M_9 = 0.3$, \n $\dot m = 0.0002$, \n $d_{\mathrm{coinc}} = 6 r_{\mathrm{S}}$, \n $d_{\mathrm{gap}} = 5 r_{\mathrm{S}}$, \n $\\vartheta_{\mathrm{gap}} = 0.02$, \n $s_1 = 0.1 r_{\mathrm{S}}$, \n $\\theta = 0.02$, \n $p_{\mathrm{B}} = 1$')
    pl.fill_between(Toy1DistanceRange,ValuesForToy1MagneticFluxDensityExpansionDilution1,ValuesForToy1MagneticFluxDensityExpansionDilution2,color='blue', alpha=0.1)
    pl.legend(loc="best", fontsize=18)
    pl.xlim(Toy1DistanceRange[0],Toy1DistanceRange[-1])
    pl.xlabel('distance $d / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('magnetic flux-density $B_{\mathrm{blob}} / \mathrm{G}$', fontsize=24)
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def _EvaluateToy1CurvatureRadiusOfDipolarFieldNearAxis():
    Toy1DistanceRange=np.linspace(1,20,20000)
    ValuesForToy1CurvatureRadiusOfDipolarFieldNearAxis1=Toy1CurvatureRadiusOfDipolarFieldNearAxis(Toy1DistanceRange,0.1)
    ValuesForToy1CurvatureRadiusOfDipolarFieldNearAxis2=Toy1CurvatureRadiusOfDipolarFieldNearAxis(Toy1DistanceRange,0.3)
    ValuesForToy1CurvatureRadiusOfDipolarFieldNearAxis3=Toy1CurvatureRadiusOfDipolarFieldNearAxis(Toy1DistanceRange,0.5)
    ValuesForToy1CurvatureRadiusOfDipolarFieldNearAxis4=Toy1CurvatureRadiusOfDipolarFieldNearAxis(Toy1DistanceRange,0.7)
    pl.figure(figsize=(12, 9), num="Curvature-radius versus distance near vacuum gap")
    pl.title('Figure 5: Curvature-radius versus distance near vacuum gap', fontsize=24)
    pl.plot(Toy1DistanceRange,ValuesForToy1CurvatureRadiusOfDipolarFieldNearAxis1, label='$\\vartheta = 0.1$')
    pl.plot(Toy1DistanceRange,ValuesForToy1CurvatureRadiusOfDipolarFieldNearAxis2, label='$\\vartheta = 0.3$')
    pl.plot(Toy1DistanceRange,ValuesForToy1CurvatureRadiusOfDipolarFieldNearAxis3, label='$\\vartheta = 0.5$')
    pl.plot(Toy1DistanceRange,ValuesForToy1CurvatureRadiusOfDipolarFieldNearAxis4, label='$\\vartheta = 0.7$')
    pl.plot(Toy1DistanceRange,np.ones(len(Toy1DistanceRange))*1000, linewidth=20,color='y', alpha=0.2, label='Approximate situation of gap')
    pl.legend(loc="best", fontsize=20)
    pl.xlim(Toy1DistanceRange[0],Toy1DistanceRange[-1])
    pl.ylim(0,300)
    pl.xlabel('distance $d / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('curvature-radius $r_{\mathrm{curv}} / r_\mathrm{S}$', fontsize=24)
    ax = pl.gca()
    ax.add_patch(pl.Polygon([[3,10],[15,100],[15,1000],[3,1000]], closed=True, facecolor='y', alpha=0.2, edgecolor='none'))
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def _EvaluateToy1RadiusOfBlob():
    Toy1DistanceRange=np.logspace(np.log10(5),4,50)
    ValuesForToy1RadiusOfBlob1=Toy1RadiusOfBlob(0.1,Toy1DistanceRange,5,0.01,M9)
    ValuesForToy1RadiusOfBlob2=Toy1RadiusOfBlob(0.1,Toy1DistanceRange,5,0.03,M9)
    ValuesForToy1RadiusOfBlob3=Toy1RadiusOfBlob(1,Toy1DistanceRange,5,0.01,M9)
    ValuesForToy1RadiusOfBlob4=Toy1RadiusOfBlob(1,Toy1DistanceRange,5,0.03,M9)
    ValuesForToy1RadiusOfBlob5=Toy1RadiusOfBlob(10,Toy1DistanceRange,5,0.01,M9)
    ValuesForToy1RadiusOfBlob6=Toy1RadiusOfBlob(10,Toy1DistanceRange,5,0.03,M9)
    pl.figure(figsize=(12, 9), num="radius of blob versus distance")
    pl.title('Figure 7: radius of blob versus distance', fontsize=24)
    pl.loglog(Toy1DistanceRange,ValuesForToy1RadiusOfBlob1, label='$d_{\mathrm{gap}} = 5 r_{\mathrm{S}}$, $M_9 = 0.3$, $s_1 = 0.1 r_{\mathrm{S}}$, $\\theta = 0.01$')
    pl.loglog(Toy1DistanceRange,ValuesForToy1RadiusOfBlob2, label='$d_{\mathrm{gap}} = 5 r_{\mathrm{S}}$, $M_9 = 0.3$, $s_1 = 0.1 r_{\mathrm{S}}$, $\\theta = 0.03$', linestyle='--', linewidth=2)
    pl.loglog(Toy1DistanceRange,ValuesForToy1RadiusOfBlob3, label='$d_{\mathrm{gap}} = 5 r_{\mathrm{S}}$, $M_9 = 0.3$, $s_1 = 1 r_{\mathrm{S}}$, $\\theta = 0.01$')
    pl.loglog(Toy1DistanceRange,ValuesForToy1RadiusOfBlob4, label='$d_{\mathrm{gap}} = 5 r_{\mathrm{S}}$, $M_9 = 0.3$, $s_1 = 1 r_{\mathrm{S}}$, $\\theta = 0.03$')
    pl.loglog(Toy1DistanceRange,ValuesForToy1RadiusOfBlob5, label='$d_{\mathrm{gap}} = 5 r_{\mathrm{S}}$, $M_9 = 0.3$, $s_1 = 10 r_{\mathrm{S}}$, $\\theta = 0.01$')
    pl.loglog(Toy1DistanceRange,ValuesForToy1RadiusOfBlob6, label='$d_{\mathrm{gap}} = 5 r_{\mathrm{S}}$, $M_9 = 0.3$, $s_1 = 10 r_{\mathrm{S}}$, $\\theta = 0.03$')
    pl.legend(loc="best", fontsize=18)
    pl.xlim(Toy1DistanceRange[0],Toy1DistanceRange[-1])
    pl.ylim(10**12,10**17)
    pl.xlabel('distance $d / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('blob-radius $R_\mathrm{blob} / \mathrm{cm}$', fontsize=24)
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def _EvaluateLRgammaCurOfEquation13():
    hGapRange=np.linspace(0.01,1,20)
    ValuesForLRgammaCurOfEquation131=LRgammaCurOfEquation13(1,M9,hGapRange,10)
    ValuesForLRgammaCurOfEquation132=LRgammaCurOfEquation13(1,M9,hGapRange,100)
    ValuesForLRgammaCurOfEquation133=LRgammaCurOfEquation13(0.1,M9,hGapRange,10)
    ValuesForLRgammaCurOfEquation134=LRgammaCurOfEquation13(0.1,M9,hGapRange,100)
    ValuesForLRgammaCurOfEquation135=LRgammaCurOfEquation13(0.01,M9,hGapRange,10)
    ValuesForLRgammaCurOfEquation136=LRgammaCurOfEquation13(0.01,M9,hGapRange,100)
    pl.figure(figsize=(12, 9), num="Lorentz-factor (curvature-radiation) versus gap height")
    pl.title('Figure 8: Lorentz-factor (curvature-radiation) versus gap height', fontsize=24)
    pl.loglog(hGapRange,ValuesForLRgammaCurOfEquation131, linestyle='--', linewidth=2, label='$B_4 = 1$, $r_{\mathrm{curv}} = 10 r_{\mathrm{S}}$, $M_9 = 0.3$')
    pl.loglog(hGapRange,ValuesForLRgammaCurOfEquation132, label='$B_4 = 1$, $r_{\mathrm{curv}} = 100 r_{\mathrm{S}}$, $M_9 = 0.3$')
    pl.loglog(hGapRange,ValuesForLRgammaCurOfEquation133, label='$B_4 = 0.1$, $r_{\mathrm{curv}} = 10 r_{\mathrm{S}}$, $M_9 = 0.3$')
    pl.loglog(hGapRange,ValuesForLRgammaCurOfEquation134, label='$B_4 = 0.1$, $r_{\mathrm{curv}} = 100 r_{\mathrm{S}}$, $M_9 = 0.3$')
    pl.loglog(hGapRange,ValuesForLRgammaCurOfEquation135, label='$B_4 = 0.01$, $r_{\mathrm{curv}} = 10 r_{\mathrm{S}}$, $M_9 = 0.3$')
    pl.loglog(hGapRange,ValuesForLRgammaCurOfEquation136, label='$B_4 = 0.01$, $r_{\mathrm{curv}} = 100 r_{\mathrm{S}}$, $M_9 = 0.3$')
    pl.legend(loc="best", fontsize=20)
    pl.xlim(hGapRange[0],hGapRange[-1])
    pl.xlabel('gap height $h / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('Lorentz-factor $\gamma_{\mathrm{cur}}$', fontsize=24)
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def _EvaluateLRgammaICOfEquation14():
    hGapRange=np.linspace(0.01,1,20)
    ValuesForLRgammaICOfEquation141=LRgammaICOfEquation14(1,hGapRange,50)
    ValuesForLRgammaICOfEquation142=LRgammaICOfEquation14(1,hGapRange,100)
    ValuesForLRgammaICOfEquation143=LRgammaICOfEquation14(0.1,hGapRange,50)
    ValuesForLRgammaICOfEquation144=LRgammaICOfEquation14(0.1,hGapRange,100)
    ValuesForLRgammaICOfEquation145=LRgammaICOfEquation14(0.01,hGapRange,50)
    ValuesForLRgammaICOfEquation146=LRgammaICOfEquation14(0.01,hGapRange,100)
    pl.figure(figsize=(12, 9), num="Lorentz-factor (IC-scattering in Thomson-regime) versus gap height")
    pl.title('Figure 9: Lorentz-factor (inverse-Compton-scattering in Thomson-regime) versus gap height', fontsize=24)
    pl.loglog(hGapRange,ValuesForLRgammaICOfEquation141, label='$B_4 = 1$, $u_{\mathrm{bb}} = 50 \mathrm{erg}/\mathrm{cm}^3$')
    pl.loglog(hGapRange,ValuesForLRgammaICOfEquation142, label='$B_4 = 1$, $u_{\mathrm{bb}} = 100 \mathrm{erg}/\mathrm{cm}^3$')
    pl.loglog(hGapRange,ValuesForLRgammaICOfEquation143, label='$B_4 = 0.1$, $u_{\mathrm{bb}} = 50 \mathrm{erg}/\mathrm{cm}^3$')
    pl.loglog(hGapRange,ValuesForLRgammaICOfEquation144, label='$B_4 = 0.1$, $u_{\mathrm{bb}} = 100 \mathrm{erg}/\mathrm{cm}^3$')
    pl.loglog(hGapRange,ValuesForLRgammaICOfEquation145, label='$B_4 = 0.01$, $u_{\mathrm{bb}} = 50 \mathrm{erg}/\mathrm{cm}^3$')
    pl.loglog(hGapRange,ValuesForLRgammaICOfEquation146, label='$B_4 = 0.01$, $u_{\mathrm{bb}} = 100 \mathrm{erg}/\mathrm{cm}^3$')
    pl.legend(loc="best", fontsize=20)
    pl.xlim(hGapRange[0],hGapRange[-1])
    pl.xlabel('gap height $h / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('Lorentz-factor $\gamma_{\mathrm{IC,Thomson}}$', fontsize=24)
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def _EvaluateToy1gammaICCooled():
    distanceRange1=np.linspace(2,700,20000)
    distanceRange2=np.linspace(15,700,20000)
    distanceRange3=np.linspace(5,700,20000)
    ValuesForToy1gammaICCooled1=Toy1gammaICCooled(M9,1,1,76,distanceRange1,2)
    ValuesForToy1gammaICCooled2=Toy1gammaICCooled(M9,0.01,0.01,76,distanceRange2,15)
    ValuesForToy1gammaICCooled3=Toy1gammaICCooled(M9,0.1,0.1,76,distanceRange3,5)
    pl.figure(figsize=(12, 9), num="Lorentz-factor (cooled via IC-scattering) versus distance")
    pl.title('Figure 10: Lorentz-factor (cooled via IC-scattering) versus distance', fontsize=24)
    pl.plot(distanceRange1,ValuesForToy1gammaICCooled1, label='Highest possible values')
    pl.plot(distanceRange2,ValuesForToy1gammaICCooled2, label='Lowest possible values')
    pl.plot(distanceRange3,ValuesForToy1gammaICCooled3, label='\n\n$M_9 = 0.3$, $B_4 = 0.1$, \n$u_{\mathrm{bb}} = 76 \mathrm{erg}/\mathrm{cm}^3$, \n$d_{\mathrm{gap}} = 5 r_{\mathrm{S}}$, $h_{\mathrm{gap}} = 0.1 r_{\mathrm{S}}$')
    pl.legend(loc="best", fontsize=20)
    pl.xlim(1,distanceRange1[-1])
    pl.xlabel('distance $d / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('Lorentz-factor $\gamma$', fontsize=24)
    pl.xscale('log')
    pl.yscale('log')
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def _EvaluateToy1gammaCurvatureCooled():
    distanceRange1=np.linspace(2,10000,200000)
    distanceRange2=np.linspace(15,10000,200000)
    distanceRange3=np.linspace(5,10000,200000)
    ValuesForToy1gammaCurvatureCooled1=Toy1gammaCurvatureCooled(M9,0.01,10,0.01,0.01,distanceRange1,2,0.04)
    ValuesForToy1gammaCurvatureCooled2=Toy1gammaCurvatureCooled(M9,1,0.1,0.03,1,distanceRange2,15,0.01)
    ValuesForToy1gammaCurvatureCooled3=Toy1gammaCurvatureCooled(M9,0.1,1,0.02,0.1,distanceRange3,5,0.025)
    pl.figure(figsize=(12, 9), num="Lorentz-factor (cooled via curvature-radiation) versus distance")
    pl.title('Figure 11: Lorentz-factor (cooled via curvature-radiation) versus distance', fontsize=24)
    pl.plot(distanceRange1,ValuesForToy1gammaCurvatureCooled1, label='Lowest possible values')
    pl.plot(distanceRange2,ValuesForToy1gammaCurvatureCooled2, label='Highest possible values')
    pl.plot(distanceRange3,ValuesForToy1gammaCurvatureCooled3, label='\n$M_9 = 0.3$, $B_4 = 0.1$, \n$s_1 = 1 r_{\mathrm{S}}$, $\\theta = 0.02$, $h_{\mathrm{gap}} = 0.1 r_{\mathrm{S}}$, \n$d_{\mathrm{gap}} = 5 r_{\mathrm{S}}$, $\\vartheta_{\mathrm{gap}} = 0.025$')
    pl.legend(loc="best", fontsize=20)
    #pl.xlim(distanceRange1[0],distanceRange1[-1])
    pl.xlabel('distance $d / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('Lorentz-factor $\gamma$', fontsize=24)
    pl.xscale('log')
    pl.yscale('log')
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def EvaluateToy1IntegrandOfEnergyLossRateTPP():
    ValuesForEpsilon = np.logspace(6,22,1000)*h # Values for the energy in units of J.
    ListOfMaxima = []
    # Linear plot:
    pl.figure(figsize=(16, 12), num="Integrand of TPP energy-loss-rate versus energy (linear)")
    pl.title('Integrand of TPP energy-loss-rate versus energy (linear)', fontsize=24)
    for gamma in np.logspace(2,16,8):
        ValuesForToy1IntegrandOfEnergyLossRateTPP = np.asarray([Toy1IntegrandOfEnergyLossRateTPP(i,gamma,ADAFNumberDensitySpectralInEnergy,ADAFSpectralLuminosity,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm) for i in ValuesForEpsilon])
        ListOfMaxima = np.append(ListOfMaxima,max(ValuesForToy1IntegrandOfEnergyLossRateTPP))
        pl.plot(ValuesForEpsilon/(me*c**2),ValuesForToy1IntegrandOfEnergyLossRateTPP, label='$\gamma = %g$' % gamma)
    Maximum = max(ListOfMaxima)
    pl.xscale('log')
    pl.legend(loc="best", fontsize=16)
    #pl.xlim(10**(-9),10**(0))
    pl.ylim(-Maximum,Maximum)
    pl.xlabel('energy / (me*c**2)', fontsize=24)
    pl.ylabel('Integrand of energy-loss-rate in 1/(m^3*J^2)', fontsize=24)
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)
    # Logarithmic plot:
    pl.figure(figsize=(16, 12), num="Integrand of TPP energy-loss-rate versus energy (logarithmic)")
    pl.title('Integrand of TPP energy-loss-rate versus energy (logarithmic)', fontsize=24)
    for gamma in np.logspace(2,16,8):
        ValuesForToy1IntegrandOfEnergyLossRateTPP = np.asarray([Toy1IntegrandOfEnergyLossRateTPP(i,gamma,ADAFNumberDensitySpectralInEnergy,ADAFSpectralLuminosity,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm) for i in ValuesForEpsilon])
        ListOfMaxima = np.append(ListOfMaxima,max(ValuesForToy1IntegrandOfEnergyLossRateTPP))
        pl.loglog(ValuesForEpsilon/(me*c**2),ValuesForToy1IntegrandOfEnergyLossRateTPP, label='$\gamma = %g$' % gamma)
    pl.legend(loc="best", fontsize=16)
    pl.xlabel('energy / (me*c**2)', fontsize=24)
    pl.ylabel('Integrand of energy-loss-rate in 1/(m^3*J^2)', fontsize=24)
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

def _EvaluateToy1SolvedgammaOverdd():
    Toy1DistanceRange, gammaRange, CurrentDistanceGap, CurrenthGap = Toy1SolvedgammaOverdd(-LRElectricFieldStrength(0.1,M9,0.1),5,0.1,1,0.025,0.02,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm) # Evaluation of Toy1SolvedgammaOverdd for Toy1ElectricFieldStrengthMax = -LRElectricFieldStrength(B4,M9,hGap).

def _EvaluatePlotToy1SolvedgammaOverdd():
    Fig = pl.figure(figsize=(12, 16), num="Figure 8: Lorentz-factor versus distance")
    Fig.suptitle("\nFigure 8: Lorentz-factor versus distance", fontsize=24)
    Plot1 = pl.subplot2grid((1, 2), (0, 0))
    Plot2 = pl.subplot2grid((1, 2), (0, 1))
    Plot1.plot(Toy1DistanceRange, gammaRange)
    Plot1.set_yscale('log')
    Plot2.loglog(Toy1DistanceRange, gammaRange, label='$\dot m = %g$, $T_\mathrm{e} = %g$' % (Dotm,ADAFTemperatureElectron))
    pl.legend(loc="best", fontsize=18)
    Plot1.set_ylabel('Lorentz-factor $\gamma$', fontsize=20)
    Plot1.set_ylim(gammaRange[-1]*0.9,max(gammaRange)*1.1)
    Plot2.set_ylim(gammaRange[-1]*0.9,max(gammaRange)*1.1)
    Plot1.set_xlim(CurrentDistanceGap-CurrenthGap/2,CurrentDistanceGap+3*CurrenthGap/2)
    Plot2.set_xlim(CurrentDistanceGap+3*CurrenthGap/2,Toy1DistanceRange[-1])
    Plot2.set_yticklabels([], fontsize=18)
    Plot1.set_xticks([CurrentDistanceGap-CurrenthGap/2,CurrentDistanceGap+CurrenthGap/2,CurrentDistanceGap+3*CurrenthGap/2])
    Plot1.set_xticklabels(['$d_\mathrm{gap}-$\n$h_\mathrm{gap}/2$','$d_\mathrm{gap}+$\n$h_\mathrm{gap}/2$','$d_\mathrm{gap}+$      \n$3h_\mathrm{gap}/2$        '], fontsize=16)
    Plot1.tick_params(labelsize=18)
    Plot2.tick_params(labelsize=18)
    Plot1.set_xlabel('distance $d$', fontsize=20)
    Plot2.set_xlabel('distance $d/r_\mathrm{S}$', fontsize=20)
    pl.subplots_adjust(wspace=0.0)

def _EvaluateEnergyLossRates():
    ValuesForEnergyLossRateTPP=np.asarray([])
    for distanceValue, gammaValue in zip(Toy1DistanceRange,gammaRange):
        ValuesForEnergyLossRateTPP = np.append(ValuesForEnergyLossRateTPP,Toy1ADAFEnergyLossRateTPP(distanceValue,gammaValue,5,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))
    ValuesForToy1EnergyLossRateCurvature=np.asarray([])
    for distanceValue, gammaValue in zip(Toy1DistanceRange,gammaRange):
        ValuesForToy1EnergyLossRateCurvature = np.append(ValuesForToy1EnergyLossRateCurvature,Toy1EnergyLossRateCurvature(distanceValue,gammaValue,5,1,0.025,0.02,M9))
    ValuesForEnergyLossRateIC=np.asarray([])
    for distanceValue, gammaValue in zip(Toy1DistanceRange,gammaRange):
        ValuesForEnergyLossRateIC = np.append(ValuesForEnergyLossRateIC,Toy1ADAFEnergyLossRateIC(distanceValue,gammaValue,5,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))
    Fig = pl.figure(figsize=(12, 16), num="Figure 9: Energy-loss-rates versus distance")
    Fig.suptitle("\nFigure 9: Energy-loss-rates versus distance", fontsize=24)
    Plot1 = pl.subplot2grid((1, 2), (0, 0))
    Plot2 = pl.subplot2grid((1, 2), (0, 1))
    Plot1.plot(Toy1DistanceRange,ValuesForEnergyLossRateTPP, label='Triplet pair-production')
    Plot1.plot(Toy1DistanceRange,ValuesForToy1EnergyLossRateCurvature, label='Curvature-radiation')
    Plot1.plot(Toy1DistanceRange,ValuesForEnergyLossRateIC, label='IC-scattering')
    Plot1.set_yscale('log')
    Plot2.loglog(Toy1DistanceRange,ValuesForEnergyLossRateTPP, label='Triplet pair-production')
    Plot2.loglog(Toy1DistanceRange,ValuesForToy1EnergyLossRateCurvature, label='Curvature-radiation')
    Plot2.loglog(Toy1DistanceRange,ValuesForEnergyLossRateIC, label='IC-scattering')
    Plot1.set_ylabel('Absolute value $\\left| \\frac{\mathrm{d}E}{\mathrm{d}t} \\right|$ of energy-loss-rate in $\\frac{\mathrm{J}}{\mathrm{s}}$', fontsize=20)
    pl.legend(loc="best", fontsize=18)
    Plot1.set_ylim(min([ValuesForEnergyLossRateTPP[-1],ValuesForToy1EnergyLossRateCurvature[-1],ValuesForEnergyLossRateIC[-1]])*0.9,max([max(ValuesForEnergyLossRateTPP),max(ValuesForToy1EnergyLossRateCurvature),max(ValuesForEnergyLossRateIC)])*1.1)
    Plot2.set_ylim(min([ValuesForEnergyLossRateTPP[-1],ValuesForToy1EnergyLossRateCurvature[-1],ValuesForEnergyLossRateIC[-1]])*0.9,max([max(ValuesForEnergyLossRateTPP),max(ValuesForToy1EnergyLossRateCurvature),max(ValuesForEnergyLossRateIC)])*1.1)
    Plot1.set_xlim(CurrentDistanceGap-CurrenthGap/2,CurrentDistanceGap+3*CurrenthGap/2)
    Plot2.set_xlim(CurrentDistanceGap+3*CurrenthGap/2,Toy1DistanceRange[-1])
    Plot2.set_yticklabels([], fontsize=18)
    Plot1.set_xticks([CurrentDistanceGap-CurrenthGap/2,CurrentDistanceGap+CurrenthGap/2,CurrentDistanceGap+3*CurrenthGap/2])
    Plot1.set_xticklabels(['$d_\mathrm{gap}-$\n$h_\mathrm{gap}/2$','$d_\mathrm{gap}+$\n$h_\mathrm{gap}/2$','$d_\mathrm{gap}+$      \n$3h_\mathrm{gap}/2$        '], fontsize=16)
    Plot1.tick_params(labelsize=18)
    Plot2.tick_params(labelsize=18)
    Plot1.set_xlabel('distance $d$', fontsize=20)
    Plot2.set_xlabel('distance $d/r_\mathrm{S}$', fontsize=20)
    pl.subplots_adjust(wspace=0.0)

def _EvaluateTeffOfEquation3dot50ByKato():
    rRange1=np.linspace(6,100,2000)
    ValuesForTeffOfEquation3dot50ByKato1=TeffOfEquation3dot50ByKato(M9,etaff,0.0002,6,rRange1)
    ValuesForTeffOfEquation3dot50ByKato2=TeffOfEquation3dot50ByKato(M9,etaff,0.0001,6,rRange1)
    rRange2=np.linspace(3,100,2000)
    ValuesForTeffOfEquation3dot50ByKato3=TeffOfEquation3dot50ByKato(M9,etaff,0.0002,3,rRange2)
    ValuesForTeffOfEquation3dot50ByKato4=TeffOfEquation3dot50ByKato(M9,etaff,0.0001,3,rRange2)
    pl.figure(figsize=(12, 9), num="Effective temperature of an accretion disk versus distance")
    pl.title('Figure 1: Effective temperature of\nan accretion disk versus distance', fontsize=24)
    pl.subplot(1,1,1)
    pl.plot(rRange1,ValuesForTeffOfEquation3dot50ByKato1, label='$\dot m = 0.0002$, $r_{\mathrm{in}} = 6 r_{\mathrm{S}}$')
    pl.plot(rRange1,ValuesForTeffOfEquation3dot50ByKato2, label='$\dot m = 0.0001$, $r_{\mathrm{in}} = 6 r_{\mathrm{S}}$')
    pl.plot(rRange2,ValuesForTeffOfEquation3dot50ByKato3, label='$\dot m = 0.0002$, $r_{\mathrm{in}} = 3 r_{\mathrm{S}}$')
    pl.plot(rRange2,ValuesForTeffOfEquation3dot50ByKato4, label='$\dot m = 0.0001$, $r_{\mathrm{in}} = 3 r_{\mathrm{S}}$')
    pl.text(30,11000,r'For $M_9 =$ %g and $\eta_{\mathrm{ff}} =$ %3g' % (M9, etaff), fontsize=20)
    pl.legend(loc="right", fontsize=20)
    pl.xlabel('radius $r / r_{\mathrm{S}}$', fontsize=24)
    pl.ylabel('effective temperature $T_{\mathrm{eff}} / \mathrm{K}$', fontsize=24)
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)
    pl.subplots_adjust(left=0.2)


