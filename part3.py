## ------------------------------------------------------------------------------------      
## Reconciling vacuum-gap-emission with electromagnetic cascades (toy-model 2)
##                                                                                           
## Implementation of various estimations (intended for the case of Mrk501)
## ------------------------------------------------------------------------------------

# This is the updated version of "Mrk501_estimates_version1.py" and now it is incorporated into the complete framework of code. 

import os
import multiprocessing

from part2 import * # Import of the file, whose path was added above.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n')
    print('\n_________________________________________________________________________________________\n')
    print('-------------------    Vacuum-gap-emission and IC-pair-cascades    ----------------------')
    print('_________________________________________________________________________________________\n\n')


## Input parameters

Toy2tauPPOnADAF = 10**(-5.0) # This is the optical depth of a photon with respect to pair-production in an ADAF-photon-field, it can be computed with my code "Comparison_of_optical_depth_version5.py". INPUT VALUE!
Toy2TestParticleInjectionRate = 7.0*10.0**33 # The total injection rate of test-particles into the gap. Unit is 1/s. INPUT VALUE!

## Computations

# Determine the limitation of the Lorentz-factor within the model of 2011ApJ...730..123L:
Toy2gammaCurv = LRgammaCurOfEquation13(B4,M9,hGap,rCurvature) # Limitation due to curvature-radiation.
Toy2gammaIC = LRgammaICOfEquation14(B4,hGap,ADAFTotalEnergyDensityNum*10) # Limitation due to IC-scattering of the ADAF-photons.

# Determine the dominating process, the maximum Lorentz-factor and the photon-energy:
if Toy2gammaIC<=Toy2gammaCurv:
    Toy2DominatingProcess = 'IC'
    Toy2gammaMax = Toy2gammaIC
    Toy2PhotonEnergy = me*c**2*Toy2gammaMax # In units of J. This assumes scattering inn the KN-regime.
else:
    Toy2DominatingProcess = 'Curvature'
    Toy2gammaMax = Toy2gammaCurv
    Toy2PhotonEnergy = ECritOfCurvatureRadiation(M9,Toy2gammaMax,rCurvature) # In units of J.

# Consider the voltage in the vacuum-gap according to 2011ApJ...730..123L:
Toy2GapVoltage = LRVoltage(B4,M9,hGap)
Toy2GapEnergy = Toy2GapVoltage*e # Energy in units of J.

# PTP = per test-particle
Toy2PTPNumberOfPhotons = Toy2GapEnergy/Toy2PhotonEnergy # The averaged number of photons, that are emitted by one test-particle.
Toy2PTPMultiplicity = (1.0-np.exp(-Toy2tauPPOnADAF))*Toy2PTPNumberOfPhotons # The number of photons, that are converted to pairs, per test-particle.
Toy2PTPNumberOfParticles = 2.0*Toy2PTPMultiplicity # The number of particles, that are produced, per test-particle.
Toy2PTPNumberOfPhotonsRemaining = Toy2PTPNumberOfPhotons-Toy2PTPMultiplicity # The averaged number of emitted photons, that are not absorbed via PP, per test-particle.

Toy2ParticleEnergy = Toy2PhotonEnergy/2.0 # The averaged energy of one pair-produced particle in units of J.

# Assume an interaction region:
Toy2DiametralExtension = rSchwarzschild(M9) # The approximate diameter of the interaction region in units of m.
Toy2InteractionVolume = Toy2DiametralExtension**3 # The corresponding volume in units of m^3.

# The injection-density per test-particle in the interaction-region:
Toy2PTPInjectionDensityParticles = Toy2PTPNumberOfParticles/Toy2InteractionVolume # Unit is 1/(m^3).
Toy2PTPInjectionDensityPhotons = Toy2PTPNumberOfPhotonsRemaining/Toy2InteractionVolume # Unit is 1/(m^3).

Toy2TestMassInjectionRate = Toy2TestParticleInjectionRate*me # The total mass, that is injected per unit time into the gap. Unit is kg/s.

# Consider the total number of species, that are injected into the interaction-region per unit space interval and per unit time.
Toy2InjectionDensityRateParticles = Toy2PTPInjectionDensityParticles*Toy2TestParticleInjectionRate # Of particles. Unit is 1/(m^3 s).
Toy2InjectionDensityRatePhotons = Toy2PTPInjectionDensityPhotons*Toy2TestParticleInjectionRate # Of photons. Unit is 1/(m^3 s).


if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('----------- Reconciling vacuum-gap-emission with IC-pair-cascades for Mrk501 -------------\n')
    print('Energy-density of ADAF-photons: Numerically: ', ADAFTotalEnergyDensityNum, 'J/m^3')
    print('                                Analytically:', ADAFTotalEnergyDensityAna, 'J/m^3\n')
    print('Lorentz-factor (curvature-radiation limited):', Toy2gammaCurv)
    print('Lorentz-factor (IC-radiation (ADAF) limited):', Toy2gammaIC)
    if Toy2gammaIC<=Toy2gammaCurv:
        print('IC-radiation dominates.\n')
    else:
        print('Curvature-radiation dominates.\n')
    print('The emitted photons have the energy', Toy2PhotonEnergy/e, 'eV, which corresponds to a dimensionless energy of', Toy2PhotonEnergy/(me*c**2),'.') 
    print('\nThe potential drop along the gap is', Toy2GapVoltage, 'V.')
    print('The energy attainable by an electron via the electrostatic potential is', Toy2GapEnergy/e, 'eV.')
    print('\nThe number of emitted photons (per test-particle) is', Toy2PTPNumberOfPhotons, '.')
    print('\nThe optical depth of a photon with respect to PP on the ADAF-photons is', Toy2tauPPOnADAF, '(INPUT PARAMETER!).')
    print('\nThe particle multiplicity is', Toy2PTPMultiplicity, '.')
    print('\nNumber of particles, that are created by one original particle after one generation:', Toy2PTPNumberOfParticles)
    print('Number of photons, that are created by one original particle after one generation and are not pair-absorbed:', Toy2PTPNumberOfPhotonsRemaining)
    print('\nThe created particles have the energy', Toy2ParticleEnergy/e, 'eV, which corresponds to a Lorentz-factor of', Toy2ParticleEnergy/(me*c**2),'.')
    print('\nThe diametral extension of the beam / interaction-volume is assumed to be', Toy2DiametralExtension/rSchwarzschild(M9), 'r_S.\n')
    print('Per test-particle injection-density: For particles:', Toy2PTPInjectionDensityParticles, '/(m^3)')
    print('                                     For photons:  ', Toy2PTPInjectionDensityPhotons, '/(m^3)')
    print('\nThe total injection-rate of test-particles into the gap is', Toy2TestParticleInjectionRate, '/s (INPUT PARAMETER!)...')
    print('... corresponding to a total mass-injection-rate into the gap of', Toy2TestMassInjectionRate, 'kg/s.\n')
    print('Total injection-density-rate: For particles:', Toy2InjectionDensityRateParticles, '/(m^3 s)')
    print('                              For photons:  ', Toy2InjectionDensityRatePhotons, '/(m^3 s)')