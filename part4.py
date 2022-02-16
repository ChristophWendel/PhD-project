## ---------------------------------------------------------------------------------------------     o   o
## Electromagnetic cascades                                                                            I
##                                                                                                   \.../
## Implementation of various formulae on synchrotron-radiation, IC-scattering and pair-production
## -------------------------------------------------------------------------------------------------------

import os
import multiprocessing

from part3 import * # Import of the file, whose path was added above.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n')
    print('\n_________________________________________________________________________________________\n')
    print('-------------------    Synchrotron-Inverse-Compton-pair-cascades    ---------------------')
    print('_________________________________________________________________________________________\n')


## Initialisation

NumberOfIteratingProcesses = 2 # The number of processes to be initialised in the EvaluateIteration and in SecondTermOfNumeratorOfEq8. For values other than 1, multiprocessing is used in these computations.
UseMultiprocessingIn4thTermOfEq1 = True # This states whether multiprocessing is used in the computation of FourthTermOfEq1.
RunIdentifier = 'Mrk501 MJD56857' # A name for the current computation-run.

UsedIterationScheme = 'PointsFromRToLIterationPointwise' # Choose the used iteration scheme. Possible values are 'IterationStepByStepPointsFromLToR' and 'PointsFromRToLIterationPointwise'. Confer "Idea for iteration schemes.png" for an explanation of the iteration schemes.
UsedAlternativeFor4thTermOfEq1 = 1 # This affects the used alternative of FourthTermOfEq1 in part2 in the case UsedDotni is NOT DotniDelta or DotniZero. Choose between 1 which will cause the use of Alternative1 or 2 for Alternative2, here. In "Chronology Modelling 3C279" in the 67th run, it was yielded that Alternative1 is slightly faster than Alternative2. INPUT VALUE!
UsedSampleIntBordersIn4thTerm = SampleIntegrationBorders11 # This specifies which version of the function SampleIntegrationBorders is actually being used in FourthTermOfEq1. In "Chronology on IC-pair-sync-cascades" in the 43th run, it was yielded that SampleIntegrationBorders9 is faster than SampleIntegrationBorders7. INPUT VALUE!
UsedSampleIntBordersOfxSync = SampleIntegrationBorders12 # This specifies which version of the function SampleIntegrationBorders is actually being used in the sampling of ValuesForxSyncSuperIt. In "Chronology on IC-pair-sync-cascades" in the 43th till 45th run, it was yielded that SampleIntegrationBorders3 is faster than 1 and 5. SampleIntegrationBorders12 is even better for the sampling of xSync, it has even less borders and just one additional border in the highest range, where nSyncPs has an exponential cut-off. INPUT VALUE!

IncludeHEPhotonEscape=True # A boolean which determines whether escape of the high-energy photons is taken into account. INPUT VALUE!
IncludeElectronEscape=True # A boolean which determines whether electron escape is taken into account. INPUT VALUE!
IncludeSynchrotron=False # A boolean which determines whether synchrotron-radiation is taken into account. INPUT VALUE!
PerformanceMode='CommonIteration' # This can attain the values 'Initialisation', 'CommonIteration' and 'SuperIteration'. 'SuperIteration' states that a super-iteration with successive common iterations including synchrotron-radiation is performed. The execution of a super-iteration might be feasible only together with 'PointsFromRToLIterationPointwise' and only with IncludeSynchrotron=True. 'CommonIteration' states that a common iteration is done, with or without synchrotron-losses in the electrons' kinetic equation but in any case without back-reaction of the SyncPs as background photons. 'Initialisation' just starts the programme and is intended for use in Pyzo. INPUT VALUE!
InjectionType = 'Chosen' # This can be 'Chosen', which means that the electron injection and the HEP injection is specified as a chosen function, or 'ResultsOfOldIt', which means that the results of a previous iteration are used as injection. In the later case all input values have to be equalised as far as necessary to the old iteration run, however, and the files to be imported have to be specified, of course. 'ResultsOfOldIt' is compatible only with 'PointsFromRToLIterationPointwise'. INPUT VALUE! 

if IncludeSynchrotron==False and PerformanceMode=='SuperIteration':
    raise InvalidSetOfParametersError('A super-iteration is reasonable only with inclusion of synchrotron-radiation.')

if IncludeSynchrotron==True and UsedIterationScheme=='IterationStepByStepPointsFromLToR':
    raise InvalidSetOfParametersError('Inclusion of synchrotron-radiation is allowed only together with "PointsFromRToLIterationPointwise".')


## Import of old iterations

if "MainProcess" in multiprocessing.current_process().name and UsedIterationScheme == 'IterationStepByStepPointsFromLToR': # The first condition prevents evaluation in child-processes.
    def EvaluateImportDataOfIteration():
        '''The string of a data-file is imported here to extract IterationCounter, and all the arrays ValuesForgaIt and ValuesForNextNElectrons of each iteration. For each iteration, the interpolated object is created again as above based on ValuesForgaIt and ValuesForNextNElectrons. The name of the file has to be inserted.'''
        OpenedLines = ImportFileToString('Run Mrk501 MJD56857, N(gamma) versus gamma from 2022-02-15 10-29 (data).dat')
        IterationCounter = eval(OpenedLines[0].replace("\n", "")[24:]) # Get the current value of IterationCounter.
        NElectronsIteratedString = '' # Create a string that will contain the dictionary NElectronsIterated.
        for Line in OpenedLines[1:]:
            NElectronsIteratedString = NElectronsIteratedString+Line # Join all remaining lines together.
        NElectronsIteratedString=NElectronsIteratedString.replace("array(", "") # The string is cleaned and ...
        NElectronsIteratedString=NElectronsIteratedString.replace(")", "")
        NElectronsIteratedString=NElectronsIteratedString.replace("\n", "")
        while "<sci" in NElectronsIteratedString:
            Start=NElectronsIteratedString.find("<sci") # ... the part-strings, that denoted the interpolated objects, are removed, so that...
            End=NElectronsIteratedString.find(">")
            Substring=NElectronsIteratedString[Start:End+1]
            NElectronsIteratedString=NElectronsIteratedString.replace(Substring, "0")
        while "  " in NElectronsIteratedString:
            NElectronsIteratedString=NElectronsIteratedString.replace("  ", " ")
        NElectronsIterated=eval(NElectronsIteratedString) # ...eval() can be used to create a dictionary again. However, the dicts in this dict contain lists of numbers, while originally ValuesForgaIt and ValuesForNextNElectrons have been numpy-arrays. Thus, take each list, convert it into an array, and save it instead of the list:
        for ItCounter in range(IterationCounter):
            for Item in ['Values for gamma','Values for next NElectrons']:
                NElectronsIterated['%s. Iteration' % ItCounter][Item] = np.asarray(NElectronsIterated['%s. Iteration' % ItCounter][Item])
            # By the way, restore the interpolated objects:
            NElectronsIterated['%s. Iteration' % ItCounter]['Interpolated object'] = interp1d(NElectronsIterated['%s. Iteration' % ItCounter]['Values for gamma'], NElectronsIterated['%s. Iteration' % ItCounter]['Values for next NElectrons'], kind='linear', bounds_error=False, fill_value=0.0)
        # Last but not least, restore the values for the next iteration:
        ValuesForgaIt = NElectronsIterated['%s. Iteration' % (IterationCounter-1)]['Values for gamma']
        ValuesForNextNElectrons = NElectronsIterated['%s. Iteration' % (IterationCounter-1)]['Values for next NElectrons']
        NElectrons = NElectronsIterated['%s. Iteration' % (IterationCounter-1)]['Interpolated object']
        LastValuesForgaIt = ValuesForgaIt
        return NElectrons, ValuesForgaIt, LastValuesForgaIt, ValuesForNextNElectrons, NElectronsIterated, IterationCounter
elif "MainProcess" in multiprocessing.current_process().name and UsedIterationScheme == 'PointsFromRToLIterationPointwise': # The first condition prevents evaluation in child-processes.
    def EvaluateImportDataOfIteration(NameOfFileToImportForNElectrons,NameOfFileToImportFornSyncPs=''):
        '''The string of a data-file is imported here to extract the current array ValuesForNElectrons and the dictionary NElectronsIterated or NElectronsSuperIt, respectively. The interpolated object is recreated again as above based on ValuesForgaIt and ValuesForNElectrons. DictOfNElectrons can either be NElectronsIterated or NElectronsSuperIt. The name of the file, in which NElectrons was stored, has to be given as argument, e.g. 'N(gamma) versus gamma from 2020-03-24 10-23 (data).dat'. If a super-iteration is imported, the name of the file, in which nSyncPs was stored, has to be inserted as well, e.g. 'nSyncPs versus xSync from 2020-02-19 19-23 (SuperIt=1, data).dat'.
        In the case of importing an unfinished super-iteration, one has to do the following: Take the latest temporary file, containing the results of the latest common iteration. Manually paste its content, together with the prefix ", c: ", where c is the number of the unfinished super-iteration-step, into the latest file, containing the results of the super-iteration. Then, manually increment the SuperItCounter in the filename by 1.'''
        if "SuperIt" in NameOfFileToImportForNElectrons:
            SuperItCounter=int(NameOfFileToImportForNElectrons.split('SuperIt=')[1].split(', data')[0])
        OpenedLines = ImportFileToString(NameOfFileToImportForNElectrons)
        DictOfNElectronsString = '' # Create a string that will contain the dictionary DictOfNElectrons.
        for Line in OpenedLines:
            DictOfNElectronsString = DictOfNElectronsString+Line # Join all remaining lines together.
        DictOfNElectronsString=DictOfNElectronsString.replace("array(", "") # The string is cleaned and ...
        DictOfNElectronsString=DictOfNElectronsString.replace(")", "")
        DictOfNElectronsString=DictOfNElectronsString.replace("\n", "")
        while "  " in DictOfNElectronsString:
            DictOfNElectronsString=DictOfNElectronsString.replace("  ", " ")
        DictOfNElectrons=eval(DictOfNElectronsString) # ...eval() can be used to recreate a real dictionary again from the string.
        # However, the dicts in this dict contain lists of numbers, while originally StackOfIteratedValues have been numpy-arrays. Thus, take each list, convert it into an array, and save it in the dict instead of the list. Afterwards, recreate the array of current values of NElectrons:
        if "SuperIt" not in NameOfFileToImportForNElectrons:
            ValuesForgaIt = np.sort(np.asarray([ga for ga in DictOfNElectrons.keys()])) # Recreate the array with the values of ga.
            for ga in ValuesForgaIt:
                DictOfNElectrons[ga]['StackOfIteratedValues'] = np.asarray(DictOfNElectrons[ga]['StackOfIteratedValues'])
            ValuesForNElectrons = np.asarray([DictOfNElectrons[ga]['StackOfIteratedValues'][-1] for ga in ValuesForgaIt])
        elif "SuperIt" in NameOfFileToImportForNElectrons:
            if 'temporary' in NameOfFileToImportForNElectrons:
                ValuesForgaIt = np.sort(np.asarray([ga for ga in DictOfNElectrons.keys()])) # Recreate the arrays with the values of ga.
                IntermediateStore = {}
                for ga in ValuesForgaIt:
                    IntermediateStore[ga] = DictOfNElectrons[ga]
                DictOfNElectrons = {SuperItCounter: {}}
                for ga in ValuesForgaIt:
                    DictOfNElectrons[SuperItCounter][ga] = IntermediateStore[ga]
                for ga in ValuesForgaIt:
                    DictOfNElectrons[SuperItCounter][ga]['StackOfIteratedValues'] = np.asarray(DictOfNElectrons[SuperItCounter][ga]['StackOfIteratedValues'])
            else:
                for c in range(1,SuperItCounter+1):
                    ValuesForgaIt = np.sort(np.asarray([ga for ga in DictOfNElectrons[c].keys()])) # Recreate the arrays with the values of ga.
                    for ga in ValuesForgaIt:
                        DictOfNElectrons[c][ga]['StackOfIteratedValues'] = np.asarray(DictOfNElectrons[c][ga]['StackOfIteratedValues'])
            ValuesForNElectrons = np.asarray([DictOfNElectrons[SuperItCounter][ga]['StackOfIteratedValues'][-1] for ga in ValuesForgaIt])
        Lowestga = ValuesForgaIt[0] # Recreate Lowestga.
        NElectronsga0 = ValuesForgaIt[-1] # Recreate NElectronsga0.
        LastValuesForgaIt = ValuesForgaIt 
        NElectrons = interp1d(ValuesForgaIt, ValuesForNElectrons, kind='linear', bounds_error=False, fill_value=0.0) # Recreate the interpolation.
        pl.figure(figsize=(18, 14), num="Cascade-equation: N(gamma) versus particle energy (Imported)") # Plot the imported results.
        ValuesForgamma = np.logspace(np.log10(Lowestga),np.log10(NElectronsga0),5000)
        pl.loglog(ValuesForgamma,ValuesForgamma*NElectrons(ValuesForgamma), label='CurrentNElectrons')
        if "SuperIt" in NameOfFileToImportForNElectrons:
            OpenedLines2 = ImportFileToString(NameOfFileToImportFornSyncPs)
            nSyncPsSuperItString = '' # Create a string that will contain the dictionary nSyncPsSuperIt.
            for Line in OpenedLines2:
                nSyncPsSuperItString = nSyncPsSuperItString+Line # Join all remaining lines together.
            nSyncPsSuperItString=nSyncPsSuperItString.replace("array(", "") # The string is cleaned and ...
            nSyncPsSuperItString=nSyncPsSuperItString.replace(")", "")
            nSyncPsSuperItString=nSyncPsSuperItString.replace("\n", "")
            while "  " in nSyncPsSuperItString:
                nSyncPsSuperItString=nSyncPsSuperItString.replace("  ", " ")
            nSyncPsSuperIt=eval(nSyncPsSuperItString) # ...eval() can be used to recreate a real dictionary again from the string.
            # However, the dicts in this dict contain lists of numbers, while originally there have been numpy-arrays. Thus, take each list, convert it into an array, and save it in the dict instead of the list. Afterwards, recreate the array of current values of nSyncPs:
            for c in range(SuperItCounter+1):
                if c in nSyncPsSuperIt.keys(): # In the case of the c-th super-iteration being unfinished, nSyncPsSuperIt has only keys until SuperItCounter-1, so evaluate only those keys that are present in nSyncPsSuperIt:
                    nSyncPsSuperIt[c]['ValuesForxSyncSuperIt'] = np.asarray(nSyncPsSuperIt[c]['ValuesForxSyncSuperIt'])
                    nSyncPsSuperIt[c]['ValuesFornSyncPs'] = np.asarray(nSyncPsSuperIt[c]['ValuesFornSyncPs'])
            if SuperItCounter in nSyncPsSuperIt.keys():
                ValuesForxSyncSuperIt = nSyncPsSuperIt[SuperItCounter]['ValuesForxSyncSuperIt']
                ValuesFornSyncPs = np.asarray(nSyncPsSuperIt[SuperItCounter]['ValuesFornSyncPs'])
            else:
                ValuesForxSyncSuperIt = nSyncPsSuperIt[SuperItCounter-1]['ValuesForxSyncSuperIt']
                ValuesFornSyncPs = np.asarray(nSyncPsSuperIt[SuperItCounter-1]['ValuesFornSyncPs'])
            nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0) # Recreate the interpolation.
            xSync0 = DeterminexSync0(NElectronsga0) # Recreate xSync0.
            xSync1 = ValuesForxSyncSuperIt[0] # Recreate xSync1.
            xSync0Used = ValuesForxSyncSuperIt[-1] # Recreate xSync0Used.
            xTotal0 = DeterminexTotal0(xSync0Used)  # Recreate xTotal0.
            pl.figure(figsize=(18, 14), num="KESP: nSyncPs(xSync) versus photon-energy xSync (Imported)") # Plot the imported results.
            ValuesForxSync = np.logspace(np.log10(ValuesForxSyncSuperIt[0]),np.log10(ValuesForxSyncSuperIt[-1]),2000)
            pl.loglog(ValuesForxSync,ValuesForxSync*nSyncPs(ValuesForxSync))
            pl.xlabel("Synchrotron-photons dimensionless energy $x_{Sync}$", fontsize=24)
            pl.ylabel("SyncPs dim.less energy times\nspectral number-density $x_{Sync} \cdot n_{SyncPs}$", fontsize=24)
            if 'temporary' in NameOfFileToImportForNElectrons:
                return NElectrons, ValuesForgaIt, LastValuesForgaIt, ValuesForNElectrons, DictOfNElectrons[SuperItCounter], NElectronsga0, Lowestga, SuperItCounter, nSyncPs, ValuesForxSyncSuperIt, ValuesFornSyncPs, nSyncPsSuperIt, xSync0, xSync1, xSync0Used, xTotal0
            else:
                return NElectrons, ValuesForgaIt, LastValuesForgaIt, ValuesForNElectrons, DictOfNElectrons, NElectronsga0, Lowestga, SuperItCounter, nSyncPs, ValuesForxSyncSuperIt, ValuesFornSyncPs, nSyncPsSuperIt, xSync0, xSync1, xSync0Used, xTotal0
        else:
            return NElectrons, ValuesForgaIt, LastValuesForgaIt, ValuesForNElectrons, DictOfNElectrons

def EvaluateImportDataOfHEPs(nHEPhotonsSumSpectralInxgaFileName):
    nHEPhotonsSumSpectralInxgaFilePath = SearchAFile0rDirectoryEverywhere(nHEPhotonsSumSpectralInxgaFileName,False)
    nHEPhotonsSumSpectralInxgaData = np.genfromtxt(os.path.join(nHEPhotonsSumSpectralInxgaFilePath,nHEPhotonsSumSpectralInxgaFileName), skip_header=1) # This reads all the data (usual points as well as upper limits) in the file into an array.
    ValuesForxgammanHEPhotons = nHEPhotonsSumSpectralInxgaData[:,0] # Select the first column, which is the energy.
    ValuesFornHEPhotonsSumSpectralInxga = nHEPhotonsSumSpectralInxgaData[:,1] # Select the second column, which is the summed number-density.
    ValuesFornHEPhotonsSpectralInxga = nHEPhotonsSumSpectralInxgaData[:,2] # Select the third column, which is the number-density of the injected HEPs.
    nHEPsInterpolated = interp1d(ValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, kind='linear', bounds_error=False, fill_value=0.0) # Create an interpolation.
    return ValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, ValuesFornHEPhotonsSpectralInxga, nHEPsInterpolated

if InjectionType=='ResultsOfOldIt':
    OldCurrentNElectrons, OldValuesForgaIt, CurrentLastValuesForgaIt, CurrentValuesForNElectrons, NElectronsIterated = EvaluateImportDataOfIteration('N(gamma) versus gamma from 2020-03-24 10-23 (data).dat')
    OldValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, ValuesFornHEPhotonsSpectralInxga, OldnHEPsInterpolated = EvaluateImportDataOfHEPs('nHEP versus energy from 2020-03-24 10-24.dat')
else: # This is just for the sake of existence.
    def OldCurrentNElectrons(ga):
        return 0
    OldValuesForgaIt=np.asarray([0,1])
    OldValuesForxgammanHEPhotons=np.asarray([0,1])
    def OldnHEPsInterpolated(xga):
        return 0


## Rough specification of the setting

# Consider a typical homogeneous broad-line region:
ExPsInnerRadius = rSchwarzschildTom(0.05,M9) # The inner radius of the external photon-field, which is assumed to be a spherical shell, in units of m.
ExPsOuterRadius = rSchwarzschildTom(0.1,M9) # The outer radius of the external photon-field in units of m. INPUT VALUE!
ExPsMeanRadius = (ExPsInnerRadius+ExPsOuterRadius)/2.0 # The mean radius of the broad line region spherical shell in units of m.
ExPsVeryOuterRadius = 2*ExPsOuterRadius # The radius of the spherical shell, that is assumed to consider a diluted ExPs-field, in units of m. INPUT VALUE!

# Consider an AGN-jet:
JetBaseRadius = rSchwarzschildTom(0.1,M9) # The radius of the jet-base in units of m.
JetHalfOpeningAngle = 0.88 # The half of the opening angle of the jet, along which the injected species are injected, in units of °. INPUT VALUE!
JetInteractionRadius = JetBaseRadius+np.tan(JetHalfOpeningAngle/360.0*2*np.pi)*ExPsMeanRadius # The radius of the jet in the distance of the centre of the external photon-field, where the interaction region is assumed to be situated, in units of m.

# Consider the actual interaction region:
InteractionRadius = rSchwarzschildTom(0.1,M9) # This is the radius of the (spherical) volume within which the cascade takes place. It is in units of m. INPUT VALUE!
InteractionVolume = (4.0*np.pi*InteractionRadius**3)/3.0 # This is the spatial volume within which the cascade takes place. It is in units of m^3. INPUT VALUE!
MeanEscapeTime = InteractionRadius/(c) # This is the time, that is in mean needed by a photon (or relativistic electron) to escape from the interaction volume. It is in units of s.

def SphericalVolume(UsedR):
    '''This returns the spatial volume in units of m^3 of a sphere of radius called UsedR.'''
    if UsedR=='InteractionRadius':
        return UsedR, 4.0*np.pi*InteractionRadius**3/3.0
    elif UsedR=='ExPsOuterRadius':
        return UsedR, 4.0*np.pi*ExPsOuterRadius**3/3.0
    elif UsedR=='ExPsVeryOuterRadius':
        return UsedR, 4.0*np.pi*ExPsVeryOuterRadius**3/3.0
    elif UsedR=='RClouds':
        return UsedR, 4.0*np.pi*RClouds**3/3.0

EmittingSolidAngle = OpeningAngleToSolidAngle(2*JetHalfOpeningAngle) # This is the solid angle, into which the radiation is focused. It is in units of sterad.
DistanceSourceEarthMpc = 149.4 # The distance from the emitting source to the detector on planet earth. It is in units of Mpc. INPUT VALUE!
DistanceSourceEarth = MpcTom(DistanceSourceEarthMpc) # The distance from the emitting source to the detector on planet earth. It is in units of m.
Redshift3C279 = 0.5362

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n\n------------------------------- Geometry of the setting ---------------------------------\n')
    print("Radius of the jet's cross-section at BLR-centre:        %3.4g m = %3.4g pc = %5.1f r_S" % (JetInteractionRadius, mToMpc(JetInteractionRadius)*1000000, JetInteractionRadius/rSchwarzschild(M9)))
    print("Half-thickness of the BLR:                              %3.4g m = %3.4g pc = %5.1f r_S" % ((ExPsOuterRadius-ExPsInnerRadius)/2, mToMpc((ExPsOuterRadius-ExPsInnerRadius)/2)*1000000, (ExPsOuterRadius-ExPsInnerRadius)/(2*rSchwarzschild(M9))))
    print("Radius of the interaction region:                       %3.4g m = %3.4g pc = %5.1f r_S" % (InteractionRadius, mToMpc(InteractionRadius)*1000000, InteractionRadius/rSchwarzschild(M9)))
    print("Volume of the interaction region:                       %3.4g m^3" % SphericalVolume('InteractionRadius')[1])

def FactorEnergyDensityToLuminosity(UsedR):
    '''This returns a quantity in units of m^3/s. Any energy-density in units of J/m^3 has to be multiplied with this factor to yield the luminosity in W, that is escaping from the region of radial size called UsedR.'''
    if UsedR=='InteractionRadius':
        return UsedR, c*4.0*np.pi*InteractionRadius**2
    elif UsedR=='ExPsOuterRadius':
        return UsedR, c*4.0*np.pi*ExPsOuterRadius**2
    elif UsedR=='ExPsVeryOuterRadius':
        return UsedR, c*4.0*np.pi*ExPsVeryOuterRadius**2
    elif UsedR=='RClouds':
        return UsedR, c*4.0*np.pi*RClouds**2

if InjectionType=='Chosen':
    FactorNumberDensityToInjectionRate = 0.0 # This is just for the sake of completeness.
    FactornHEPsToFlux = c*4.0*np.pi*ExPsOuterRadius**2/(EmittingSolidAngle*DistanceSourceEarth**2*me*c**2) # This is the factor, which is multiplied with a spectral number-density to yield a flux (cf. part7). It is in units of m/(s*J). Here, it is assumed that that radiation, which would leave a spherical source, is emitted only into the EmittingSolidAngle. INPUT VALUE!
elif InjectionType=='ResultsOfOldIt':
    FactorNumberDensityToInjectionRate = c*EmittingSolidAngle*ExPsOuterRadius**2/SphericalVolume('InteractionRadius')[1] # If a number-density of an old iteration is used as injection of the new iteration, the imported number-density has to be multiplied with this factor to yield the injection rate. This factor has to be in units of 1/s and it seems reasonable that it is something like c*Area/Volume. INPUT VALUE!
    FactornHEPsToFlux = c*ExPsVeryOuterRadius**2/(DistanceSourceEarth**2*me*c**2) # Here it is assumed, that the radiation is emitted only into the EmittingSolidAngle and that it leaves the source at the distance ExPsVeryOuterRadius. INPUT VALUE!
    if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
        print("Conversion factor for number-density to injection-rate: %3.4g /s" % FactorNumberDensityToInjectionRate)
if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print("Conversion factor for number-density to flux-density:   %3.4g m/(s*J)" % FactornHEPsToFlux)


## Injection of primary particles  

# Various physically reasonable possibilities for the injected particles are given now. The distributions in ga are not extending unlimited but are cutted-off at certain values of ga via a Heaviside-function.

InjectedalphaPL = -2.0    # The power-law-index. INPUT VALUE!
InjectedalphaExp = -0.0003    # The exponent's coefficient. INPUT VALUE!
InjectedTemperature = 10.0**11.5       # The temperature in units of K! INPUT VALUE!
InjectedTheta = Theta(InjectedTemperature)  # The thermal energy in units of electron rest-mass energy.
InjectedGaussianLocation = 3.4*10**(12)/(me*(c**2)/e)  # The value of ga, where the Gaussian distribution peaks, i. e. the mean. INPUT VALUE! Assume that particles with energy 3.4 TeV are injected.
InjectedGaussianWidth = 0.23*InjectedGaussianLocation       # The width of the distribution. INPUT VALUE!
if InjectionType=='ResultsOfOldIt':
    Injectedga0 = OldValuesForgaIt[-1] # The upper cut-off of the distribution of the injection-rate. INPUT VALUE!
    Injectedga1 = OldValuesForgaIt[0] # The lower cut-off of the distribution. INPUT VALUE!
elif InjectionType=='Chosen':
    Injectedga0 = InjectedGaussianLocation+3.0*InjectedGaussianWidth # The upper cut-off of the distribution of the injection-rate. INPUT VALUE!
    Injectedga1 = InjectedGaussianLocation-3.0*InjectedGaussianWidth # The lower cut-off of the distribution. INPUT VALUE!
# The Gaussian's coefficient and the Gaussian's norm are specified, now. Consider now that the relation Norm=Coefficient*(Width*np.sqrt(2*np.pi)) holds. Thus, if Coefficient is given, Norm can be determined and if Norm is given, Coefficient can be determined.
UseInjectedGaussianCoefficientManualInput = False # A boolean, that states whether a self chosen value for Coefficient (True) or a self chosen value for Norm (False) is used. INPUT VALUE!
InjectedGaussianNormManualInput = 3.3*10.0**4 # This is the norm of a Gaussian distribution, if InjectedGaussianNorm/(InjectedGaussianWidth*np.sqrt(2*np.pi)) is used as its coefficient. INPUT VALUE!
InjectedGaussianCoefficientManualInput = 10.0**(5.0) # A self chosen coefficient for dot N_i being a Gaussian distribution. It is only applied, if...  INPUT VALUE!
if UseInjectedGaussianCoefficientManualInput == True: # ... this evaluates to true.
    InjectedGaussianCoefficient = InjectedGaussianCoefficientManualInput
    InjectedGaussianNorm = InjectedGaussianCoefficient*(InjectedGaussianWidth*np.sqrt(2.0*np.pi))
elif UseInjectedGaussianCoefficientManualInput == False: # Otherwise, a value, that normalises the Gaussian to InjectedGaussianNorm, is applied.
    InjectedGaussianNorm = InjectedGaussianNormManualInput
    InjectedGaussianCoefficient = InjectedGaussianNorm/(InjectedGaussianWidth*np.sqrt(2.0*np.pi))
InjectedExpCoefficient = 10.0**(7.0) # The coefficient for dot N_i being an exponential-function. INPUT VALUE!
InjectedPLCoefficient = 10.0**(14.0) # The coefficient for dot N_i being a power-law. INPUT VALUE!
InjectedHeaviCoefficient = 10.0**(7.0) # The coefficient for dot N_i being a Heaviside-function. INPUT VALUE!
InjectedRelMaxCoefficient = 10000.0 # The coefficient for dot N_i being a Maxwell-Juettner-distribution. It is in units of 1/(s*m^3). INPUT VALUE!

# The following quantities correspond to dot N_i(ga). (ga is the dimensionless energy of the particles.) According to my analysis, it is an injection-rate, i. e. a spectral number-density per time-interval. Thus, it is in units of 1/(s*m^3).

def DotNiExp(ga):  # An exponential-function is used here.
    return InjectedExpCoefficient*np.exp(InjectedalphaExp*ga)*Heavi(Injectedga1,ga,Injectedga0)
# Plotting-test:    X=np.logspace(2,5,301);Y=DotNiExp(X);pl.loglog(X,Y)#;pl.xscale('log')

def DotNiPL(ga):  # A power-law is used here.
    return InjectedPLCoefficient*ga**(InjectedalphaPL)*Heavi(Injectedga1,ga,Injectedga0)
# Plotting-test:    X=np.logspace(5,8,301);Y=DotNiPL(X);pl.loglog(X,Y)

def DotNiGauss(ga):    # A Gaussian distribution is used here. It is not necessarily normalised. 
    return InjectedGaussianCoefficient*np.exp(-0.5*((ga-InjectedGaussianLocation)/InjectedGaussianWidth)**2.0)*Heavi(Injectedga1,ga,Injectedga0)
# Plotting-test:    X=np.logspace(5,7,901);Y=DotNiGauss(X);pl.loglog(X,Y)

def DotNiGaussPlusPL(ga):    # A Gaussian distribution is used here. It is not necessarily normalised. 
    return (InjectedPLCoefficient*ga**InjectedalphaPL+InjectedGaussianCoefficient*np.exp(-0.5*((ga-InjectedGaussianLocation)/InjectedGaussianWidth)**2.0))*Heavi(Injectedga1,ga,Injectedga0)
# Plotting-test:    X=np.logspace(np.log10(Injectedga1),np.log10(Injectedga0),901);Y=DotNiGaussPlusPL(X);pl.loglog(X,Y)

def DotNiRelMax(ga): # A Maxwell-Juettner-distribution is used here. Without the coefficient InjectedRelMaxCoefficient and without the Heaviside-function, it is normalised to 1, that is its integral is 1. With the coefficient InjectedRelMaxCoefficient and without the Heaviside-function, it is normalised to InjectedRelMaxCoefficient, that is its integral is InjectedRelMaxCoefficient. The mean gamma of a Maxwell-Juettner-distribution is 3*InjectedTheta, i. e. the integral over ga*DotNiRelMax(ga) from 1 to infinity divided by InjectedRelMaxCoefficient is equal to 3*InjectedTheta.
    return InjectedRelMaxCoefficient*ga**2.0*np.sqrt(1.0-ga**(-2.0))*np.exp(-ga/InjectedTheta)/(InjectedTheta*ModifiedBesselSecondKind(2.0,1.0/InjectedTheta))*Heavi(Injectedga1,ga,Injectedga0)
# Plotting-test:    X=np.logspace(4,8,301);Y=DotNiRelMax(X);pl.loglog(X,Y);pl.ylim(10**(-10),10)

def DotNiHeavi(ga):     # A Heaviside-function is used here. Below Injectedga1 and above Injectedga0 it is equal to zero. Between Injectedga1 and Injectedga0 it is equal to InjectedHeaviCoefficient. At Injectedga1 and at Injectedga0 it is equal to 0.5*InjectedHeaviCoefficient.
    return InjectedHeaviCoefficient*Heavi(Injectedga1,ga,Injectedga0)
# Plotting-test:    X=np.logspace(2,5,301);Y=DotNiHeavi(X);pl.loglog(X,Y)

def DotNiZero(ga):     # A vanishing distribution is used here. 
    return 0.0

if InjectionType=='Chosen':
    UsedDotNi = DotNiGauss  # The used distribution for the injection-rate. INPUT VALUE!
elif InjectionType=='ResultsOfOldIt':
    def UsedDotNi(ga):
        return OldCurrentNElectrons(ga)*FactorNumberDensityToInjectionRate # The used distribution for the injection-rate. It is also possible to use the imported electrons' number-density, here. This has to be converted into an injection rate by multiplying FactorNumberDensityToInjectionRate. In this case, however, one has to redefine Injectedga1 and Injectedga0, too.

## Injection of primary HE-photons

# Various physically reasonable possibilities for the injected high-energy photons are given now. The distributions in xga are not extending unlimited but are cutted-off at certain values of xga. However, this is not realised via a Heaviside-function, here. Therefore it has to be considered later, if for example integration or plotting is to be done.

HEPhotonsalphaPL = -2.0    # In the case DotniPL, this is the power-law-index. In the case DotniPLSmoothlyBroken, this is the constant part of the power-law index. In the case DotniPLBroken, this is the power-law index of the first branch. INPUT VALUE!
HEPhotonsbetaPL = -0.05 # In the case DotniPLSmoothlyBroken, this is the slope/varying part of the power-law-index, i.e. the curvature of the power-law in log-log-plotting. In the case DotniPLBroken, this is the power-law index of the second branch. INPUT VALUE!
HEPhotonsPLSharpnessIndex = 2.0 # A number that determines the sharpness/smoothness of the break in a broken power-law. Only non-zero values are admitted and usually only positive values are reasonable. High values mean a narrow break and low values cause a wide break. INPUT VALUE!
HEPhotonsExpCutOff = -0.0001    # The exponent's coefficient. INPUT VALUE!
HEPhotonsPlanckTemperature = 10.0**12.0       # The temperature in units of K! INPUT VALUE!
HEPhotonsTheta = Theta(HEPhotonsPlanckTemperature)  # The thermal energy in units of electron rest-mass energy.
HEPhotonsGaussianLocation = 3.0*10**(12)/(me*(c**2)/e)  # The value of xga, where the Gaussian distribution peaks, i. e. the mean. INPUT VALUE!
HEPhotonsGaussianWidth = 0.23*HEPhotonsGaussianLocation       # The width of the Gaussian distribution. INPUT VALUE!
HEPhotonsxga0Delta = 2.0*10**16/(me*(c**2)/e)   # The location of the Delta-peak. INPUT VALUE! Assume that 20 PeV-photons are injected. These could stem from pion decay.
if InjectionType=='ResultsOfOldIt':
    HEPhotonsxga0NonDelta = OldValuesForxgammanHEPhotons[-1] # The upper cut-off of the distribution of the injection-rate. INPUT VALUE!
    HEPhotonsxga1NonDelta = OldValuesForxgammanHEPhotons[0] # The lower cut-off of the distribution. INPUT VALUE!
elif InjectionType=='Chosen':
    HEPhotonsxga0NonDelta = HEPhotonsGaussianLocation+3.0*HEPhotonsGaussianWidth#1.0*10**7 # The upper cut-off of the distribution of the injection-rate. INPUT VALUE!
    HEPhotonsxga1NonDelta = HEPhotonsGaussianLocation-3.0*HEPhotonsGaussianWidth#0.01*HEPhotonsxga0NonDelta # The lower cut-off of the distribution. INPUT VALUE!
HEPhotonsxgaBreak = 1.0*10.0**6    # The location of the break in a broken power-law. INPUT VALUE!
# The Gaussian's coefficient and the Gaussian's norm are specified, now. Consider now that the relation Norm=Coefficient*(Width*np.sqrt(2*np.pi)) holds. Thus, if Coefficient is given, Norm can be determined and if Norm is given, Coefficient can be determined.
UseHEPhotonsGaussianCoefficientManualInput = False # A boolean, that states whether a self chosen value for Coefficient (True) or a self chosen value for Norm (False) is used. INPUT VALUE!
HEPhotonsGaussianNormManualInput = 6.2*10**3 # This is the norm of a Gaussian distribution, if HEPhotonsGaussianNorm/(HEPhotonsGaussianWidth*np.sqrt(2*np.pi)) is used as its coefficient. INPUT VALUE!
HEPhotonsGaussianCoefficientManualInput = 10.0**(5.0) # A self chosen coefficient for dot n_i being a Gaussian distribution. It is only applied, if...  INPUT VALUE!
if UseHEPhotonsGaussianCoefficientManualInput == True: # ... this evaluates to true.
    HEPhotonsGaussianCoefficient = HEPhotonsGaussianCoefficientManualInput
    HEPhotonsGaussianNorm = HEPhotonsGaussianCoefficient*(HEPhotonsGaussianWidth*np.sqrt(2.0*np.pi))
elif UseHEPhotonsGaussianCoefficientManualInput == False: # Otherwise, a value, that normalises the Gaussian to HEPhotonsGaussianNorm, is applied.
    HEPhotonsGaussianNorm = HEPhotonsGaussianNormManualInput
    HEPhotonsGaussianCoefficient = HEPhotonsGaussianNorm/(HEPhotonsGaussianWidth*np.sqrt(2.0*np.pi))
HEPhotonsExpCoefficient = 10.0**(10.0) # The coefficient for dot n_i being an exponential-function. INPUT VALUE!
HEPhotonsPLCoefficient = 4*10.0**(10.0) # The coefficient for dot n_i being a power-law. INPUT VALUE!
DotniDelta = 0.001  # A Dirac-Delta-function can be used as the population of injected HE-photons and this quantity is the Dirac-Delta-function's coefficient. In this case only the coefficient for dot n_i is implemented. The delta(xga-xga0)-part is omitted as it drops away anyway via the substitution in the integral. INPUT VALUE!    

# The following quantities correspond to dot n_i(x_ga). xga is the dimensionless energy of the photons. According to my analysis, it is a injection-rate, i. e. a spectral number-density per time-interval. Thus, it is in units of 1/(s*m^3).

def DotniExp(xga):  # An exponential-function is used here.
    return HEPhotonsExpCoefficient*np.exp(-xga/HEPhotonsExpCutOff)
# Plotting-test:    X=np.logspace(3,5,31);Y=DotniExp(X);pl.plot(X,Y);pl.yscale('log')

def DotniPL(xga):  # A power-law is used here.
    return HEPhotonsPLCoefficient*xga**(HEPhotonsalphaPL)
# Plotting-test:    X=np.logspace(np.log10(HEPhotonsxga1NonDelta),np.log10(HEPhotonsxga0NonDelta),31);Y=DotniPL(X);pl.loglog(X,Y)

def DotniPLExpCutOff(xga): # A power-law with an exponential cut-off. A Heaviside-function is not needed here, as it is already implemented in DotniPL.
    return DotniPL(xga)*np.exp(-xga/HEPhotonsExpCutOff)
# Plotting-test:    X=np.logspace(-1,5,31);Y=DotniPLExpCutOff(X);pl.loglog(X,Y)

def DotniPLSmoothlyBroken(xga):  # A power-law with logarithmically changing spectral index is used here. At xga=HEPhotonsxga1NonDelta it adopts the value of the usual power-law. At higher values of xga, the spectral index changes according to the choice of HEPhotonsbetaPL. If HEPhotonsbetaPL < 0 (> 0), the power-law gets softer (harder) with increasing xga. The choice of the base of the logarithm doesnt matter, as it can be absorbed in HEPhotonsbetaPL. This function is similar to what is sometimes called a LogParabola.
    return HEPhotonsPLCoefficient*(xga/HEPhotonsxga1NonDelta)**(HEPhotonsalphaPL+HEPhotonsbetaPL*np.log10(xga/HEPhotonsxga1NonDelta))
# Plotting-test: X=np.logspace(np.log10(HEPhotonsxga1NonDelta),np.log10(HEPhotonsxga0NonDelta),31);Y=DotniPLSmoothlyBroken(X);pl.loglog(X,Y)

def DotniEplogpar(xga):
    '''This was used by Josefa and is similar to a logparabola'''
    return 10**(HEPhotonsbetaPL*(np.log10((xga/HEPhotonsxga1NonDelta)/HEPhotonsxgaBreak))**2)/(xga/HEPhotonsxga1NonDelta)**2
# Plotting-test: X=np.logspace(np.log10(HEPhotonsxga1NonDelta),np.log10(HEPhotonsxga0NonDelta),31);Y=DotniEplogpar(X);pl.loglog(X,Y)

def DotniPLBroken(xga):  # A broken power-law, i.e. two power-laws, whose spectral indices blend into each other within a given range, that is specified by HEPhotonsxgaBreak and HEPhotonsPLSharpnessIndex. This is from 2019 Meyer, Scargle & Blandford.
    return HEPhotonsPLCoefficient*((xga/HEPhotonsxga1NonDelta)**HEPhotonsalphaPL)*(1.0+(xga/HEPhotonsxgaBreak)**HEPhotonsPLSharpnessIndex)**((HEPhotonsbetaPL-HEPhotonsalphaPL)/HEPhotonsPLSharpnessIndex)
# Plotting-test: X=np.logspace(np.log10(HEPhotonsxga1NonDelta),np.log10(HEPhotonsxga0NonDelta),310);Y=DotniPLBroken(X);pl.loglog(X,Y)

def DotniPLBroken2(xga):  # A broken power-law, i.e. two power-laws, whose spectral indices blend into each other, that is specified by HEPhotonsxgaBreak. This is based on 2017MNRAS.469..255G. It has however the caveat that the slope of the first power-law is always that power-law index which is bigger and the slope of the second power-law is always that power-law index which is smaller.
    return HEPhotonsPLCoefficient*((xga/HEPhotonsxgaBreak)**(HEPhotonsalphaPL))/(1.0+(xga/HEPhotonsxgaBreak)**(HEPhotonsalphaPL-HEPhotonsbetaPL))
# Plotting-test: X=np.logspace(np.log10(HEPhotonsxga1NonDelta),np.log10(HEPhotonsxga0NonDelta),310);Y=DotniPLBroken2(X);pl.loglog(X,Y)

def DotniPlanck(xga):  # A Planck-spectrum is used here, even though this seems unphysically. For reasonable temperatures, no significant population in the HE-regime is caused.  
    return (8.0*np.pi*me**3.0*c**3.0*xga**2.0) / ( h**3.0 * (np.exp(xga*me*c**2.0/(kB*HEPhotonsPlanckTemperature)) - 1.0) )
# Plotting-test:    X=np.logspace(3,5,100);Y=DotniPlanck(X);pl.loglog(X,Y)

def DotniGauss(xga):    # A Gaussian distribution is used here. It is not normalised. 
    return HEPhotonsGaussianCoefficient*np.exp(-0.5*((xga-HEPhotonsGaussianLocation)/HEPhotonsGaussianWidth)**2)
# Plotting-test:    X=np.logspace(3,5,61);Y=DotniGauss(X);pl.loglog(X,Y)

def DotniZero(xga): # Just a vanishing distribution. 
    return 0.0

if InjectionType=='Chosen':
    UsedDotni = DotniZero # The used distribution for the injection-rate. INPUT VALUE!
elif InjectionType=='ResultsOfOldIt':
    def UsedDotni(xga):
        return OldnHEPsInterpolated(xga)*FactorNumberDensityToInjectionRate # The used distribution for the injection-rate. One can also use the imported HEPs' number-density, here. This has to be converted into an injection rate by multiplying with FactorNumberDensityToInjectionRate. Then, however, one has to redefine HEPhotonsxga0NonDelta and HEPhotonsxga1NonDelta, too.

if UsedDotni == DotniDelta:
    HEPhotonsxga0 = HEPhotonsxga0Delta       # This makes the following code easier
    HEPhotonsxga1 = HEPhotonsxga0Delta       # This makes the following code easier
    NameOfUsedDotni = 'DotniDelta'
else:
    HEPhotonsxga0 = HEPhotonsxga0NonDelta    # This makes the following code easier
    HEPhotonsxga1 = HEPhotonsxga1NonDelta    # This makes the following code easier
    if isinstance(UsedDotni, types.FunctionType):
        NameOfUsedDotni = UsedDotni.__name__

# Determine the power-density of the injected photons, which is in units of 1/(m^3*s):
if UsedDotni != DotniDelta:
    def HEPhotonsPowerDensitySpectralInxga(xga,Dotni):
        return xga*Dotni(xga)

if UsedDotni == DotniDelta:
    def HEPhotonsPowerDensity(Dotni,RelativeError=1): # The second argument is only for syntactical reasons.
        return HEPhotonsxga0Delta*DotniDelta
else:
    def HEPhotonsPowerDensity(Dotni,RelativeError=integrate.quad.__defaults__[3]):
        if InjectionType=='Chosen':
            LowestIntegrationBorder = HEPhotonsxga1
            BiggestIntegrationBorder = HEPhotonsxga0
            if LowestIntegrationBorder >= BiggestIntegrationBorder:
                return 0.0
            else:
                NumberOfBorders = max(6,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
                ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
                HEPhotonsPowerDensityResult = 0.0
                for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                    HEPhotonsPowerDensityResult += integrate.quad(HEPhotonsPowerDensitySpectralInxga, LeftBorder, RightBorder, args=(Dotni), epsrel=RelativeError)[0]
                return HEPhotonsPowerDensityResult
        elif InjectionType=='ResultsOfOldIt':
            ListOfIntegrationBorders = OldValuesForxgammanHEPhotons
            HEPhotonsPowerDensityResult = 0.0
            for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                HEPhotonsPowerDensityResult += integrate.quad(HEPhotonsPowerDensitySpectralInxga, LeftBorder, RightBorder, args=(Dotni), epsrel=RelativeError)[0]
            return HEPhotonsPowerDensityResult


## Soft external background photon-fields 

# Physically possible soft, external photon-fields are implemented now. Previously these photons have been denoted by SBPs (soft background photons) or LEPs (low-energy photons). With the inclusion of synchrotron-photons, two populations of low-energy photons are available, the synchrotron-photons (SyncPs) and the soft external photons (ExPs). The umbrella term for both these will now be lowenergy photons (LEPs).
# The spectral number-densities of the ExPs in units 1/m^3 are considered now.
# In version 1, the cut-off of the background photon-fields was realised via Heaviside-functions in the spectral number-densities. During this, it was integrated over the whole range of photon energy x, as it is to be done according to eqauation A1. Now, as one notices that the contribution to the integral of the part in the cutted-off region of the spectral number-density is vanishing anyway, one can neglect this part in the integral, that is one can shorten the integration-range to the non-vanishing part of the spectral number-density. Then again, one can drop the Heaviside-functions.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n\n---------------- Specifying the soft external background-photon-field -------------------\n')

ImportedLineDataFile = 'Emission lines for Mrk 501 case.csv' # The file, that will be imported as input emission line data. INPUT VALUE!
LineData=ImportFileToString(ImportedLineDataFile) # Import the file, containing the table, including the wavelengths and the relative contributions of lines that are used as the soft background photons.
LineWavelengthnm = np.asarray([]) # An array of the lines' wavelengths in units of nm.
LineFluxContribution = np.asarray([]) # An array of the lines' relative contributions, corresponding to the respective wavelength.
LineName = []
for i in range(1,len(LineData)):
    LineWavelengthnm = np.append(LineWavelengthnm,eval(LineData[i].replace("\n", "").split(';')[0]))
    LineFluxContribution = np.append(LineFluxContribution,eval(LineData[i].replace("\n", "").split(';')[1]))
    LineName = np.append(LineName,LineData[i].replace("\n", "").split(';')[2].split(' #')[0])

alphaPL = -1.7    # The power-law-index. INPUT VALUE!
alphaExp = -200.0    # The exponent's coefficient. INPUT VALUE!
PlanckTemperature = 0.35*10.0**4.0       # The temperature in units of K! INPUT VALUE!
PlanckTheta = Theta(PlanckTemperature)  # The thermal energy in units of electron rest-mass energy 
x0NonDelta = 12*PlanckTheta   # The upper cut-off of the spectral number-density. INPUT VALUE!
x0Delta = h/(121.5*10**(-9)*me*c)   # The situation of the Delta-peak. E. g. the Lyman-alpha energy (121.5 nm) or the Balmer-alpha energy (656.3 nm). INPUT VALUE!
x0MultiDelta = h/(LineWavelengthnm*10**(-9)*me*c) # The situations of a set of Delta-peaks, which has to be given as a numpy-array. INPUT VALUE!
x1 = 0.02*PlanckTheta    # The lower cut-off of the spectral number-density. INPUT VALUE!
ExpCoefficient = 10.0**30.0 # The coefficient for n_0 being an exponential-function. INPUT VALUE!
PLCoefficient = 10.0**26.0 # The coefficient for n_0 being a power-law. INPUT VALUE!
UseDeltaCoefficientManualInput = True # A boolean, that states whether a self chosen value (True) or a value (False), which was predefined by Zdziarski, is to be used as DeltaCoefficient. INPUT VALUE!
DeltaCoefficientManualInput = 1.0*10.0**20.0 # This is an auxiliary quantity for the following function. It describes DeltaCoefficient in the cases, where no special normalisation is used. INPUT VALUE! 
# Now, define the value of DeltaCoefficient, which is the norm of n0 being a Dirac-Delta-function. The reason for these normalisations could be that it ensures that the cascade remains saturated and linear.
if UsedDotNi==DotNiZero and UsedDotni==DotniDelta and UseDeltaCoefficientManualInput == False:
    DeltaCoefficient = DotniDelta/(sigmaT*c) # Here, the normalisation K/(sigmaT*c) according to (page 794 of) Zdziarski is used, where K is the norm of either the injected HE-photon-distribution...
elif UsedDotNi==DotNiGauss and UsedDotni==DotniZero and UseDeltaCoefficientManualInput == False:
    DeltaCoefficient = HEPhotonsGaussianNorm/(sigmaT*c) # ... or of the injected particle-distribution.
else:
    DeltaCoefficient = DeltaCoefficientManualInput
CoefficientKLines = 9.7*10**12 # This is a measure of the total number-density of the emission line photons. It is denoted by K_lines in my documentations. It is in units of 1/m^3. INPUT VALUE!
MultiDeltaCoefficient = CoefficientKLines*LineFluxContribution/x0MultiDelta # The set of coefficients of the Dirac-Delta-functions, that are located at the values of x0MultiDelta with the corresponding index. These coefficients are measured in 1/m^3. LineFluxContribution gives the relative contributions to the flux. Dividing this by the energies of the respective line gives the relative contributions to the number-density. The factor CoefficientKLines gives the normalisation of the soft background photon-field.
AFsADAFContribution = 0.5*ComplementaryErrorFunction(2.0*(np.log10(Dotm)+1.3)) # This was parametrised based on 1997ApJ...489..865E Esin, McClintock, Narayan. It describes the ADAF-contribution to the entire accretion flow's photon-field. At Dotm below 0.001, it is =1, above Dotm=1, it is =0 and in between it decreases from 1 to 0, with being =0.5 at Dotm=10^(-1.3). It is dimensionless.
AFsPlanckContribution = 1.0-AFsADAFContribution # Describes an accretion flow's Planck-contribution.
ADAFx1 = 10**7*h/(me*c**2) # The energy of the lowest border (in dimensionless units), where the ADAF is assumed to be cut-off. 10^7 Hz is motivated by the frequency, at which radio waves are absorbed by ionised media.
ADAFxMin = ApplyADAFIntelligentStorage('ADAFnuMin')*h/(me*c**2) # The value of the dimensionless energy x, where the ADAF has a kink, due to the change from ADAFSpectralLuminosityCyclosynchrotronBelownuMin to ADAFSpectralLuminosityCyclosynchrotronOf25.
ADAFxPeak = ApplyADAFIntelligentStorage('ADAFnuPeak')*h/(me*c**2) # The value of the dimensionless energy x, where the ADAF has a kink, due to the change from ADAFSpectralLuminosityCyclosynchrotronOf25 to ADAFSpectralLuminosityComptonised.
ADAFxComptonisedCutOff = 3*kB*ADAFTemperatureElectron/(me*c**2) # The value of the dimensionless energy x, where the ADAF has a cut-off, due to the cut-off of ADAFSpectralLuminosityComptonised.
ADAFxBremsstrahlungCutOff = 300*kB*ADAFTemperatureElectron/(me*c**2) # Mathematically, the bremsstrahlung-spectrum is extended up to infinity. As an integrand, it causes problems, however, if it becomes too small. Therefore, this value of x is assumed to be a cut-off of the bremsstrahlung-contribution. It was found in the computation of ADAFTotalEnergyDensity such that the result is very similar to the analytical one.

if ADAFxBremsstrahlungCutOff <= 0.1: # The case satisfying Zdz's assumption x<<1.
    ADAFx0 = ADAFxBremsstrahlungCutOff
else: # The case where the ADAF reaches energies too high with respect to the assumption.
    ADAFx0 = 0.1
    if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
        print("\nAttention: The ADAF's upper cut-off had to be reduced to satisfy the assumption x<<1.")

# x is the dimensionless energy of all soft photons (LEPs). 
def n0Exp(x):  # A cutted-off exponential-function is used here.
    return ExpCoefficient*np.exp(alphaExp*x)
# Plotting-test:    X=np.logspace(-6,0,100);Y=n0Exp(X);pl.loglog(X,Y);pl.xscale('log')

def n0PL(x):  # A possible soft photon-field. This gives the spectral number-density of the soft photons in units 1/m^3. A cutted-off power-law is used here.
    return PLCoefficient*x**(alphaPL)
# Plotting-test:    X=np.logspace(-6,0,100);Y=n0PL(X);pl.loglog(X,Y)

n0PlanckCoefficient = (8.0*np.pi*me**3.0*c**3.0)/(h**3.0) # This can make the evaluation of the following function faster.
def n0Planck(x):  # A possible soft photon-field. This gives the spectral number-density of the soft photons in units 1/m^3. A cutted-off Planck-spectrum is used here. 
    return (n0PlanckCoefficient*x**2.0) / (np.exp(x/PlanckTheta) - 1.0)
# Plotting-test:   X=np.logspace(np.log10(x1),np.log10(x0NonDelta),100);Y=n0Planck(X);pl.loglog(X,Y)

def n0Delta(x):  # A possible soft photon-field. This gives the number-density of the soft photons in units 1/m^3. A Dirac-Delta-function is used here. Only the coefficient is implemented. The delta(x-x0)-part is omitted as it drops away anyway via the substitution in the integral. 
    return DeltaCoefficient

def n0MultiDelta(x):  # A possible soft photon-field. This gives the number-density of the soft photons in units 1/m^3. A set of Dirac-Delta-functions is used here. Each Dirac-Delta-peak has its own location and coefficient. Only the coefficient is implemented. The delta(x-x0)-part is omitted as it drops away anyway via the substitution in the integral. 
    return float(MultiDeltaCoefficient[x==x0MultiDelta]) # The part inside of the squared brackets yields a boolean mask, which if applied to MultiDeltaCoefficient determines that coefficient that is corresponding to x. If x adopts a value that is not included in x0MultiDelta, an error is raised.

def n0AF(x): # The photon-field of an accretion flow.
    return AFsPlanckContribution*n0Planck(x) + AFsADAFContribution*ADAFNumberDensitySpectralInx(x,ADAFSpectralLuminosity,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    def Evaluaten0AF():
        ValuesForx = np.logspace(np.log10(ADAFx1),np.log10(ADAFxBremsstrahlungCutOff),100)
        ValuesForn0AF = np.asarray([n0AF(i) for i in ValuesForx])
        pl.figure(figsize=(16, 12), num="Spectral number-density of an accretion flow versus dimensionless energy x")
        pl.title('Spectral number-density of an accretion flow versus dimensionless energy x', fontsize=24)
        pl.loglog(ValuesForx,ValuesForn0AF, label="$\dot m = 0.0001$")
        pl.xlabel('Dimensionless energy $x$', fontsize=24)
        pl.ylabel('Spectral number-density in $m^{-3}$', fontsize=24)
        ax = pl.gca()
        pl.legend(loc="best", fontsize=22)
        ax.xaxis.set_tick_params(labelsize=24)
        ax.yaxis.set_tick_params(labelsize=24)

Usedn0 = n0MultiDelta  # The used distribution for the soft photon background. INPUT VALUE!

if Usedn0 == n0Delta:
    x0 = x0Delta    # This makes the following code easier.
elif Usedn0 == n0MultiDelta:
    x0 = max(x0MultiDelta)    # This makes the following code easier.
elif Usedn0==n0Exp or Usedn0==n0PL or Usedn0==n0Planck:
    x0 = x0NonDelta        # This makes the following code easier.
elif Usedn0==n0AF:
    x0 = ADAFx0        # This makes the following code easier.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n\n------------------------- Check for validity of the input values ------------------------\n')

    if ADAFx1 >= ADAFxMin or ADAFxMin >= ADAFxPeak or ADAFxPeak >= ADAFxComptonisedCutOff: # This is the within the scope of the current parameters maximally possible value of epsilonCMOfEquation2_2. If it is >1 pair-production sets in.
        raise InvalidOrderOfADAFParametersError()
    else: # The case ADAFx1 < ADAFxMin < ADAFxPeak < ADAFxComptonisedCutOff.
        print("The ADAF's borders are well-ordered.")
    
    if UsedDotNi != DotNiZero:
        if Injectedga1 < 10.0: # The injected particles have to be relativistic.
            raise InvalidSetOfParametersError('Invalid energy of injected particles.')
        else:
            print('The injected particles are relativistic.')
        if Injectedga0*x0 < 1: # For the Klein-Nishina-regime to be satisfied, it has to be ga*x>=1.
            raise InvalidSetOfParametersError('Klein-Nishina-regime not satisfied.')
        else:
            print('Klein-Nishina-limit is given.')
    
    if UsedDotni != DotniZero:
        if HEPhotonsxga1 < 10.0: # The injected photons' energy has to be above 10*me*c^2.
            raise InvalidSetOfParametersError('Invalid energy of injected photons.')
        else:
            print("The injected photons' energy satisfies xga>>1.")
        if HEPhotonsxga0*x0 < 1: # For pair-production to happen, it has to be xga*x>=1.
            raise InvalidSetOfParametersError('Pair-production threshold not satisfied.')
        else:
            print('The pair-production threshold is satisfied.')


## Total number-density of soft background-photons

# It is to be integrated over x from x_1 to x_0. The resulting quantity is the total number-density of all the soft photons. Its unit is again 1/m^3. The total number-density of an ADAF was already defined in part1.
if Usedn0 == n0Delta:
    def TotalNumberDensity(n0):   # Using the substitution-feature of the Dirac-Delta-function. 
        return DeltaCoefficient
elif Usedn0 == n0MultiDelta:
    def TotalNumberDensity(n0):   # Using the substitution-feature for the set Dirac-Delta-functions. 
        return np.sum(MultiDeltaCoefficient)
elif Usedn0==n0Exp or Usedn0==n0PL or Usedn0==n0Planck:
    def TotalNumberDensity(n0):   # Ordinary Integration. 
        IntegralOfn0 = integrate.quad(n0, x1, x0NonDelta)[0]
        return IntegralOfn0


## Energy-density of soft background-photons

def IntegrandOfDenominatorOfEq8(x,n0):
    '''The integrand of the denominator of equation 8 is x*n0(x) according to Zdz., which essentially is the spectral energy-density of the soft photon-field. It is in units of 1/m^3.'''
    return x*n0(x)

# It is to be integrated along x from 0 to x_0, that is effectively from x_1 to x_0. The resulting quantity is the total energy-density of all the external, soft photons. Its unit is again 1/m^3. For the case of an ADAF, the energy-density was again already defined in part1.
if Usedn0 == n0Delta:
    def IntegralOfDenominatorOfEq8(n0):   # Using the substitution-feature of the Dirac-Delta-function. 
        return IntegrandOfDenominatorOfEq8(x0Delta,n0)
elif Usedn0 == n0MultiDelta:
    def IntegralOfDenominatorOfEq8(n0):   # Using the substitution-feature for the set Dirac-Delta-functions. 
        IntegralOfDenominatorOfEq8Sum = 0.0
        for i in x0MultiDelta:
            IntegralOfDenominatorOfEq8Sum += IntegrandOfDenominatorOfEq8(i,n0)
        return IntegralOfDenominatorOfEq8Sum
elif Usedn0==n0Exp or Usedn0==n0PL or Usedn0==n0Planck:
    def IntegralOfDenominatorOfEq8(n0):   # Ordinary Integration. 
        IntegralOfDenominatorOfEq8Result = integrate.quad(IntegrandOfDenominatorOfEq8, x1, x0NonDelta, args=(n0))[0]
        return IntegralOfDenominatorOfEq8Result

ExPsEnergyDensity = IntegralOfDenominatorOfEq8(Usedn0)*me*c**2 # The total energy-density in units of J/m^3 of the currently used external, soft photon-field (ExPs).

if Usedn0 == n0MultiDelta:
    EnergyDensitySingleDeltas = np.asarray([IntegrandOfDenominatorOfEq8(i,Usedn0) for i in x0MultiDelta]) # This is an array which contains the energy-density values for each single Dirac-Delta-function-peak, in units of 1/m^3.
    EnergyDensitySingleDeltaContributions = EnergyDensitySingleDeltas/IntegralOfDenominatorOfEq8(Usedn0) # This dimensionless array gives the ratios (relative contributions), with which each single line contributes to the total energy-density.

def PlanckEnergyDensitySpectralInx(x):
    '''This might be the same as IntegrandOfDenominatorOfEq8 in the case n0Planck. It is not obvious why it was defined twice.'''
    return x*n0Planck(x)

PlanckEnergyDensity = integrate.quad(PlanckEnergyDensitySpectralInx, x1, x0NonDelta)[0]*me*c**2 # The total energy-density of the Planck-distribution in units of J/m^3.

def AFEnergyDensitySpectralInx(x):
    return x*n0AF(x)

def AFEnergyDensity(LowestIntegrationBorder,BiggestIntegrationBorder):
    if LowestIntegrationBorder >= BiggestIntegrationBorder:
        return 0.0
    else:
        NumberOfBorders = max(6,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        for PotentialBorder in [ADAFxMin,ADAFxPeak,ADAFxComptonisedCutOff]:
            if LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder and (PotentialBorder not in ListOfIntegrationBorders):
                ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
        ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
        #print('List of integration borders of AFEnergyDensity', ListOfIntegrationBorders)
        AFEnergyDensityResult = 0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            AFEnergyDensityResult += integrate.quad(AFEnergyDensitySpectralInx, LeftBorder, RightBorder)[0]
        return AFEnergyDensityResult*me*c**2 # The unit of this is J/m^3.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n\n------------------------------------ Accretion flow -------------------------------------\n')
    
    print('Energy-density of pure ADAF-emission (numerically integrated):   ', ADAFTotalEnergyDensityNum, 'J/m^3')
    print('Energy-density of pure ADAF-emission (analytically integrated):  ', ADAFTotalEnergyDensityAna, 'J/m^3')
    print('Luminosity of pure ADAF-emission (numerically integrated):   ', c*ADAFSurface(ADAFrOuter,ADAFrInner,M9)*ADAFTotalEnergyDensityNum, 'J/m^3')
    print('Luminosity of pure ADAF-emission (analytically integrated):  ', c*ADAFSurface(ADAFrOuter,ADAFrInner,M9)*ADAFTotalEnergyDensityAna, 'J/m^3')
        
    print('                                  weighted Planck-contribution:  ', AFsPlanckContribution*PlanckEnergyDensity, 'J/m^3')
    print('                                    weighted ADAF-contribution:  ', AFsADAFContribution*ADAFTotalEnergyDensityAna, 'J/m^3')
    print('Energy-density of AF-emission (weighted sum ADAF + Planck):      ', AFsPlanckContribution*PlanckEnergyDensity + AFsADAFContribution*ADAFTotalEnergyDensityAna, 'J/m^3')
    
    print('Energy-density of AF-emission (weighted, completely integrated): ', AFEnergyDensity(ADAFx1,ADAFxBremsstrahlungCutOff), 'J/m^3')
    print('Fraction of energy-density of AF-emission that is omitted by resetting the upper cut-off:', AFEnergyDensity(0.1,ADAFxBremsstrahlungCutOff)/AFEnergyDensity(ADAFx1,ADAFxBremsstrahlungCutOff))

    print('\n\n------------------------ Soft background-photons and BLR-photons ------------------------\n')
    
    print('Energy-density of pure Planck-distribution (integrated):             ', PlanckEnergyDensity, 'J/m^3')
    print('Energy-density of pure Planck-distribution (Stefan-Boltzmann-law)):  ', BBTotalEnergyDensity(PlanckTemperature), 'J/m^3\n')
    
    if Usedn0==n0Exp or Usedn0==n0PL or Usedn0==n0Planck or Usedn0 == n0Delta or Usedn0 == n0MultiDelta:
        print('Test for IntegralOfDenominatorOfEq8:\n    Total energy-density of soft external background-photons:        ', ExPsEnergyDensity, 'J/m^3')
    print('Total number-density of soft external background-photons:             %g 1/m^3' % TotalNumberDensity(Usedn0))
    if Usedn0 == n0MultiDelta:
        print('Total number-density of external H I Lyman-alpha-photons:             %g 1/m^3\n' % MultiDeltaCoefficient[3])


## The external photons' luminosity

# Consider a typical homogeneous broad-line region:
if Usedn0 == n0Delta:
    if InjectionType=='Chosen':
        ExPsLuminosityHomog = TotalNumberDensity(Usedn0)*(4.0*np.pi*ExPsOuterRadius**2*c*x0Delta*me*c**2) # An estimate for the luminosity of external photons in units of W.
    elif InjectionType=='ResultsOfOldIt':
        ExPsLuminosityHomog = TotalNumberDensity(Usedn0)*(4.0*np.pi*ExPsVeryOuterRadius**2*c*x0Delta*me*c**2) # The same estimate for a bigger region.
elif Usedn0 == n0MultiDelta:
    if InjectionType=='Chosen':
        # BLREnergyDensitySingleDeltas = EnergyDensitySingleDeltaContributions*BLRLuminosity/(4.0*np.pi*ExPsOuterRadius**2*c*me*c**2) # This array, which is in units of 1/m^3, yields the absolute contributions of each single line to the total energy-density.
        # BLRn = np.sum(BLREnergyDensitySingleDeltas/x0MultiDelta) # The total number-density of BLR-photons in units of 1/m^3.
        ExPsLuminosityHomog = ExPsEnergyDensity*FactorEnergyDensityToLuminosity('ExPsOuterRadius')[1] # The approximate luminosity of the entire homogeneous external photon-field in units of W.
        # LyAlphan = LyAlphaLuminosity/(4.0*np.pi*ExPsOuterRadius**2*c*x0MultiDelta[3]*me*c**2) # An estimate for the total number-density of Lyman-alpha-photons in units of 1/m^3.
        HILyAlphaLuminosityHomog = EnergyDensitySingleDeltas[3]*me*c**2*FactorEnergyDensityToLuminosity('ExPsOuterRadius')[1] # The approximate luminosity of the homogeneous H I Lyman-alpha-photon-field in units of W.
        Around135nmLuminosityHomog = np.sum(EnergyDensitySingleDeltas[7:10])*me*c**2*(4.0*np.pi*ExPsOuterRadius**2*c) # The approximate luminosity of the homogeneous photon-field around 135 nm in units of W. Lines from 135-135/2 nm till 135+135/2 nm are included.
        UVLuminosityHomog = np.sum(EnergyDensitySingleDeltas[3:])*me*c**2*(4.0*np.pi*ExPsOuterRadius**2*c) # The approximate luminosity of the homogeneous UV-photon-field in units of W. Lines above 10 nm are included.
    elif InjectionType=='ResultsOfOldIt':
        ExPsLuminosityHomog = ExPsEnergyDensity*FactorEnergyDensityToLuminosity('ExPsVeryOuterRadius')[1] # The same, just for a bigger region.
        HILyAlphaLuminosityHomog = EnergyDensitySingleDeltas[3]*me*c**2*FactorEnergyDensityToLuminosity('ExPsVeryOuterRadius')[1] # The same, just for a bigger region.
        Around135nmLuminosityHomog = np.sum(EnergyDensitySingleDeltas[7:10])*me*c**2*(4.0*np.pi*ExPsVeryOuterRadius**2*c) # The same, just for a bigger region.
        UVLuminosityHomog = np.sum(EnergyDensitySingleDeltas[3:])*me*c**2*(4.0*np.pi*ExPsVeryOuterRadius**2*c) # The same, just for a bigger region.

def RadiusLuminosityRelationKaspi2007(Luminosity):
    '''This determines the radius of the BLR in units of m from the 135 nm-continuum-luminosity in units of W. The radius is the C IV time lag radius. This is from 2007ApJ...659..997K.'''
    return 10**15*0.20*(Luminosity/10**37)**0.55

def RadiusLuminosityRelationKaspi2021(Luminosity):
    '''This determines the radius of the BLR in units of m from the 135 nm-continuum-luminosity in units of W. The radius is the C IV time lag radius. This is from 2021ApJ...915..129K.'''
    return 10**15*0.25*(Luminosity/10**37)**0.45

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print("\n------------------------ Soft background-photons' luminosity ------------------------\n")
    
    if InjectionType=='Chosen':
        print('Assuming homogeneous external photon-field of size %3.4g m = %3.4g pc = %3.4g r_S:' % (ExPsOuterRadius, mToMpc(ExPsOuterRadius)*1000000, ExPsOuterRadius/rSchwarzschild(M9)))
    elif InjectionType=='ResultsOfOldIt':
        print('Assuming homogeneous external photon-field of size %3.4g m = %3.4g pc = %3.4g r_S:' % (ExPsVeryOuterRadius, mToMpc(ExPsVeryOuterRadius)*1000000, ExPsVeryOuterRadius/rSchwarzschild(M9)))
    if Usedn0 == n0Delta or Usedn0 == n0MultiDelta:
        print('    Corresponding total luminosity:                      %g W' % ExPsLuminosityHomog)
        print('        Corresponding C IV radius:                       %g m' % RadiusLuminosityRelationKaspi2021(ExPsLuminosityHomog))
    if Usedn0 == n0MultiDelta:
        print('    Corresponding H I Lyman-alpha-luminosity:            %g W' % HILyAlphaLuminosityHomog)
        print('        Corresponding C IV radius:                       %g m' % RadiusLuminosityRelationKaspi2021(HILyAlphaLuminosityHomog))
        print('    Corresponding 135 nm-luminosity:                     %g W' % Around135nmLuminosityHomog)
        print('        Corresponding C IV radius:                       %g m' % RadiusLuminosityRelationKaspi2021(Around135nmLuminosityHomog))
        print('    Corresponding UV-luminosity:                         %g W' % UVLuminosityHomog)
        print('        Corresponding C IV radius:                       %3.4g m = %3.4g r_S' % (RadiusLuminosityRelationKaspi2021(UVLuminosityHomog),RadiusLuminosityRelationKaspi2021(UVLuminosityHomog)/rSchwarzschild(M9)))

# Consider a number NClouds of ExPs-emitting clouds of radial size RClouds:
NClouds=10 # The number of ionised gas clouds.
RClouds = rSchwarzschildTom(0.061,M9) # The radial size of the clouds in units of m.
HILyAlphaLuminosityOneCloud = FactorEnergyDensityToLuminosity('RClouds')[1]*EnergyDensitySingleDeltas[3]*me*c**2 # The H I Lyman-alpha-luminosity of one single cloud of ExPs in units of W.
ExPsLuminosityOneCloud = FactorEnergyDensityToLuminosity('RClouds')[1]*ExPsEnergyDensity # As ExPsEnergyDensity is in J/m^3, this is in W.
HILyAlphaLuminosityAllClouds = NClouds*HILyAlphaLuminosityOneCloud # The H I Lyman-alpha-luminosity of all the clouds in units of W.
ExPsLuminosityAllClouds = NClouds*ExPsLuminosityOneCloud # The total ExPs' luminosity of all the clouds in W.
if "MainProcess" in multiprocessing.current_process().name and (Usedn0 == n0Delta or Usedn0 == n0MultiDelta):
    print('\nAssuming cloudy external photon-field with %s clouds of radial size %3.4g m = %3.4g pc = %3.4g r_S:' % (NClouds, RClouds, mToMpc(RClouds)*1000000, RClouds/rSchwarzschild(M9)))
    print('    Corresponding total luminosity per cloud:            %3.4g W' % ExPsLuminosityOneCloud)
    print('    Corresponding H I Lyman-alpha-luminosity per cloud:  %3.4g W' % HILyAlphaLuminosityOneCloud)
    print('    Corresponding total luminosity:                      %3.4g W' % ExPsLuminosityAllClouds)
    print('    Corresponding H I Lyman-alpha-luminosity:            %3.4g W' % HILyAlphaLuminosityAllClouds)

HILyAlphaLuminosityOneCloudInteractionRadius = FactorEnergyDensityToLuminosity('InteractionRadius')[1]*EnergyDensitySingleDeltas[3]*me*c**2 # The H I Lyman-alpha-luminosity in units of W of one single cloud of ExPs with radius InteractionRadius.
ExPsLuminosityOneCloudInteractionRadius = FactorEnergyDensityToLuminosity('InteractionRadius')[1]*ExPsEnergyDensity # As ExPsEnergyDensity is in J/m^3, this is in W. Again, assuming R=InteractionRadius.
if "MainProcess" in multiprocessing.current_process().name and (Usedn0 == n0Delta or Usedn0 == n0MultiDelta):
    print('\nAssuming cloudy external photon-field with clouds of radial size R=InteractionRadius:')
    print('    Corresponding total luminosity per cloud:            %3.4g W' % ExPsLuminosityOneCloudInteractionRadius)
    print('    Corresponding H I Lyman-alpha-luminosity per cloud:  %3.4g W' % HILyAlphaLuminosityOneCloudInteractionRadius)


## Energy-loss-rate and cooling time-scales

# The energy-loss-rate in the extreme Klein-Nishina-limit is computed now, with equation 2.57 of B&G1970. It is to be integrated along the non-vanishing range of SpectralNumberDensity, that is effectively from x_1 to x_0. The resulting quantity is the total energy-loss-rate of an electron on the external, soft photons. Its unit is again J/s. For the case of an ADAF, the energy-loss-rate is determined afterwards.
if Usedn0 == n0Delta:
    def EnergyLossRateICKN(gamma,n0):   # Using the substitution-feature of the Dirac-Delta-function.
        return 3.0*sigmaT*me*c**3*IntegrandOfEnergyLossRateICKN(x0Delta,gamma,n0)/8
elif Usedn0 == n0MultiDelta:
    def EnergyLossRateICKN(gamma,n0):   # Using the substitution-feature for the set Dirac-Delta-functions.
        EnergyLossRateICKNSum = 0.0
        for i in x0MultiDelta:
            EnergyLossRateICKNSum += IntegrandOfEnergyLossRateICKN(i,gamma,n0)
        return 3.0*sigmaT*me*c**3*EnergyLossRateICKNSum/8
elif Usedn0==n0Exp or Usedn0==n0PL or Usedn0==n0Planck:
    def EnergyLossRateICKN(gamma,n0):   # Ordinary Integration. 
        EnergyLossRateICKNResult = integrate.quad(IntegrandOfEnergyLossRateICKN, x1, x0NonDelta, args=(gamma,n0))[0]
        return 3.0*sigmaT*me*c**3*EnergyLossRateICKNResult/8

def CoolingTimeScale(gamma,EnergyLossRate,EnergyDensity):
    '''This is the cooling time scale of a process whose energy-loss-rate is EnergyLossRate. It is in units of s.'''
    return gamma*me*c**2/EnergyLossRate(gamma,EnergyDensity)

CoolingLengthICThomson = CoolingTimeScale(0.01/x0,EnergyLossRateICThomson,ExPsEnergyDensity)*c # The cooling length of electrons in the Thomson-range in units of m.
CoolingLengthICKN = CoolingTimeScale(1000/x0,EnergyLossRateICKN,Usedn0)*c # The cooling length of electrons in the KN-range in units of m.

## Sample values and labels

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    gamma = 30000.0       # The used value for ga in the scenario of inverse-Compton-scattering and pair-production. INPUT VALUE!
    NumOfDecadesOfgamma = np.log10(gamma)      # This is just an auxiliary quantity, which gives the number of decades, that are included in the interval, which is going from 1 to gamma.
    
    xgamma = 1000000.0    # The used value for xga in the scenario of inverse-Compton-scattering and pair-production. INPUT VALUE!
    NumOfDecadesOfxgamma = np.log10(xgamma)      # This is just an auxiliary quantity, which gives the number of decades, that are included in the interval, which is going from 1 to xgamma.
    
    if Usedn0 == n0Planck:
        LabelForFixedgamma = r"$\gamma = 10^{%.2f}$,  $n_0(x)$: Planck-spectrum with $T = 10^{%.2f}$K and $x_0 = {%d}\cdot \mathrm{PlanckTheta}$ and $x_1 = 10^{%.2f}x_0$" % (NumOfDecadesOfgamma,np.log10(PlanckTemperature),x0NonDelta/PlanckTheta,np.log10(x1/x0NonDelta))       # For pretty labelling of figures.
        ExPsEnergyDensity = PlanckEnergyDensity
    elif Usedn0 == n0PL:
        LabelForFixedgamma = r"$\gamma = 10^{%.2f}$,  $n_0(x) = A \cdot x^{\alpha}$ with ${\alpha} = {%.2f}$ and $x_0 = 10^{%.2f}$ and $x_1 = 10^{%.2f}x_0$" % (NumOfDecadesOfgamma,alphaPL,np.log10(x0NonDelta),np.log10(x1/x0NonDelta))   # For pretty labelling of figures.
    elif Usedn0 == n0Exp:
        LabelForFixedgamma = r"$\gamma = 10^{%.2f}$,  $n_0(x) = A \cdot \exp(\alpha \cdot x)$ with $\alpha = {%.2f}$ and $x_0 = 10^{%.2f}$ and $x_1 = 10^{%.2f}x_0$" % (NumOfDecadesOfgamma,alphaExp,np.log10(x0NonDelta),np.log10(x1/x0NonDelta))   # For pretty labelling of figures.
    elif Usedn0 == n0Delta:
        LabelForFixedgamma = r"$\gamma = 10^{%.2f}$,  $n_0(x) = A \cdot \delta(x-x_0)$ with $x_0 = 10^{%.2f}$" % (NumOfDecadesOfgamma,np.log10(x0Delta))   # For pretty labelling of figures.
    
    if Usedn0 == n0Planck:
        LabelForFixedxgamma = r"$x_{\gamma} = 10^{%.2f}$,  $n_0(x)$: Planck-spectrum with $T = 10^{%.2f}$ and $x_0 = {%.2f}\cdot \mathrm{PlanckTheta}$ and $x_1 = 10^{%.2f}x_0$" % (NumOfDecadesOfxgamma,np.log10(PlanckTemperature),x0NonDelta/PlanckTheta,np.log10(x1/x0NonDelta))       # For pretty labelling of figures.
    elif Usedn0 == n0PL:
        LabelForFixedxgamma = r"$x_{\gamma} = 10^{%.2f}$,  $n_0(x) = A \cdot x^{\alpha}$ with ${\alpha} = {%.2f}$ and $x_0 = 10^{%.2f}$ and $x_1 = 10^{%.2f}x_0$" % (NumOfDecadesOfxgamma,alphaPL,np.log10(x0NonDelta),np.log10(x1/x0NonDelta))   # For pretty labelling of figures.
    elif Usedn0 == n0Exp:
        LabelForFixedxgamma = r"$x_{\gamma} = 10^{%.2f}$,  $n_0(x) = A \cdot \exp(\alpha \cdot x)$ with $\alpha = {%.2f}$ and $x_0 = 10^{%.2f}$ and $x_1 = 10^{%.2f}x_0$" % (NumOfDecadesOfxgamma,alphaExp,np.log10(x0NonDelta),np.log10(x1/x0NonDelta))   # For pretty labelling of figures.
    elif Usedn0 == n0Delta:
        LabelForFixedxgamma = r"$x_{\gamma} = 10^{%.2f}$,  $n_0(x) = A \cdot \delta(x-x_0)$ with $x_0 = 10^{%.2f}$" % (NumOfDecadesOfxgamma,np.log10(x0Delta))   # For pretty labelling of figures.

## Escape

def ElectronEscapeRate(ga):
    '''This gives the spectral probability-rate for losses of electrons due to escape from the interaction region. ga is the energy of the electrons. It is in units of 1/s.'''
    return 1.0/MeanEscapeTime

def HEPhotonEscapeRate(xga):
    '''This gives the spectral probability-rate for losses of high-energetic photons due to escape from the interaction region. xga is the energy of the photons. It is in units of 1/s.'''
    return 1.0/MeanEscapeTime

def SynchrotronPhotonEscapeRate(xSync):
    '''This gives the spectral probability-rate for losses of synchrotron-photons due to escape from the interaction region. xSync is the energy of the photons. It is in units of 1/s.'''
    return 1.0/MeanEscapeTime


## Implementation of various formulae on synchrotron-radiation

NonRelativisticGyroFrequency = e*B4/me # The cyclotron-frequency, in units of Hz.

def TimeScaleSynchrotron(B4,ga,alpha):
    '''This determines the time-scale during which synchrotron-radiation is dominant over curvature-radiation according to 2015MNRAS.447.1164V. B4 might be the magnetic flux-density in 10**4 Gauss, ga is the initial Lorentz-factor and alpha is the initial pitch-angle.'''
    return (5.0)/(B4**2.0*ga*(np.sin(alpha)**2))

def TimeScaleSynchrotronAveraged(B4,ga):
    '''This determines the averaged time-scale during which synchrotron-radiation is dominant over curvature-radiation according to 2015MNRAS.447.1164V. B4 might be the magnetic flux-density in 10**4 Gauss, ga is the initial Lorentz-factor and it was averaged over the initial pitch-angle alpha.'''
    return TimeScaleSynchrotron(B4,ga,np.pi/2.0)*3.0/2.0

SynchrotronEnergyLossRateCoefficient = (4.0*sigmaT)/(3.0*me*c) # Abbreviation for the subsequent definition.

def SynchrotronEnergyLossRateAbsoluteValue(ga):
    '''This is the absolute value of the energy-loss-rate of an ultra-relativistic electron (or positron) due to synchrotron-radiation in the magnetic field of energy-density MagneticEnergyDensity. This equation is based on equation 8.9 of the book "High Energy Astrophysics" by Malcolm Longair, up to the division by m*c^2 to yield the loss-rate of the Lorentz-factor (normalised energy). It was averaged over the pitch-angle for an isotropic distribution. The units are 1/s.'''
    return SynchrotronEnergyLossRateCoefficient*MagneticEnergyDensity*ga**2

def SynchrotronSpectralProductionRate(ga,NElectrons,ValuesForga=None):
    '''This is an additional term in the electrons' kinetic equation. It describes the increase of the spectral number-density per time at energy ga due to synchrotron-radiation energy-losses of high-energetic electrons. It is in units of 1/(m^3*s).
    It can be evaluated at any ga. In this case ValuesForga mustn't be specified. However, in case of the iteration it is evaluated only at the sampling points, which in this case have to be given via the argument ValuesForga.'''
    if isinstance(ValuesForga, np.ndarray): # In this case, ga is an element of ValuesForga.
        IndexOfThisga = np.where(ValuesForga==ga)[0][0] # the index of ga in ValuesForga.
        if IndexOfThisga==(len(ValuesForga)-1): # In this case, ga is the biggest value of ValuesForga, hence an even higher value of ga does not exist.
            return 0 # The synchrotron energy-loss from higher values of ga is = 0, because NElectrons is vanishing there.
        else: # In this case, ga is not the biggest value of ValuesForga.
            IndexOfNextHigherga=IndexOfThisga+1 # The index of the value of gamma immedeately above the variable ga.
            Denominator = ValuesForga[IndexOfNextHigherga]-ValuesForga[IndexOfThisga]
            Numerator = NElectrons(ValuesForga[IndexOfNextHigherga])*SynchrotronEnergyLossRateAbsoluteValue(ValuesForga[IndexOfNextHigherga]) # Evaluation at next higher ga.
            return Numerator/Denominator
    else:
        Denominator = 0.05*ga # This is chosen quite arbitrarily.
        Numerator = NElectrons(ga+Denominator)*SynchrotronEnergyLossRateAbsoluteValue(ga+Denominator) # Evaluation at ga+Denominator, because this corresponds to the next higher ga.
        return Numerator/Denominator

def SynchrotronSpectralLossRate(ga,ValuesForga=None):
    '''This describes the probability-rate for loss of electrons via synchrotron-radiation, i.e. the decrease per time at energy ga. It is in units of 1/s.
    It can be evaluated at any ga. In this case ValuesForga mustn't be specified. However, in case of the iteration it is evaluated only at the sampling points, which in this case have to be given via the argument ValuesForga.'''
    if isinstance(ValuesForga, np.ndarray): # In this case, ga is an element of ValuesForga.
        IndexOfThisga = np.where(ValuesForga==ga)[0][0] # the index of ga in ValuesForga.
        if IndexOfThisga==(len(ValuesForga)-1): # In this case, ga is the biggest value of ValuesForga, hence an even higher value of ga does not exist.
            IndexOfNextLowerga=IndexOfThisga-1 # The index of the value of gamma immedeately below the variable ga.
            Denominator = ValuesForga[IndexOfThisga]-ValuesForga[IndexOfNextLowerga] # Instead of (ValuesForga[IndexOfNextHigherga]-ValuesForga[IndexOfThisga]), which should be used, the next lower gamma-difference is used. 
        else: # In this case, ga is not the biggest value of ValuesForga.
            IndexOfNextHigherga=IndexOfThisga+1 # The index of the value of gamma immedeately above the variable ga.
            Denominator = ValuesForga[IndexOfNextHigherga]-ValuesForga[IndexOfThisga]
        Numerator = SynchrotronEnergyLossRateAbsoluteValue(ga) # Evaluation at ga.
        return Numerator/Denominator
    else:
        Denominator = 0.05*ga # This is chosen quite arbitrarily.
        Numerator = SynchrotronEnergyLossRateAbsoluteValue(ga) # Evaluation at ga.
        return Numerator/Denominator

SynchrotronCriticalFrequencyCoefficient = 3.0*np.pi/(2.0*4.0)*NonRelativisticGyroFrequency # Abbreviation for the subsequent definition. The pi/4 comes from averaging sin(pitch-angle) over an isotropic distribution of electrons. This is in units of Hz.

def SynchrotronCriticalFrequency(ga):
    '''This is the critical frequency of synchrotron-radiation in units of Hz, of an electron of Lorentz-factor ga. It is equation 8.55 of the book "High Energy Astrophysics" by Malcolm Longair, averaged over an isotropic distribution of pitch-angles.'''
    return SynchrotronCriticalFrequencyCoefficient*ga**2

SynchrotronCriticalxSyncCoefficient = SynchrotronCriticalFrequencyCoefficient*h/(me*c**2) # Abbreviation for the subsequent definition. This is dimensionless.

def SynchrotronCriticalxSync(ga):
    '''This is just the conversion of frequency to dimensionless energy, to save speed up computations and to make coding shorter.'''
    return SynchrotronCriticalxSyncCoefficient*ga**2

def SynchrotronEmissivityIntegrand(Argument):
    return ModifiedBesselSecondKindRealOrder(5.0/3.0,Argument)

SynchrotronEmissivityCoefficient = (2.0*np.sqrt(3.0)*e**3*B4)/(3.0*8.0*np.pi**2*epsilon0*c*me*h) # Abbreviation for the subsequent definition. The 2/3 comes from averaging over an isotropic distribution of electrons. An additional division by c**2*me comes from the conversion to dimensionless energy. An additional multiplication with (c**2*me)/h comes from the conversion of spectral in frequency to spectral in dimensionless energy. The two c**2*me cancel each other. The unit of this quantity is 1/s.

def SynchrotronEmissivity(ga,xSync,RelativeError=integrate.quad.__defaults__[3]):
    '''This gives the dimensionless energy that is emitted as synchrotron-radiation per unit time and per unit dimensionless energy by one electron (positron). It is averaged over an isotropic distribution of pitch-angles. It is based on equation 8.58 of the book "High Energy Astrophysics" by Malcolm Longair, where it is called an emissivity. (In contrast, according to Rieger an emissivity includes multiplication with the spectral number-density of electrons and integration along their energy.) The unit of this quantity is 1/s.'''
    RelativeEnergy = xSync/SynchrotronCriticalxSync(ga) # This is denoted with x by Longair.
    LowestIntegrationBorder = RelativeEnergy
    BiggestIntegrationBorder = 1000*RelativeEnergy # In truth this should be infinity. It might be a delicate choice.
    ListOfIntegrationBorders = SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
    SynchrotronEmissivityIntegral = 0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        SynchrotronEmissivityIntegral += integrate.quad(SynchrotronEmissivityIntegrand, LeftBorder, RightBorder, epsrel=RelativeError)[0]
    return SynchrotronEmissivityCoefficient*RelativeEnergy*SynchrotronEmissivityIntegral

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n\n----------------- ---------- Consider synchrotron-radiation -----------------------------\n')
    print('Test for SynchrotronEnergyLossRateAbsoluteValue:', SynchrotronEnergyLossRateAbsoluteValue(gamma),'1/s')
    print('Test for SynchrotronCriticalFrequency:          ', SynchrotronCriticalFrequency(gamma),'Hz')
    print('Test for SynchrotronCriticalxSync:              ', SynchrotronCriticalxSync(gamma))
    print('Test for SynchrotronEmissivity:                 ', SynchrotronEmissivity(gamma,SynchrotronCriticalxSync(gamma)),'1/s')

    def EvaluateSynchrotronSpectralProductionRate():
        ValuesForgammaLog = np.logspace(0.5*Lowestga,2*CurrentNElectronsga0,5000)  # An almost continuous sampling-range for ga, i.e. not only over the sampling points.
        ValuesForSynchrotronSpectralProductionRateContinuous = [SynchrotronSpectralProductionRate(i,CurrentNElectrons) for i in ValuesForgammaLog]
        ValuesForSynchrotronSpectralProductionRateDiscrete = [SynchrotronSpectralProductionRate(i,CurrentNElectrons,CurrentValuesForgaIt) for i in CurrentValuesForgaIt]
        pl.figure(figsize=(12, 9), num="Cascade-equation: Synchrotron-term versus final electron energy")
        pl.loglog(ValuesForgammaLog*x0,ValuesForSynchrotronSpectralProductionRateContinuous, '-')
        pl.loglog(CurrentValuesForgaIt*x0,ValuesForSynchrotronSpectralProductionRateDiscrete, 'r.')
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        #pl.legend(loc="best", fontsize=13)
        pl.xlabel("Final electron energy $\gamma' \cdot x_0$", fontsize=16)
        pl.ylabel("Synchrotron-term in s$^{-1} \cdot $m$^{-3}$", fontsize=16)
        pl.savefig("Cascade-equation - Synchrotron-term versus electron energy from %s.svg" % GetCurrentDate())
    
    def EvaluateSynchrotronEmissivity():
        ValuesForxSyncLog = np.logspace(np.log10(0.0001*SynchrotronCriticalxSync),np.log10(30*CriticalxSync),5000)
        ValuesForSynchrotronEmissivity = [SynchrotronEmissivity(gamma,xSync) for xSync in ValuesForxSyncLog]
        pl.figure(figsize=(12, 9), num="Cascade-equation: Synchrotron-emissivity versus photon energy")
        pl.loglog(ValuesForxSyncLog,ValuesForSynchrotronEmissivity, label="$B = %s \mathrm{T}, \; \gamma = %s$" % (B4, gamma))
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        pl.legend(loc="best", fontsize=13)
        pl.xlabel("Dimensionless photon energy $x$", fontsize=16)
        pl.ylabel("Synchrotron-emissivity in s$^{-1}$", fontsize=16)
        pl.savefig("Cascade-equation - Synchrotron-emissivity versus photon energy from %s.svg" % GetCurrentDate())
        pl.figure(figsize=(12, 9), num="Cascade-equation: Synchrotron-emissivity versus photon energy (normalised)")
        pl.loglog(ValuesForxSyncLog/CriticalxSync,ValuesForSynchrotronEmissivity/SynchrotronEmissivityCoefficient)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        pl.xlabel("Dimensionless photon energy divided by dimensionless critical photon energy $x / x_{\mathrm{crit}}$", fontsize=16)
        pl.ylabel("Normalised synchrotron-emissivity", fontsize=16)
        pl.savefig("Cascade-equation - Synchrotron-emissivity versus photon energy (normalised) from %s.svg" % GetCurrentDate())


## Considering various inverse-Compton-scattering formulae according to 1988ApJ...335..786Z

def EOfA2(ga,xOrxSync):
    return ga*xOrxSync # Equation A2, this is just a definition, it is in units 1. ga is the incoming electron's energy gamma. This definition is equally valid for xSync instead of x. So, xOrxSync is either x or xSync.

def EAstOfA3(ga,gaP):
    return (ga/gaP - 1.0)/4.0      # Equation A3, this is just a definition, it is in units 1. gaP is the scattered electron's energy.

def EAstOfA3Alt(ga,xga):
    return (ga/(ga-xga) - 1.0)/4.0      # Alternative to equation A3. It is in units 1. xga is the scattered photon's energy. 

def rOfA4(ga,gaP):
    return 0.5*(ga/gaP + gaP/ga) # Equation A4, this is just a definition, it is in units 1.

def gaPminOf6(x0OrxSync0Used,ga):        # This is equation 6 and describes the minimally possible value of gaP for a certain x0 and a certain ga. It is identical with the later defined LowerBorderOfNonVanRangeOfCOfA1 (which actually is obsolete now) and in the case n0!=n0Delta it describes the lower border of the non-vanishing range of COfA1. This definition is equally valid for xSync0Used instead of x0. So, x0OrxSync0Used is either x0 or xSync0Used. This is called gamma'_{IC, min} in the PhD thesis.
    return ga/(1.0 + 4.0*x0OrxSync0Used*ga)

def GammaLimit(x0OrxSync0Used,gaP):      # This is a limit for the value of ga in dependence on x0 and gaP. According to my analysis, the following limitations hold: Firstly ga has to be smaller than this limit for the case gaP<1/(4*x0), whilst secondly ga has to be bigger than this limit for the case gaP>1/(4*x0). Additionally ga has always to be bigger than gaP. This makes the second limitation superfluous as ga>gaP is more stringent than ga>GammaLimit. Again, this definition is equally valid for xSync0Used instead of x0. This is called gamma_{IC, max} in the PhD thesis.
    return 1.0/(1.0/gaP-4.0*x0OrxSync0Used)

def GammaMinOnxga(x0OrxSync0Used,xga):      # This is the lower border of ga in dependence on xga for a certain x0. For big enough absolute values of xga, it is asymtotically to xga+1/(4*x0). This definition is equally valid for xSync0Used instead of x0. This is called gamma_{IC, th} in the PhD thesis.
    return xga/2.0*(1.0+np.sqrt(1.0+1.0/(x0OrxSync0Used*xga)))

def xmaxOfA20(ga,x0OrxSync0Used):    # The maximum kinematically allowed value of the scattered photon-energy according to equation A20, which can either be yielded from GammaMinOnxga or from gaPminOf6. This definition is equally valid for xSync0Used instead of x0. This is called x_{gamma, max} in the PhD thesis.
    return (4.0*ga**2.0*x0OrxSync0Used)/(1.0 + 4.0*ga*x0OrxSync0Used)

# Obsolete comment: It was realised, that (for n0!=n0Delta) COfA1 is non-vanishing, only if the lower integration boundary is smaller than x0. This is equivalent with gaP > 1/(4*x0+1/ga). So, for latter convenience define the lower boundary of the non-vanishing range in gaP of the function COfA1:
def LowerBorderOfNonVanRangeOfCOfA1(x0OrxSync0Used,ga):
    return 1.0/(4.0*x0OrxSync0Used+1.0/ga)    # This is equation 6. Now, from version 9 on, this definition is obsolete. The above definition of gaPminOf6 replaces it. This definition is equally valid for xSync0Used instead of x0.

IntegrandOfA1Coefficient = 3.0*sigmaT*c/4.0 # Abbreviation for the following function and for IntegrandOfB1.

def PartOfIntegrandOfA1(xOrxSync,ga,gaP):
    '''This is the part behind n_0 of the integrand of eq. A1 of Zdz. This definition is equally valid for xSync instead of x.'''
    EAstOverE = EAstOfA3(ga,gaP)/EOfA2(ga,xOrxSync) # An abbreviation for the subsequent expression, where it could save evaluation time. Added in v2_2.
    return IntegrandOfA1Coefficient/(EOfA2(ga,xOrxSync)*ga) * ( rOfA4(ga,gaP) + (2.0-rOfA4(ga,gaP))*EAstOverE - 2.0*EAstOverE**2.0 - 2.0*EAstOverE*np.log(1.0/EAstOverE) )

def PartOfIntegrandOfA1InFrontOfBrackets(xOrxSync,ga):
    '''This is the part behind n_0 and in front of the brackets of the integrand of eq. A1 of Zdz. This definition is equally valid for xSync instead of x.'''
    return IntegrandOfA1Coefficient/(EOfA2(ga,xOrxSync)*ga)

def IntegrandOfA1(x,ga,gaP,n0):     # The part of equation A1 behind dx. According to a simple dimension analysis, this is in units 1/s. According to my analysis, it describes the double-spectral probability-rate for IC-scattering events of an electron with energy ga to final energy gaP off an external, soft photon background of spectral number-density n0. Thus it is a probability per unit time dt, per unit photon-energy dx and per unit final electron-energy dgaP.
    if gaP < gaPminOf6(x,ga): # This case has to be treated separately, to prevent that negative values are adopted.
        return 0.0
    elif gaP < ga: # The case gaP = gaPminOf6(x,ga) is included here, because the integrand is equal to 0 at ga=gamma and gaP=gaPminOf6(x,gamma).
        return n0(x) * PartOfIntegrandOfA1(x,ga,gaP) # In v2_2, the explicit definition was substituted by this function.
    elif gaP == ga:
        return n0(x) * PartOfIntegrandOfA1InFrontOfBrackets(x,ga) # In v2_2, the explicit definition was substituted by this function.
    else: # The integrand cannot be evaluated at gaP>ga.
        return 0.0

def IntegrandOfA1WithSynchrotron(xSync,ga,gaP,nSyncPs):     # The analogon to the part of equation A1 behind dx, just with the synchrotron-photon-field. This is in units 1/s. It describes the double-spectral probability-rate for IC-scattering events of an electron with energy ga to final energy gaP off a soft synchrotron-photon background of spectral number-density nSyncPs. Thus it is a probability per unit time dt, per unit SyncP-energy dxSync and per unit final electron-energy dgaP. nSyncPs(xSync) is an interpolated object with kinks.
    if gaP < gaPminOf6(xSync,ga): # This case has to be treated separately, to prevent that negative values are adopted.
        return 0.0
    elif gaP < ga: # The case gaP = gaPminOf6(xSync,ga) is included here, because the integrand is equal to 0 at ga=gamma and gaP=gaPminOf6(xSync,gamma).
        return nSyncPs(xSync) * PartOfIntegrandOfA1(xSync,ga,gaP) # In v2_2, the explicit definition was substituted by this function.
    elif gaP == ga:
        return nSyncPs(xSync) * PartOfIntegrandOfA1InFrontOfBrackets(xSync,ga) # In v2_2, the explicit definition was substituted by this function.
    else: # The integrand cannot be evaluated at gaP>ga.
        return 0.0

# Equation A1. According to a simple dimension analysis, this is in units 1/s. According to my analysis, it describes the spectral probability-rate for IC-scattering events of an electron with energy ga to final energy gaP off a soft photon background of spectral number-density n0. Thus it is a probability per unit time dt and per unit final electron-energy dgaP. The kinematically motivated lower integration border is EAstOfA3(ga,gaP)/ga.
if Usedn0 == n0Delta:
    def COfA1WithoutSynchrotron(ga,gaP,n0,RelativeError=integrate.quad.__defaults__[3]):   # Defining COfA1WithoutSynchrotron separately for the case of n0 being a Dirac-Delta-function is a trick to bypass the integration via making use of the substitution-feature of the Dirac-Delta-function.
        return IntegrandOfA1(x0,ga,gaP,n0) # The case x0<EAstOfA3(ga,gaP)/ga needn't be excluded here via an if-test because it is already treated and set =0 in IntegrandOfA1.
elif Usedn0 == n0MultiDelta:
    def COfA1WithoutSynchrotron(ga,gaP,n0,RelativeError=integrate.quad.__defaults__[3]):   # Defining COfA1WithoutSynchrotron separately for the case of n0 being a set of Dirac-Delta-functions again makes use of bypassing the integration via making use of the substitution-feature. This is done for each Delta-peak and then the sum of them is determined.
        COfA1WithoutSynchrotronSum = 0.0
        for i in x0MultiDelta:
            COfA1WithoutSynchrotronSum += IntegrandOfA1(i,ga,gaP,n0) # The case x0<EAstOfA3(ga,gaP)/ga needn't be excluded here via an if-test because it is already treated and set =0 in IntegrandOfA1.
        return COfA1WithoutSynchrotronSum
elif Usedn0==n0Exp or Usedn0==n0PL or Usedn0==n0Planck:  # In the following cases, the integration has to be performed numerically.
    def COfA1WithoutSynchrotron(ga,gaP,n0,RelativeError=integrate.quad.__defaults__[3]):
        LowestIntegrationBorder = max(x1,EAstOfA3(ga,gaP)/ga)
        BiggestIntegrationBorder = x0
        if LowestIntegrationBorder >= BiggestIntegrationBorder:
            return 0.0
        else:
            NumberOfBorders = max(6,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
            ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
            COfA1WithoutSynchrotronResult = 0.0
            for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                COfA1WithoutSynchrotronResult += integrate.quad(IntegrandOfA1, LeftBorder, RightBorder, args=(ga,gaP,n0), epsrel=RelativeError)[0]
            return COfA1WithoutSynchrotronResult
elif Usedn0==n0AF: # This case is still experimental.
    def COfA1WithoutSynchrotron(ga,gaP,n0,RelativeError=integrate.quad.__defaults__[3]):
        LowestIntegrationBorder = max(ADAFx1,EAstOfA3(ga,gaP)/ga)
        BiggestIntegrationBorder = x0
        if LowestIntegrationBorder >= BiggestIntegrationBorder:
            return 0.0
        else:
            NumberOfBorders = max(6,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
            ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
            for PotentialBorder in [ADAFxMin,ADAFxPeak,ADAFxComptonisedCutOff]:
                if LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder and (PotentialBorder not in ListOfIntegrationBorders):
                    ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
            ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
            COfA1WithoutSynchrotronResult = 0.0
            for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                COfA1WithoutSynchrotronResult += integrate.quad(IntegrandOfA1, LeftBorder, RightBorder, args=(ga,gaP,n0), epsrel=RelativeError)[0]
            return COfA1WithoutSynchrotronResult
    # In contrast to version 1, the cut-off of the spectral number-densities is not included via Heaviside-functions but via changed integration borders. Effectively, it is integrated only over the intersection of the original integration range with the non-vanishing part of the spectral number-density. So the new lower integration border is the maximum of the old lower integration border and the lower cut-off of the spectral number-density while the new upper integration border is the upper cut-off of the spectral number-density. The old lower integration border is dependent on ga and on gaP. Hence, it can be bigger than the new upper integration border. In this case, it must be prevented to exceed the new upper integration border. 
    # Comment of version 16: The if-test has been changed such that it immediately returns 0 in the case LowerIntBorder >= x0.
    # Comment of version 19: The integration range is divided into part-ranges.
    # Comment of v2_3: Was renamed from COfA1 to COfA1WithoutSynchrotron.

def COfA1WithSynchrotron(ga,gaP,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''The analogue to equation A1 just with the synchrotron-radiation as soft photons. This is in units 1/s. It describes the spectral probability-rate for IC-scattering events of an electron with energy ga to final energy gaP off a soft SyncP-background of spectral number-density nSyncPs. Thus it is a probability per unit time dt and per unit final electron-energy dgaP. The kinematically motivated lower integration border is EAstOfA3(ga,gaP)/ga. However, nSyncPs is non-vanishing only in between of ValuesForxSyncSuperIt, hence one needs to integrate only along this range.'''
    LowestIntegrationBorder = max(ValuesForxSyncSuperIt[0],EAstOfA3(ga,gaP)/ga) # Actually, the first argument is the current xSync1.
    BiggestIntegrationBorder = ValuesForxSyncSuperIt[-1] # Actually, ValuesForxSyncSuperIt[-1] is the current xSync0Used and is the (physically) correct upper border.
    if LowestIntegrationBorder >= BiggestIntegrationBorder:
        return 0.0
    else:
        if LowestIntegrationBorder<0.01*BiggestIntegrationBorder: # Comment of v2_5: By evaluating EvaluatePartOfIntegrandOfA1 for reasonable parameters of ga and gaP, one can recognise that PartOfIntegrandOfA1 always has its peak slightly right of the lower border (EAstOfA3(ga,gaP)/ga) of its non-vanishing range. Above the peak, it is a power-law with slope of approx -1. If PartOfIntegrandOfA1 is multiplied with nSyncPs, the slope of the IntegrandOfA1 is even reduced. Therefore, for those cases where the peak of PartOfIntegrandOfA1 is low enough, it is not necessary to extend the integration range up to ValuesForxSyncSuperIt[-1]. Hence, it is chosen that if the lower border is lower than one tenth of the upper border, then the upper border is reduced by factor 0.5 and if the lower border is lower than one hundredth of the upper border, then the upper border is reduced by factor 0.1.
            NumericalBiggestIntegrationBorder=0.1*BiggestIntegrationBorder
            ListOfIntegrationBorders = np.asarray([LowestIntegrationBorder,NumericalBiggestIntegrationBorder])
        elif LowestIntegrationBorder<0.1*BiggestIntegrationBorder:
            NumericalBiggestIntegrationBorder=0.5*BiggestIntegrationBorder
            ListOfIntegrationBorders = np.asarray([LowestIntegrationBorder,NumericalBiggestIntegrationBorder])
        else:
            NumericalBiggestIntegrationBorder=BiggestIntegrationBorder
            ListOfIntegrationBorders = np.asarray([LowestIntegrationBorder,NumericalBiggestIntegrationBorder])
        for PotentialBorder in ValuesForxSyncSuperIt[1:-1]: # At the elements of ValuesForxSyncSuperIt, nSyncPs has kinks.
            if LowestIntegrationBorder<PotentialBorder<NumericalBiggestIntegrationBorder:
                ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
        ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
        COfA1WithSynchrotronResult = 0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            COfA1WithSynchrotronResult += integrate.quad(IntegrandOfA1WithSynchrotron, LeftBorder, RightBorder, args=(ga,gaP,nSyncPs), epsrel=RelativeError)[0]
        return COfA1WithSynchrotronResult

if IncludeSynchrotron and PerformanceMode!='CommonIteration':
    def COfA1(ga,gaP,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
        return COfA1WithoutSynchrotron(ga,gaP,n0,RelativeError)+COfA1WithSynchrotron(ga,gaP,nSyncPs,ValuesForxSyncSuperIt,RelativeError)
else:
    def COfA1(ga,gaP,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]): # For the usage of COfA1 (and all functions that are built on it) in the case without synchrotron-back-reaction, one should just use nSyncPs=None and ValuesForxSyncSuperIt=None.
        return COfA1WithoutSynchrotron(ga,gaP,n0,RelativeError)

def COfA1DependingOnxgaWithoutSynchrotron(ga,xga,n0,RelativeError=integrate.quad.__defaults__[3]):
    '''This is a contribution to COfA1DependingOnxga. It is in units 1/s. It describes the spectral probability-rate for IC-scattering events of an electron with energy ga off the external soft photon background of spectral number-density n0 to a scattered photon of energy xga. Thus it is a probability per unit time dt and per unit final photon energy xga.'''
    return COfA1WithoutSynchrotron(ga,ga-xga,n0,RelativeError)

def COfA1DependingOnxgaWithSynchrotron(ga,xga,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is a contribution to COfA1DependingOnxga. It is in units 1/s. It describes the spectral probability-rate for IC-scattering events of an electron with energy ga off a soft SyncP-background of spectral number-density nSyncPs to a scattered photon of energy xga. Thus it is a probability per unit time dt and per unit final photon energy xga.'''
    return COfA1WithSynchrotron(ga,ga-xga,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def COfA1DependingOnxga(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    return COfA1(ga,ga-xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)     # This is a variation of equation A1. It is in units 1/s. It might describe the spectral probability-rate for IC-scattering events of an electron with energy ga off a soft photon background of spectral number-density n0 (plus a soft SyncP-background of spectral number-density nSyncPs) to a scattered photon of energy xga. Thus it is a probability per unit time dt and per unit final photon energy xga.

def IntegrandOfA21WithoutSynchrotron(gaP,ga,n0,RelativeError=integrate.quad.__defaults__[3]):      # ga and gaP are reversed, because in the following integral it is to be integrated over gaP, so this has to be the first argument. 
    return COfA1WithoutSynchrotron(ga,gaP,n0,RelativeError)
    
def IntegrandOfA21WithSynchrotron(gaP,ga,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):      # ga and gaP are reversed, because in the following integral it is to be integrated over gaP, so this has to be the first argument. 
    return COfA1WithSynchrotron(ga,gaP,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def IntegrandOfA21(gaP,ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):      # ga and gaP are reversed, because in the following integral it is to be integrated over gaP, so this has to be the first argument. 
    return COfA1(ga,gaP,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def COfA21WithoutSynchrotron(ga,n0,RelativeError=integrate.quad.__defaults__[3]):     # Equation A21. According to a simple dimension analysis, this is in units 1/s, too. According to my analysis, it describes the probability-rate for IC-scattering events of an electron with energy ga to any final energy off solely a soft photon background of spectral number-density n0. Thus it is a probability per unit time dt.
# According to Zdz. one should integrate from 1 to ga. If one does so, the following conclusions can be drawn: For n0=n0Delta it works well but output becomes negative, thus seems unphysical. Thus, in contrast to equation A21, it is integrated from gaPminOf6 to ga, as this seems physically more reasonable. As gaPminOf6 is the site of the zero of COfA1WithoutSynchrotron the result has to be positive. For n0=n0Exp, the integration seems to have problems, if one integrates over the site, at which COfA1WithoutSynchrotron jumps from 0 to a non-vanishing value. For n0=n0Planck it is similar: The integration has problems if the integration range extends beyond the jumping point. For n0=n0PL it is quite obscure: The integration seems to have problems for ga=10**8, if the integration range includes both sites gaP approx 2000 and gaP approx 1.2*10**8.
    LowestIntegrationBorder = gaPminOf6(x0,ga)
    BiggestIntegrationBorder = ga
    ListOfIntegrationBorders = [LowestIntegrationBorder,BiggestIntegrationBorder]
    PotentialAdditionalBorders = [] # The integrand can have kinks, due to additional contributions to COfA1. These additional contributions are caused by discontinuities along x.
    if Usedn0==n0MultiDelta:
        PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[gaPminOf6(i,ga) for i in x0MultiDelta])
    else:
        PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[gaPminOf6(x0,ga)])
    for PotentialBorder in PotentialAdditionalBorders: # Include these points into the list of integration-borders, if they are in the interior of the range.
        if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
            ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)        
    ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
    # Now, subdivide each of the arising ranges and take the divisions as additional borders:
    AdditionalBorders = []
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        AdditionalBorders=np.append(AdditionalBorders,SampleIntegrationBorders10(LeftBorder,RightBorder)[0][1:-1])
    for PotentialBorder in AdditionalBorders: # Include these borders into the list of integration-borders, if they are not included yet.
        if PotentialBorder not in ListOfIntegrationBorders:
            ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
    ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
    COfA21WithoutSynchrotronResult = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        COfA21WithoutSynchrotronResult += integrate.quad(IntegrandOfA21WithoutSynchrotron, LeftBorder, RightBorder, args=(ga,n0,RelativeError), epsrel=RelativeError)[0]   # ...Then, the integration over the entire range is determined via the sum of integrations over the small ranges. 
    return COfA21WithoutSynchrotronResult
    # Comment of version 15: Sometimes even for n0Delta, the integration-range has to be subdivided. This was done in version 16 by removing the if-test.
    # Comment of version 19: During pion decay modelling, an irregularity became obvious, which was due to an error in COfA21WithoutSynchrotron. Therefore the minimum number of integration borders was raised to 7.

def COfA21WithSynchrotron(ga,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''The analogue to equation A21 just for the case of a pure synchrotron background. This is in units 1/s, too, and describes the probability-rate for IC-scattering events of an electron with energy ga to any final energy off a soft SyncP-background nSyncPs. Thus it is a probability per unit time dt.'''
    LowestIntegrationBorder = gaPminOf6(ValuesForxSyncSuperIt[-1],ga) # Actually, the first argument is the current xSync0Used.
    BiggestIntegrationBorder = ga
    ListOfIntegrationBorders = SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
    COfA21WithSynchrotronResult = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        COfA21WithSynchrotronResult += integrate.quad(IntegrandOfA21WithSynchrotron, LeftBorder, RightBorder, args=(ga,nSyncPs,ValuesForxSyncSuperIt,RelativeError), epsrel=RelativeError)[0]
    return COfA21WithSynchrotronResult

def COfA21(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''Equation A21 with the used soft photons (hence soft background photons or soft photons plus SyncPs). This is in units 1/s, too. It describes the probability-rate for IC-scattering events of an electron with energy ga to any final energy off a soft photon background of spectral number-density n0 (plus a soft SyncP-background nSyncPs). Thus it is a probability per unit time dt. For the usage of COfA21 (and all functions that are built on it) in the case without synchrotron, one should just use nSyncPs=None and ValuesForxSyncSuperIt=None.'''
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        LowestIntegrationBorder = min(gaPminOf6(x0,ga),gaPminOf6(ValuesForxSyncSuperIt[-1],ga)) # This is valid only if the two summands of the integrand vanish outside of their definition range.
    else:
        LowestIntegrationBorder = gaPminOf6(x0,ga)
    BiggestIntegrationBorder = ga
    ListOfIntegrationBorders = [LowestIntegrationBorder,BiggestIntegrationBorder]    
    PotentialAdditionalBorders = [] # The integrand can have kinks, due to additional contributions to COfA1. These additional contributions are caused by discontinuities along x.
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[gaPminOf6(ValuesForxSyncSuperIt[-1],ga)])
    if Usedn0==n0MultiDelta:
        PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[gaPminOf6(i,ga) for i in x0MultiDelta])
    else:
        PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[gaPminOf6(x0,ga)])
    for PotentialBorder in PotentialAdditionalBorders: # Include these points into the list of integration-borders, if they are in the interior of the range.
        if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
            ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)        
    ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)    
    # Now, subdivide each of the arising ranges and take the divisions as additional borders:
    AdditionalBorders = []
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        AdditionalBorders=np.append(AdditionalBorders,SampleIntegrationBorders10(LeftBorder,RightBorder)[0][1:-1])
    for PotentialBorder in AdditionalBorders: # Include these borders into the list of integration-borders, if they are not included yet.
        if PotentialBorder not in ListOfIntegrationBorders:
            ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
    ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
    COfA21Result = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        COfA21Result += integrate.quad(IntegrandOfA21, LeftBorder, RightBorder, args=(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError), epsrel=RelativeError)[0]
    return COfA21Result

# This is an alternative implementation of equation A21 with the used soft photons (hence soft background photons or soft photons plus SyncPs). This is physically the same as COfA21, however, it does not integrate IntegrandOfA21 (which essentially is COfA1), but it is the sum of both part-integrals COfA21WithoutSynchrotron and COfA21WithSynchrotron.
if IncludeSynchrotron and PerformanceMode!='CommonIteration':
    def COfA21Alternative(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
        return COfA21WithoutSynchrotron(ga,n0,RelativeError)+COfA21WithSynchrotron(ga,nSyncPs,ValuesForxSyncSuperIt,RelativeError)
else:
    def COfA21Alternative(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
        return COfA21WithoutSynchrotron(ga,n0,RelativeError)

def IntegrandOfDotgammaWithoutSynchrotron(gaP,ga,n0,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the absolute value of the change in gamma per unit interval gaP and per unit time interval. This is in units of 1/s. In other words, this is the absolute value of the integrand of Zdz.'s eq. 3. The following integral it is to be integrated over gaP, so this has to be the first argument.'''
    return (ga-gaP)*COfA1WithoutSynchrotron(ga,gaP,n0,RelativeError)
    
def IntegrandOfDotgammaWithSynchrotron(gaP,ga,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the absolute value of the change in gamma per unit interval gaP and per unit time interval. This is in units of 1/s. The following integral it is to be integrated over gaP, so this has to be the first argument.'''
    return (ga-gaP)*COfA1WithSynchrotron(ga,gaP,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def IntegrandOfDotgamma(gaP,ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the absolute value of the change in gamma per unit interval gaP and per unit time interval. This is in units of 1/s. The following integral it is to be integrated over gaP, so this has to be the first argument.'''
    return (ga-gaP)*COfA1(ga,gaP,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def DotgammaWithoutSynchrotron(ga,n0,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the absolute value of the change in gamma per unit time interval for IC-scattering events of an electron with energy ga to any final energy off solely a soft external photon background of spectral number-density n0. This is in units of 1/s. This is similar to Zdz.'s eq. 3, just with the lower integration border changed to the minimally possible gaP, which means that all (not only the continuous) energy-losses are integrated.'''
    LowestIntegrationBorder = gaPminOf6(x0,ga)
    BiggestIntegrationBorder = ga
    ListOfIntegrationBorders = [LowestIntegrationBorder,BiggestIntegrationBorder]
    PotentialAdditionalBorders = [] # The integrand can have kinks, due to additional contributions to COfA1. These additional contributions are caused by discontinuities along x.
    if Usedn0==n0MultiDelta:
        PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[gaPminOf6(i,ga) for i in x0MultiDelta])
    else:
        PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[gaPminOf6(x0,ga)])
    for PotentialBorder in PotentialAdditionalBorders: # Include these points into the list of integration-borders, if they are in the interior of the range.
        if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
            ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)        
    ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
    # Now, subdivide each of the arising ranges and take the divisions as additional borders:
    AdditionalBorders = []
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        AdditionalBorders=np.append(AdditionalBorders,SampleIntegrationBorders10(LeftBorder,RightBorder)[0][1:-1])
    for PotentialBorder in AdditionalBorders: # Include these borders into the list of integration-borders, if they are not included yet.
        if PotentialBorder not in ListOfIntegrationBorders:
            ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
    ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
    DotgammaWithoutSynchrotronResult = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        DotgammaWithoutSynchrotronResult += integrate.quad(IntegrandOfDotgammaWithoutSynchrotron, LeftBorder, RightBorder, args=(ga,n0,RelativeError), epsrel=RelativeError)[0]
    return DotgammaWithoutSynchrotronResult

def DotgammaWithSynchrotron(ga,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the absolute value of the change in gamma per unit time interval for IC-scattering events of an electron with energy ga to any final energy off solely a synchrotron-photon background of spectral number-density nSync. This is in units of 1/s.'''
    LowestIntegrationBorder = gaPminOf6(ValuesForxSyncSuperIt[-1],ga)
    BiggestIntegrationBorder = ga
    ListOfIntegrationBorders = SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
    DotgammaWithSynchrotronResult = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        DotgammaWithSynchrotronResult += integrate.quad(IntegrandOfDotgammaWithSynchrotron, LeftBorder, RightBorder, args=(ga,nSyncPs,ValuesForxSyncSuperIt,RelativeError), epsrel=RelativeError)[0]
    return DotgammaWithSynchrotronResult

def Dotgamma(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the absolute value of the change in gamma per unit time interval for IC-scattering events of an electron with energy ga to any final energy off a soft external photon background and a SyncP-background. This is in units of 1/s.'''
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        LowestIntegrationBorder = min(gaPminOf6(x0,ga),gaPminOf6(ValuesForxSyncSuperIt[-1],ga)) # This is valid only if the two summands of the integrand vanish outside of their definition range.
    else:
        LowestIntegrationBorder = gaPminOf6(x0,ga)
    BiggestIntegrationBorder = ga
    ListOfIntegrationBorders = [LowestIntegrationBorder,BiggestIntegrationBorder]    
    PotentialAdditionalBorders = [] # The integrand can have kinks, due to additional contributions to COfA1. These additional contributions are caused by discontinuities along x.
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[gaPminOf6(ValuesForxSyncSuperIt[-1],ga)])
    if Usedn0==n0MultiDelta:
        PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[gaPminOf6(i,ga) for i in x0MultiDelta])
    else:
        PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,[gaPminOf6(x0,ga)])
    for PotentialBorder in PotentialAdditionalBorders: # Include these points into the list of integration-borders, if they are in the interior of the range.
        if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
            ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)        
    ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)    
    # Now, subdivide each of the arising ranges and take the divisions as additional borders:
    AdditionalBorders = []
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        AdditionalBorders=np.append(AdditionalBorders,SampleIntegrationBorders10(LeftBorder,RightBorder)[0][1:-1])
    for PotentialBorder in AdditionalBorders: # Include these borders into the list of integration-borders, if they are not included yet.
        if PotentialBorder not in ListOfIntegrationBorders:
            ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
    ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
    DotgammaResult = 0.0
    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
        DotgammaResult += integrate.quad(IntegrandOfDotgamma, LeftBorder, RightBorder, args=(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError), epsrel=RelativeError)[0]
    return DotgammaResult

# This is an alternative implementation of Dotgamma with the used soft photons (hence soft background photons or soft photons plus SyncPs). This is physically the same as Dotgamma, however, it does not integrate IntegrandOfDotgamma, but it is the sum of both part-integrals DotgammaWithoutSynchrotron and DotgammaWithSynchrotron.
if IncludeSynchrotron and PerformanceMode!='CommonIteration':
    def DotgammaAlternative(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
        return DotgammaWithoutSynchrotron(ga,n0,RelativeError)+DotgammaWithSynchrotron(ga,nSyncPs,ValuesForxSyncSuperIt,RelativeError)
else:
    def DotgammaAlternative(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
        return DotgammaWithoutSynchrotron(ga,n0,RelativeError)

def ElectronLossRate(ga,n0,nSyncPs,ValuesForxSyncSuperIt,ValuesForga=None,RelativeError=integrate.quad.__defaults__[3]):
    '''This is in units of 1/s and gives the spectral probability-rate of events that effectively make electrons disappear from the energy ga. It includes IC-down-scattering (included in every case for the case of soft background photons and, depending on the value of IncludeSynchrotron, also for the case of a SyncP-background) as well as escape from the interaction region as well as energy-losses due to emission of synchrotron-radiation.'''
    if IncludeElectronEscape==False:
        if IncludeSynchrotron==False:
            return COfA21(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)
        elif IncludeSynchrotron==True:
            return COfA21(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)+SynchrotronSpectralLossRate(ga,ValuesForga)
    elif IncludeElectronEscape==True:
        if IncludeSynchrotron==False:
            return COfA21(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)+ElectronEscapeRate(ga)
        elif IncludeSynchrotron==True:
            return COfA21(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)+ElectronEscapeRate(ga)+SynchrotronSpectralLossRate(ga,ValuesForga)

def cOfA19(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):                   # Equation A19. According to my analysis, this describes a normalised spectral probability for IC-scattering events of an electron with energy ga off a soft photon background of spectral number-density n0 (plus nSyncPs) to a scattered photon of energy xga. Thus it is a probability per unit final photon energy dxga.
    return COfA1DependingOnxga(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)/COfA21(ga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def NormalisedSpectralICScatteringProbability(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,ValuesForga=None,RelativeError=integrate.quad.__defaults__[3]):        # This describes the same as cOfA19, however in contrast to cOfA19, it is normalised to the total electron-loss-rate, not only to the IC-scattering-rate.
    return COfA1DependingOnxga(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)/ElectronLossRate(ga,n0,nSyncPs,ValuesForxSyncSuperIt,ValuesForga,RelativeError)

def OpticalDepthOfElectrons(Process,*PositionalArgumentsOfProcess):
    '''This is the optical depth of the relativistic electrons with respect to a loss-process Process. It is dimensionless.
    It was yielded based on the relation loss-probability-rate-of-process = c*dtau/dlength, where tau is the optical depth and length the path-length. Solving for dtau and integrating yields tau = loss-probability-rate-of-process*MeanEscapeLength/c, where MeanEscapeLength is the mean distance an electron has to pass until it leaves the interaction region. Using c = MeanEscapeLength/MeanEscapeTime gives the following expression.
    Process can be any of the contributions to ElectronLossRate or ElectronLossRate itself, hence IC-scattering on the soft background-photons (Process=COfA21WithoutSynchrotron, arguments have to be ga and n0), IC-scattering on the synchrotron-photons (Process=COfA21WithSynchrotron, arguments have to be ga, nSyncPs and ValuesForxSyncSuperIt), IC-scattering on the soft background-photons (plus the SyncPs) (Process=COfA21, arguments have to be ga, n0, nSyncPs and ValuesForxSyncSuperIt), escape from the region (Process=ElectronEscapeRate, arguments have to be ga), energy-loss via synchrotron-radiation (Process=SynchrotronSpectralLossRate, arguments have to be ga) or all these (Process=ElectronLossRate, arguments have to be ga, n0, nSyncPs and ValuesForxSyncSuperIt).'''
    return Process(*PositionalArgumentsOfProcess)*MeanEscapeTime # This essentially is the equation tau = loss-probability-rate-of-process * MeanEscapeTime.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    
    def EvaluategaPminOf6():
        #log-log plot 
        x0=0.00001
        print('1/(4*x0) =', 1.0/(4.0*x0))
        GammaRange=np.logspace(0,9,1500)
        pl.figure(num="IC-scattering: Lower border of final electron energy gamma'")
        pl.loglog(GammaRange,gaPminOf6(x0,GammaRange),label="$\gamma'_{min} = \gamma/(1 + 4 \cdot x_0 \cdot \gamma)$")
        pl.loglog(GammaRange,np.ones(len(GammaRange))*1.0/(4.0*x0),label="$1/(4 \cdot x_0)$")
        pl.loglog(GammaRange,GammaRange,label="$\gamma$")
        pl.legend(loc="lower right", fontsize=14)
        pl.xlabel("Electron energy $\gamma$", fontsize=16)
        pl.ylabel("$\gamma'$", fontsize=16)
        pl.ylim(1.0, 1.0/(4.0*x0)*2.0)
        pl.text(2*10**5,100,r'For $x_0 = 10^{%.2f}$' % np.log10(x0), fontsize=16)
        
        # lin-lin plot
        ExtraSpace = matplotlib.patches.Rectangle((0, 0), 0.1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
        pl.rc('font', family='serif')
        for i in [0,1]:
            x0=[0.01,0.0025][i]
            Color=['red','blue'][i]
            ColouredSpace = matplotlib.patches.Rectangle((0, 0), 1, 1, fc="%s" % Color, fill=True, edgecolor='none', linewidth=0, alpha=.2)
            if i==0:
                AdditionalLinewidth=0.6
                AdditionalLinewidth2=0.9
            elif i==1:
                AdditionalLinewidth=0.0
                AdditionalLinewidth2=0.0
            print('1/(4*x0) =', 1.0/(4.0*x0))
            GammaLowerRange=np.linspace(-500,-1/(4*x0),1500,endpoint=False)
            GammaUpperRange=np.linspace(-0.999/(4*x0),1000,1500)
            GammaTotalRange=np.append(GammaLowerRange,GammaUpperRange)
            AllowedGammaRange=GammaTotalRange[GammaTotalRange>0]
            ValuesForgaPminOf6=gaPminOf6(x0,GammaUpperRange)
            AllowedValuesForgaPminOf6=ValuesForgaPminOf6[ValuesForgaPminOf6>0]
            pl.figure(figsize=(11.8, 11), num="IC-scattering: gamma'_{IC min}")
            Aa=pl.plot(GammaLowerRange,gaPminOf6(x0,GammaLowerRange),label="$\gamma ' = \gamma'_{\mathrm{IC,\,min}}(\gamma,\,x)$", color='%s' % Color, linewidth=1.5+AdditionalLinewidth)[0]
            pl.plot(GammaUpperRange,ValuesForgaPminOf6, color='%s' % Color, linewidth=1.5+AdditionalLinewidth)
            Bb=pl.plot(GammaTotalRange,np.ones(len(GammaTotalRange))*1.0/(4.0*x0),label="$\gamma ' = 1/(4 \, x)$", color='%s' % Color, linestyle=':', linewidth=2+AdditionalLinewidth)[0]
            #Cc=pl.plot(GammaTotalRange,GammaTotalRange,label="$\gamma ' = \gamma$", color='%s' % Color, linestyle='--', linewidth=1.7+AdditionalLinewidth+AdditionalLinewidth2)[0]
            #Dd=pl.plot([-1/(4*x0),-1/(4*x0)],[-10000,10000],label="$\gamma  = -1/(4 \, x)$", linestyle='--', color='%s' % Color, linewidth=1.7+AdditionalLinewidth)[0]
            pl.fill_between(AllowedGammaRange, AllowedValuesForgaPminOf6, AllowedGammaRange, alpha=.2, color='%s' % Color)
            if i==0:
                Ee=pl.plot(GammaTotalRange,GammaTotalRange,label="$\gamma' = \gamma$", linewidth=1.0, color='gray')[0]
                L1=[Ee, ExtraSpace, ExtraSpace, Aa, Bb, ColouredSpace]
                L2=[Ee.get_label(), '   ', '$x=%.4f:$' % x0, Aa.get_label(), Bb.get_label(), '$\mathrm{Allowed\,Range}$']
            elif i==1:
                L1.extend([ExtraSpace,ExtraSpace, Aa, Bb, ColouredSpace])
                L2.extend(['   ', '$x=%.4f:$' % x0, Aa.get_label(), Bb.get_label(), '$\mathrm{Allowed\,Range}$'])
                pl.legend(L1,L2, loc='best', fontsize=18)
        pl.xlabel("Original electron energy $\gamma$", fontsize=26)
        pl.ylabel("Final electron energy $\gamma'$", fontsize=26)
        pl.xlim(-220,380)
        pl.ylim(-170,440)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=26)
        ax.yaxis.set_tick_params(labelsize=26)
        pl.subplots_adjust(top=0.98, bottom=0.09, left=0.13, right=0.98)
        pl.savefig("IC scattering - gamma'_ICmin from %s.pdf" % GetCurrentDate())
    
    def EvaluateGammaLimit():
        x0=0.00001
        print('1/(4*x0) =', 1.0/(4.0*x0))
        GammaPRange=np.linspace(1.0,40.0*1.0/(4.0*x0),300000)
        ValuesForGammaLimit=GammaLimit(x0,GammaPRange)
        pl.figure(figsize=(12.5, 6), num="IC-scattering: Limit of original electron energy gamma")
        pl.subplot(1,3,1)
        pl.plot(GammaPRange,ValuesForGammaLimit,label="$\gamma_{lim} = \gamma'/(1 - 4 \cdot x_0 \cdot \gamma')$")
        pl.plot(GammaPRange,GammaPRange,label="$\gamma'$")
        pl.ylabel("$\gamma$", fontsize=16)
        pl.xlim(GammaPRange[0],GammaPRange[0.001*len(GammaPRange)])
        pl.ylim(ValuesForGammaLimit[0],ValuesForGammaLimit[0.001*len(GammaPRange)])
        pl.text(100,900,r'For $x_0 = 10^{%.2f}$' % np.log10(x0), fontsize=16)
        pl.subplot(1,3,2)
        pl.plot(GammaPRange,ValuesForGammaLimit,label="$\gamma_{lim} = \gamma'/(1 - 4 \cdot x_0 \cdot \gamma')$")
        pl.plot(GammaPRange,GammaPRange,label="$\gamma'$")
        pl.legend(loc="upper right", fontsize=14)
        pl.xlabel("Final electron energy $\gamma'$", fontsize=16)
        pl.xlim(1.0/(4.0*x0)-0.08*1.0/(4.0*x0),1.0/(4.0*x0)+0.08*1.0/(4.0*x0))
        pl.ylim(-0.01*np.max(ValuesForGammaLimit),0.01*np.max(ValuesForGammaLimit))
        pl.xticks(np.linspace(1.0/(4.0*x0)-0.08*1.0/(4.0*x0),1.0/(4.0*x0)+0.08*1.0/(4.0*x0), 3))
        pl.subplot(1,3,3)
        pl.plot(GammaPRange,ValuesForGammaLimit)
        pl.plot(GammaPRange,-np.ones(len(GammaPRange))*1.0/(4.0*x0),'r',label="$- 1/(4 \cdot x_0)$")
        pl.legend(loc="upper right", fontsize=14)
        pl.xlim(GammaPRange[0.1*len(GammaPRange)],GammaPRange[-1])
        pl.ylim(ValuesForGammaLimit[0.1*len(GammaPRange)],0)
        pl.xticks(np.linspace(GammaPRange[0.1*len(GammaPRange)]-1,GammaPRange[-1], 3))
        pl.subplots_adjust(left=0.08,wspace=0.4)
    
    def EvaluateGammaMinOnxga():
        print('1/x0 =', 1.0/x0)
        xGammaRange=np.append(np.linspace(-2.0/x0,-1.0/x0,301),np.linspace(0.001,1.0/x0,301))
        ValuesForGammaMinOnxga=GammaMinOnxga(x0,xGammaRange)
        pl.figure(figsize=(12, 9), num="IC-scattering: Minimum allowed value of gamma versus final photon energy")
        pl.subplot(1,1,1)
        pl.plot(xGammaRange,ValuesForGammaMinOnxga,linestyle='None',marker='.',markersize=1, label='$\gamma_{min}(x_\gamma)$')
        pl.plot(xGammaRange,xGammaRange+1.0/(4.0*x0),label='$x_\gamma + 1/(4 \cdot x_0)$')
        pl.plot(xGammaRange,0.5*np.sqrt(xGammaRange/x0),label='$0.5 \cdot \sqrt{x_\gamma / x_0}$ (limit for $x_\gamma x_0 << 1$)')
        pl.text(-2.0/x0,1.0/x0,r'For $1/(4 \cdot x_0) = %d $' % (1.0/(4.0*x0)), fontsize=16)
        pl.legend(loc="right", fontsize=14)
        pl.xlabel('Final photon energy $x_\gamma$', fontsize=16)
        pl.ylabel('Original electron energy $\gamma$', fontsize=16)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        pl.ylim(xGammaRange[0]+1.0/(4.0*x0),ValuesForGammaMinOnxga[-1])
        pl.subplots_adjust(left=0.2)
        
        ExtraSpace = matplotlib.patches.Rectangle((0, 0), 0.1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
        pl.rc('font', family='serif')
        pl.figure(figsize=(9.5, 8.7), num="IC-scattering: Minimum allowed value of gamma (gamma_{IC, th}) versus final photon energy")
        for i in [0,1]:
            x0=[0.01,0.0025][i]
            Color=['red','blue'][i]
            ColouredSpace = matplotlib.patches.Rectangle((0, 0), 1, 1, fc="%s" % Color, fill=True, edgecolor='none', linewidth=0, alpha=.2)
            print('1/x0 =', 1.0/x0)
            xGammaLowerRange=np.linspace(-2000,-1.0/x0,2001)
            xGammaUpperRange=np.linspace(0.001,2000,2001)
            xGammaTotalRange=np.linspace(-2000,2000,4001)
            ValuesForGammaMinOnxgaL=GammaMinOnxga(x0,xGammaLowerRange)
            ValuesForGammaMinOnxgaU=GammaMinOnxga(x0,xGammaUpperRange)
            Aa=pl.plot(xGammaLowerRange,ValuesForGammaMinOnxgaL, label='$\gamma = \gamma_{\mathrm{IC,\,th}}(x_\gamma,\,x)$', linewidth=1.9, color='%s' % Color)[0]
            pl.plot(xGammaUpperRange,ValuesForGammaMinOnxgaU, linewidth=1.9, color='%s' % Color)
            Bb=pl.plot(xGammaTotalRange,xGammaTotalRange+1.0/(4.0*x0),label='$\gamma = x_\gamma + 1/(4 \, x)$', linewidth=1.9, color='%s' % Color, linestyle='--')[0]
            #Cc=pl.plot(xGammaUpperRange,0.5*np.sqrt(xGammaUpperRange/x0),label='$0.5 \, \sqrt{x_\gamma / x} \,$ (limit for $x_\gamma x \ll 1$)', linewidth=1.9, color='%s' % Color, linestyle=':')[0]
            pl.fill_between(xGammaUpperRange, ValuesForGammaMinOnxgaU, 100000, alpha=.2, color='%s' % Color)
            if i==0:
                Dd=pl.plot(xGammaTotalRange,xGammaTotalRange,label='$\gamma = x_\gamma$', linewidth=1.0, color='gray')[0]
                L1=[Dd, ExtraSpace, ExtraSpace, Aa, Bb, ColouredSpace]
                L2=[Dd.get_label(), '   ', '$x=%.4f:$' % x0, Aa.get_label(), Bb.get_label(), '$\mathrm{Allowed\,Range}$']
            elif i==1:
                L1.extend([ExtraSpace,ExtraSpace, Aa, Bb, ColouredSpace])
                L2.extend(['   ', '$x=%.4f:$' % x0, Aa.get_label(), Bb.get_label(), '$\mathrm{Allowed\,Range}$'])
                pl.legend(L1,L2, loc='best', fontsize=18)
        pl.xlabel('Final gamma-ray photon energy $x_\gamma$', fontsize=26)
        pl.ylabel('Original electron energy $\gamma$', fontsize=26)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=26)
        ax.yaxis.set_tick_params(labelsize=26)
        pl.xlim(-450,299)
        pl.ylim(-299,350)#pl.ylim(xGammaLowerRange[0]+1.0/(4.0*x0),ValuesForGammaMinOnxgaU[-1])
        pl.subplots_adjust(top=0.98, bottom=0.11, left=0.16, right=0.97)
        pl.savefig("IC scattering - gamma_ICth from %s.pdf" % GetCurrentDate())
    
    def EvaluatexmaxOfA20():
        print('1/x0 =', 1.0/x0)
        GammaRange=np.linspace(-2.0/x0,2.0/x0,401)
        GammaRange=GammaRange[GammaRange!=-1.0/(4.0*x0)]
        ValuesForxmaxOfA20=xmaxOfA20(GammaRange,x0)
        pl.figure(figsize=(12, 9), num="IC-scattering: Maximum allowed value of x_gamma versus original electron energy")
        pl.subplot(1,1,1)
        pl.plot(GammaRange,ValuesForxmaxOfA20, label='$x_{\gamma,max}(\gamma)$')
        pl.plot(GammaRange,GammaRange-1.0/(4.0*x0),label='$\gamma - 1/(4 \cdot x_0)$')
        pl.scatter(0, xmaxOfA20(0,x0), 20, color='red')  # Draw a point
        pl.annotate(r'$(0,x_{\gamma,max}(0)=0)$',xy=(0, xmaxOfA20(0,x0)), xycoords='data',xytext=(+10, +50), textcoords='offset points', fontsize=16,arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
        pl.text(1.0/x0,0,r'For $1/(4 \cdot x_0) = %d $' % (1.0/(4.0*x0)), fontsize=16)
        pl.legend(loc="best", fontsize=14)
        pl.xlabel('Original electron energy $\gamma$', fontsize=16)
        pl.ylabel('Final photon energy $x_\gamma$', fontsize=16)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        pl.ylim(1.5*ValuesForxmaxOfA20[0],1.5*ValuesForxmaxOfA20[-1])
        pl.xlim(GammaRange[0],GammaRange[-1])
        pl.subplots_adjust(left=0.2)
        
        ExtraSpace = matplotlib.patches.Rectangle((0, 0), 0.1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
        pl.rc('font', family='serif')
        for i in [0,1]:
            x0=[0.01,0.0025][i]
            Color=['red','blue'][i]
            ColouredSpace = matplotlib.patches.Rectangle((0, 0), 1, 1, fc="%s" % Color, fill=True, edgecolor='none', linewidth=0, alpha=.2)
            print('1/x0 =', 1.0/x0)
            GammaLowerRange=np.linspace(-600,-1.0/(4*x0),1001,endpoint=False)
            GammaUpperRange=np.linspace(-0.9999/(4*x0),1000,1001)
            GammaTotalRange=np.append(GammaLowerRange,GammaUpperRange)
            AllowedGammaRange=GammaTotalRange[GammaTotalRange>0]
            ValuesForxmaxOfA20L=xmaxOfA20(GammaLowerRange,x0)
            ValuesForxmaxOfA20U=xmaxOfA20(GammaUpperRange,x0)
            AllowedValuesForxmaxOfA20=xmaxOfA20(AllowedGammaRange,x0)
            pl.figure(figsize=(9.5, 8.7), num="IC-scattering: Maximum allowed value of x_gamma (x_{gamma max}) versus original electron energy")
            Aa=pl.plot(GammaLowerRange,ValuesForxmaxOfA20L, label='$x_{\gamma} = x_{\gamma,\,\mathrm{max}}(\gamma,\,x)$', linewidth=1.9, color='%s' % Color)[0]
            pl.plot(GammaUpperRange,ValuesForxmaxOfA20U, linewidth=1.9, color='%s' % Color)
            Bb=pl.plot(GammaTotalRange,GammaTotalRange-1.0/(4.0*x0),label='$x_{\gamma} = \gamma - 1/(4 \, x)$', linewidth=1.9, color='%s' % Color, linestyle='--')[0]
            pl.fill_between(AllowedGammaRange, 0, AllowedValuesForxmaxOfA20, alpha=.2, color='%s' % Color)
            pl.scatter(0, xmaxOfA20(0,x0), 30, color='black')  # Draw a point
            pl.annotate(r'$(0,\,x_{\gamma,\,\mathrm{max}}(0,\,x)=0)$',xy=(0, xmaxOfA20(0,x0)), xycoords='data',xytext=(-5, +120), textcoords='offset points', fontsize=16,arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
            if i==0:
                Dd=pl.plot(xGammaTotalRange,xGammaTotalRange,label='$x_\gamma = \gamma$', linewidth=1.0, color='gray')[0]
                L1=[Dd, ExtraSpace, ExtraSpace, Aa, Bb, ColouredSpace]
                L2=[Dd.get_label(), '   ', '$x=%.4f:$' % x0, Aa.get_label(), Bb.get_label(), '$\mathrm{Allowed\,Range}$']
            elif i==1:
                L1.extend([ExtraSpace,ExtraSpace, Aa, Bb, ColouredSpace])
                L2.extend(['   ', '$x=%.4f:$' % x0, Aa.get_label(), Bb.get_label(), '$\mathrm{Allowed\,Range}$'])
                pl.legend(L1,L2, loc='best', fontsize=18)
        pl.xlabel('Original electron energy $\gamma$', fontsize=26)
        pl.ylabel('Final gamma-ray photon energy $x_\gamma$', fontsize=26)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=26)
        ax.yaxis.set_tick_params(labelsize=26)
        pl.ylim(-450,299)
        pl.xlim(-299,350)
        # pl.ylim(1.5*ValuesForxmaxOfA20L[0],1.5*ValuesForxmaxOfA20U[-1])
        # pl.xlim(GammaTotalRange[0],GammaTotalRange[-1])
        pl.subplots_adjust(top=0.98, bottom=0.12, left=0.17, right=0.99)
        pl.savefig("IC scattering - x_gammamax from %s.pdf" % GetCurrentDate())
    
    def EvaluatePartOfIntegrandOfA1(Testga,TestgaP):
        # If one applies this for example for Testga=10**5,TestgaP=20 one sees, that the part of the integrand behind n0 is positive only for values of x >approx 10**(-2) and negative below this.
        # Look at PartOfIntegrandOfA1:
        ValuesForx=np.logspace(-7,-1,100)
        ValuesForPartOfIntegrandOfA1=PartOfIntegrandOfA1(ValuesForx,Testga,TestgaP)
        print(ValuesForPartOfIntegrandOfA1)
        PartOfIntegrandOfA1Max=max(ValuesForPartOfIntegrandOfA1)
        pl.figure(num="Part of integrand of equation A1")
        pl.plot(ValuesForx,ValuesForPartOfIntegrandOfA1, label="$\gamma = %.2f$, $\gamma' = %.2f$" % (Testga, TestgaP))
        pl.scatter(EAstOfA3(Testga,TestgaP)/Testga,0.0) # This is the lowest integration border of COfA1.
        pl.xscale('log')
        pl.ylim(-PartOfIntegrandOfA1Max,PartOfIntegrandOfA1Max)
        pl.xlabel("Soft photon energy $x$")
        pl.legend(loc="upper left")
        # Get the location of the zero: 
        def FuncForRootDetermination(x):
            return PartOfIntegrandOfA1(x,Testga,TestgaP)
        Root=optimize.fsolve(FuncForRootDetermination, 10.0**(-20.0))
        print('Zero at x=%e' % Root)
        # The integrand does behave well for all gaP except gaP=ga. One can show, that the part of the integrand in the squared brackets converges to 1 in the limit gaP-->ga. Thus consider the part in front of the squared brackets. 
        # Look at PartOfIntegrandOfA1InFrontOfBrackets: 
        ValuesForPartOfIntegrandOfA1InFrontOfBrackets=PartOfIntegrandOfA1InFrontOfBrackets(ValuesForx,Testga)
        print(ValuesForPartOfIntegrandOfA1InFrontOfBrackets)
        pl.plot(ValuesForx,ValuesForPartOfIntegrandOfA1InFrontOfBrackets, label='In front of brackets')
        # Compare the two expressions:
        DifferenceOfBothExpressions=ValuesForPartOfIntegrandOfA1InFrontOfBrackets-ValuesForPartOfIntegrandOfA1
        print('Difference of both expressions:\n', DifferenceOfBothExpressions)
        # The more gaP nears ga, the smaller gets the difference. Thus, for gaP=ga, (3.0*sigmaT*c)/(4.0*EOfA2(ga,x)*ga) is a good approximation to the whole part of the integrand behind n0. This is the reason for the if-test in the definition of IntegrandOfA1.
    
    # COfA1 can be plotted in dependence on gaP for a certain value of ga and for a certain soft photon distribution:
    def EvaluateCOfA1():
        LowestPlottingBorder = gaPminOf6(x0,gamma)
        BiggestPlottingBorder = gamma
        ValuesForgammaPLog = np.logspace(np.log10(LowestPlottingBorder),np.log10(BiggestPlottingBorder),NumOfDecadesOfgamma*1000)
        # Firstly, plot the various contributions to COfA1:
        pl.figure(figsize=(12, 9), num="IC-scattering: Spectral probability versus final electron energy")
        ValuesForCOfA1WithoutSynchrotron = np.asarray([COfA1WithoutSynchrotron(gamma,i,Usedn0) for i in ValuesForgammaPLog])
        pl.loglog(ValuesForgammaPLog,ValuesForCOfA1WithoutSynchrotron,label='$\gamma = %4.3g /(4 \cdot x_0)$, Only ExPs' % (gamma*4*x0))
        if IncludeSynchrotron:
            ValuesForCOfA1WithSynchrotron = np.asarray([COfA1WithSynchrotron(gamma,i,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgammaPLog])
            pl.loglog(ValuesForgammaPLog,ValuesForCOfA1WithSynchrotron,label='$\gamma = %4.3g /(4 \cdot x_0)$, Only SyncPs' % (gamma*4*x0))        
            ValuesForCOfA1 = np.asarray([COfA1(gamma,i,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgammaPLog])
            pl.loglog(ValuesForgammaPLog,ValuesForCOfA1,label='$\gamma = %4.3g /(4 \cdot x_0)$, Sum' % (gamma*4*x0))
        pl.loglog(np.ones(2)*1/(4*x0),np.asarray([10**(-30),10**(30)]),label="$\gamma' = 1/(4 \cdot x_0)$")
        PotentialJumps = [] # The function will have kinks, due to additional contributions to COfA1. These additional contributions are caused by discontinuities along x.
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            PotentialJumps = np.append(PotentialJumps,[gaPminOf6(ValuesForxSyncSuperIt[-1],gamma)])
        if Usedn0==n0MultiDelta:
            PotentialJumps = np.append(PotentialJumps,[gaPminOf6(i,gamma) for i in x0MultiDelta])
        else:
            PotentialJumps = np.append(PotentialJumps,[gaPminOf6(x0,gamma)])
        for i in range(len(PotentialJumps)):
            gaP = PotentialJumps[i]
            if LowestPlottingBorder < gaP < BiggestPlottingBorder:
                pl.scatter(gaP, COfA1WithoutSynchrotron(gamma,gaP,Usedn0), 20, color='blue')
        #pl.xscale('log')
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=20)
        ax.yaxis.set_tick_params(labelsize=20)
        pl.legend(loc="best", fontsize=16)
        pl.xlabel("Final electron Lorentz-factor $\gamma'$", fontsize=20)
        pl.ylabel("Spectral IC-scattering probability-rate $C(\gamma,\gamma')$ in $s^{-1}$", fontsize=20)
        #pl.savefig("IC-scattering - Spectral probability versus final electron energy from %s.svg" % GetCurrentDate())
        # Secondly, plot the various contributions to COfA1, normalised to COfA21:
        pl.figure(figsize=(12, 9), num="IC-scattering: Normalised, spectral probability versus final electron energy")
        if not(IncludeSynchrotron):
            ValuesForCOfA1WithoutSynchrotronNormalised = ValuesForCOfA1WithoutSynchrotron/(COfA21(gamma,Usedn0,None,None)*np.asarray([1.0 for i in ValuesForgammaPLog]))
        elif IncludeSynchrotron:
            ValuesForCOfA1WithoutSynchrotronNormalised = ValuesForCOfA1WithoutSynchrotron/(COfA21(gamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt)*np.asarray([1.0 for i in ValuesForgammaPLog]))
            ValuesForCOfA1WithSynchrotronNormalised = ValuesForCOfA1WithSynchrotron/(COfA21(gamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt)*np.asarray([1.0 for i in ValuesForgammaPLog]))
            pl.loglog(ValuesForgammaPLog,ValuesForCOfA1WithSynchrotronNormalised,label='$\gamma = %4.3g /(4 \cdot x_0)$, Only SyncPs' % (gamma*4*x0))        
            ValuesForCOfA1Normalised = ValuesForCOfA1/(COfA21(gamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt)*np.asarray([1.0 for i in ValuesForgammaPLog]))
            pl.loglog(ValuesForgammaPLog,ValuesForCOfA1Normalised,label='$\gamma = %4.3g /(4 \cdot x_0)$, Sum' % (gamma*4*x0))
        pl.loglog(ValuesForgammaPLog,ValuesForCOfA1WithoutSynchrotronNormalised,label='$\gamma = %4.3g /(4 \cdot x_0)$, Only ExPs' % (gamma*4*x0))
        pl.loglog(np.ones(2)*1/(4*x0),np.asarray([10**(-20),10**(20)]),label="$\gamma' = 1/(4 \cdot x_0)$")
        #pl.xscale('log')
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=30)
        ax.yaxis.set_tick_params(labelsize=30)
        pl.legend(loc="best", fontsize=22)
        pl.xlabel("Final electron Lorentz-factor $\gamma'$", fontsize=30)
        pl.ylabel("Normalised, spectral IC-scattering\n probability $C(\gamma,\gamma')$", fontsize=30)
        pl.subplots_adjust(top=0.95, bottom=0.14, left=0.15, right=0.96)
        pl.savefig("IC-scattering - Normalised, spectral probability versus final electron energy from %s.svg" % GetCurrentDate())
    # Obsolete comment: For n0!=n0Delta, the function COfA1WithoutSynchrotron always vanishes for all gammaP below gammaP \approx 1/(4*x0), (almost) independently on the value of gamma (i. e. for reasonably big gamma). The reason is the following: For reasonably big gamma, the lower integration border of COfA1WithoutSynchrotron is approx 1/(4*gammaP), which exceeds the upper border (which in the end is x0) for gammaP smaller than 1/(4*x0). This is in accordance with equation 6.
    # For gammaP bigger than gammaP \approx 1/(4*x1), the effective lower integration border is x1, and thus constant. Hence, the change of the integral with changing gammaP solely depends on the change of the integrand with changing gammaP. So, look at the integrand! The integrand moves towards lower x with increasing gammaP. Thus, the integral decreases with increasing gammaP. The integrand is pathological for gammaP=gamma, so, it is understandable that COfA1WithoutSynchrotron is "not a number" for gammaP=gamma. To prevent this, the if-Test was introduced in the definition of IntegrandOfA1. 
    # Comment of version 6: For gamma=10**4 and x0Delta = 10**(-3) the following is the case for n0 being Dirac-Delta-function-like: For gammaP below approximately 100 the value of the function is extremely negative and above 100 it is positive, but comparatively small in amount. For gamma=10/x0 and PlanckTemperature = 10**6 and x0NonDelta=30*PlanckTheta the following is the case for n0 being Planck- or power-law-like: For gammaP below approximately 45 the value of the function is zero. This site is called "lower border of non-vanishing range of COfA1WithoutSynchrotron" below. Above this site, the function behaves polite and moderate.  
    # Comment of version 9: In the case of gamma=10**4, x0Delta = 10**(-3) and n0 being Dirac-Delta-function-like the zero at about 100 is identical with the minimally allowed gaP. I checked this for other numerical cases, too. So, the region in gaP with negative values of COfA1WithoutSynchrotron is just the region below the allowed range of gaP.
    
    # COfA1 can also be plotted in dependence on ga for a certain value of gaP and for a certain soft photon distribution:
    def EvaluateCOfA1DependingOnga():
        for gammaP in [0.2/x0,2/x0]:
            LowestPlottingBorder = gammaP # This is equal to the LowestIntegrationBorder of the ThirdTerm.
            BiggestPlottingBorder = BiggestIntegrationBorderOfThirdTermOfEq1(ValuesForxSyncSuperIt,gammaP,CurrentLastValuesForgaIt[-1])*0.9999999999 # This is equal to the BiggestIntegrationBorder of the ThirdTerm.
            ValuesForgammaLog = np.logspace(np.log10(LowestPlottingBorder),np.log10(BiggestPlottingBorder),7000)
            pl.figure(figsize=(12, 9), num="IC-scattering: Spectral probability versus original electron energy")
            if not(IncludeSynchrotron):
                ValuesForCOfA1WithoutSynchrotron = [COfA1WithoutSynchrotron(i,gammaP,Usedn0) for i in ValuesForgammaLog]
                #ValuesForCOfA1WithoutSynchrotronNormalised = [COfA1WithoutSynchrotron(i,gammaP,Usedn0)/COfA21(i,Usedn0,None,None) for i in ValuesForgammaLog]        
            elif IncludeSynchrotron:
                ValuesForCOfA1WithoutSynchrotron = [COfA1WithoutSynchrotron(i,gammaP,Usedn0) for i in ValuesForgammaLog]
                #ValuesForCOfA1WithoutSynchrotronNormalised = [COfA1WithoutSynchrotron(i,gammaP,Usedn0)/COfA21(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgammaLog]
                ValuesForCOfA1WithSynchrotron = [COfA1WithSynchrotron(i,gammaP,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgammaLog]
                #ValuesForCOfA1WithSynchrotronNormalised = [COfA1WithSynchrotron(i,gammaP,Usedn0,nSyncPs,ValuesForxSyncSuperIt)/COfA21(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgammaLog]
                pl.loglog(ValuesForgammaLog*x0,ValuesForCOfA1WithSynchrotron,label="$\gamma' \cdot x_0 = %s$, Only SyncPs" % (gammaP*x0))    
                ValuesForCOfA1 = [COfA1(i,gammaP,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgammaLog]
                #ValuesForCOfA1Normalised = [COfA1(i,gammaP,Usedn0,nSyncPs,ValuesForxSyncSuperIt)/COfA21(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgammaLog]
                pl.loglog(ValuesForgammaLog*x0,ValuesForCOfA1,label="$\gamma' \cdot x_0 = %s$, Sum" % (gammaP*x0))        
            pl.loglog(ValuesForgammaLog*x0,ValuesForCOfA1WithoutSynchrotron,label="$\gamma' \cdot x_0 = %s$, Only ExPs" % (gammaP*x0))        
            pl.scatter(gammaP*x0, COfA1WithoutSynchrotron(gammaP,gammaP,Usedn0), 20, color='red')
            PotentialJumps = [] # The function will have kinks, due to additional contributions to COfA1. These additional contributions are caused by discontinuities along x.
            if IncludeSynchrotron and PerformanceMode!='CommonIteration':
                PotentialJumps = np.append(PotentialJumps,[GammaLimit(ValuesForxSyncSuperIt[-1],gammaP)])
            if Usedn0==n0MultiDelta:
                PotentialJumps = np.append(PotentialJumps,[GammaLimit(i,gammaP) for i in x0MultiDelta])
            else:
                PotentialJumps = np.append(PotentialJumps,[GammaLimit(x0,gammaP)])
            for i in range(len(PotentialJumps)):
                ga = PotentialJumps[i]
                if LowestPlottingBorder < ga < BiggestPlottingBorder:
                    if i==1:
                        pl.scatter(ga*x0, COfA1(ga,gammaP,Usedn0,nSyncPs,ValuesForxSyncSuperIt), 20, color='blue')
                    else:
                        pl.scatter(ga*x0, COfA1WithoutSynchrotron(ga,gammaP,Usedn0), 20, color='blue')
        #pl.xscale('log')
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=20)
        ax.yaxis.set_tick_params(labelsize=20)
        pl.legend(loc="best", fontsize=16)
        pl.xlabel("Original electron Lorentz-factor $\gamma$ times $x_0$", fontsize=20)
        pl.ylabel("Spectral IC-scattering probability $C(\gamma,\gamma')$", fontsize=20)
        pl.savefig("IC-scattering - Spectral probability versus original electron energy from %s.svg" % GetCurrentDate())
    
    # COfA1DependingOnxga can be plotted in dependence on xga for various values of ga:
    def EvaluateCOfA1DependingOnxgaForCertainga():
        ValuesForxgaLog = np.logspace(0,np.log10(xmaxOfA20(gamma,x0)),5000)  # A sampling-range for xga. Actually xga ranges from 0 to xmaxOfA20. Logarithmic division of the interval is used.
        pl.figure(figsize=(12, 9), num="IC-scattering: Spectral probability-rate versus final photon energy")
        if not(IncludeSynchrotron):
            ValuesForCOfA1DependingOnxgaForCertainga = [COfA1DependingOnxga(gamma,i,Usedn0,None,None) for i in ValuesForxgaLog]
        elif IncludeSynchrotron:
            ValuesForCOfA1DependingOnxgaForCertainga = [COfA1DependingOnxga(gamma,i,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForxgaLog]
        pl.loglog(ValuesForxgaLog,ValuesForCOfA1DependingOnxgaForCertainga)    
        pl.xlabel("Final photon energy $x_{\gamma}$")
        pl.ylabel("Spectral IC-scattering probability-rate $C(\gamma,\gamma-x_{\gamma})$ in s$^{-1}$")
    
    # COfA1DependingOnxga can be plotted in dependence on ga for various values of xga:
    def EvaluateCOfA1DependingOnxgaForCertainxga():
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            LowestPlottingBorder=min(GammaMinOnxga(x0,xgamma),GammaMinOnxga(ValuesForxSyncSuperIt[-1],xgamma)) # ValuesForxSyncSuperIt[-1] is the current xSync0Used. Essentially, this is the lowest integration border of IntegralOfBracketsOfEq1.
        else:
            LowestPlottingBorder=GammaMinOnxga(x0,xgamma)
        BiggestPlottingBorder=10**3*LowestPlottingBorder
        ValuesForgaLog = np.logspace(np.log10(LowestPlottingBorder),np.log10(BiggestPlottingBorder),500)  # A sampling-range for ga. Actually ga ranges from GammaMinOnxga(x0,xga) to infinity. Logarithmic division of the interval is used.
        pl.figure(figsize=(12, 9), num="IC-scattering: Spectral probability-rate versus initial electron energy")
        ValuesForCOfA1DependingOnxgaForCertainxgaWithoutSynchrotron = [COfA1DependingOnxgaWithoutSynchrotron(i,xgamma,Usedn0) for i in ValuesForgaLog]
        pl.loglog(ValuesForgaLog,ValuesForCOfA1DependingOnxgaForCertainxgaWithoutSynchrotron, label='$x_{\gamma} = %4.3g /{x_0}$, Only ExPs' % (xgamma*x0))
        #gaPoint = 1.0*xgamma+0.004*10**7 # For n0Planck with x0=10^-4 and values of xgamma higher than about 10^5, a sharp maximum is situated here, which becomes sharper with increasing xgamma.
        if not(IncludeSynchrotron):
            ValuesForCOfA1DependingOnxgaForCertainxga = [COfA1DependingOnxga(i,xgamma,Usedn0,None,None) for i in ValuesForgaLog]
            #pl.scatter(gaPoint, COfA1DependingOnxga(gaPoint,xgamma,Usedn0,None,None), 40, color='red')  # Draw a point
        elif IncludeSynchrotron:
            ValuesForCOfA1DependingOnxgaForCertainxgaWithSynchrotron = [COfA1DependingOnxgaWithSynchrotron(i,xgamma,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgaLog]
            pl.loglog(ValuesForgaLog,ValuesForCOfA1DependingOnxgaForCertainxgaWithSynchrotron, label='$x_{\gamma} = %4.3g /{x_0}$, Only SyncPs' % (xgamma*x0))
            ValuesForCOfA1DependingOnxgaForCertainxga = [COfA1DependingOnxga(i,xgamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgaLog]
            #pl.scatter(gaPoint, COfA1DependingOnxga(gaPoint,xgamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt), 40, color='red')  # Draw a point
        pl.loglog(ValuesForgaLog,ValuesForCOfA1DependingOnxgaForCertainxga, label='$x_{\gamma} = %4.3g /{x_0}$, Sum' % (xgamma*x0))
        PotentialJumps = [] # The function will have kinks, due to additional contributions to COfA1. These additional contributions are caused by discontinuities along x.
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            PotentialJumps = np.append(PotentialJumps,[GammaMinOnxga(ValuesForxSyncSuperIt[-1],xgamma)])
        if Usedn0==n0MultiDelta:
            PotentialJumps = np.append(PotentialJumps,[GammaMinOnxga(i,xgamma) for i in x0MultiDelta])
        else:
            PotentialJumps = np.append(PotentialJumps,[GammaMinOnxga(x0,xgamma)])
        for i in range(len(PotentialJumps)):
            ga = PotentialJumps[i]
            if LowestPlottingBorder < ga < BiggestPlottingBorder:
                if i==1:
                    pl.scatter(ga, COfA1DependingOnxga(ga,xgamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt), 20, color='blue')
                else:
                    pl.scatter(ga, COfA1DependingOnxgaWithoutSynchrotron(ga,xgamma,Usedn0), 20, color='blue')
        pl.legend(loc="best", fontsize=16)
        pl.xlabel("Initial electron energy $\gamma$")
        pl.ylabel("Spectral IC-scattering probability-rate $C(\gamma,\gamma-x_{\gamma})$ in s$^{-1}$")
        pl.savefig("IC-scattering - Spectral probability-rate versus initial electron energy from %s.svg" % GetCurrentDate())

    def EvaluatecOfA19AndNormalisedSpectralICScatteringProbability():    
        ValuesForxgaLin = np.linspace(0,xmaxOfA20(gamma,x0),30)  # A sampling-range for xga. Actually xga ranges from 0 to xmaxOfA20. In contrast to above, linear division of the interval is now used.
        if not(IncludeSynchrotron):
            ValuesForcOfA19 = np.array([cOfA19(gamma,i,Usedn0,None,None) for i in ValuesForxgaLin])
            ValuesForNormalisedSpectralICScatteringProbability = np.array([NormalisedSpectralICScatteringProbability(gamma,i,Usedn0,None,None) for i in ValuesForxgaLin])
        elif IncludeSynchrotron:
            ValuesForcOfA19 = np.array([cOfA19(gamma,i,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForxgaLin])
            ValuesForNormalisedSpectralICScatteringProbability = np.array([NormalisedSpectralICScatteringProbability(gamma,i,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForxgaLin])
        pl.figure(figsize=(12, 9), num="IC-scattering: Normalised spectral probability versus final photon energy")
        ScalingFactor = xmaxOfA20(gamma,x0) # This is to rescale the x-axis. However, to conserve the normalisation, one has to rescale the y-axis, too. If one axis is stretched, the other one has to be shrunk, and vice versa. 
        pl.plot(ValuesForxgaLin/ScalingFactor,ValuesForcOfA19*ScalingFactor, label='Normalised solely to total IC-scattering-rate')
        pl.plot(ValuesForxgaLin/ScalingFactor,ValuesForNormalisedSpectralICScatteringProbability*ScalingFactor, label='Normalised to total electron-loss-rate')
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        pl.legend(loc="upper left", fontsize=14)
        pl.xlabel("Final photon energy $x_{\gamma}/x_{max}$", fontsize=16)
        pl.ylabel("Normalised spectral IC-scattering probability $c(\gamma,x_{\gamma}) \cdot x_{max}$", fontsize=16)
        pl.ylim(0,5)

def EvaluateLowerBorderOfNonVanRangeOfCOfA1AndCOfA21AndElectronLossRate(ValuesForxSyncSuperIt,ValuesFornSyncPs,ValuesForgaIt,SuperItCounter):
    print('\nCalling EvaluateLowerBorderOfNonVanRangeOfCOfA1AndCOfA21AndElectronLossRate:')
    if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray:
        nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
    else:
        nSyncPs = None
    ValuesForgammaLog = np.logspace(np.log10(ValuesForgaIt[0]),np.log10(ValuesForgaIt[-1]),61)  # A sampling-range for ga. ga ranges over reasonable values.
    #print('Sample values for gamma:\n', ValuesForgammaLog)
    ValuesForLowerBorderOfNonVanRangeOfCOfA1 = LowerBorderOfNonVanRangeOfCOfA1(x0,ValuesForgammaLog) # Values for the lower border of the non-vanishing range of COfA1.
    #print('Corresponding values for the lower border of the non-vanishing range of COfA1 / minimally possible values of gamma\':\n', ValuesForLowerBorderOfNonVanRangeOfCOfA1)
    pl.rc('font', family='serif')
    Figure = pl.figure(figsize=(18, 10), num="IC-scattering: Probability-rate versus electron energy")
    LeftyAxis = Figure.add_subplot(111)
    # On the left y-axis, the probability-rate of scattering-events is drawn:
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        ValuesForElectronLossRate = [ElectronLossRate(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=0.01) for i in ValuesForgammaLog]
        ValuesForCOfA21 = [COfA21(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=0.01) for i in ValuesForgammaLog]
        ValuesForCOfA21WithSynchrotron = [COfA21WithSynchrotron(i,nSyncPs,ValuesForxSyncSuperIt,RelativeError=0.01) for i in ValuesForgammaLog]
    else:
        ValuesForElectronLossRate = [ElectronLossRate(i,Usedn0,None,None,RelativeError=0.01) for i in ValuesForgammaLog]
        ValuesForCOfA21 = [COfA21(i,Usedn0,None,None,RelativeError=0.01) for i in ValuesForgammaLog]
    ValuesForCOfA21WithoutSynchrotron = [COfA21WithoutSynchrotron(i,Usedn0,RelativeError=0.01) for i in ValuesForgammaLog]
    ValuesForCOfA21Alternative = [2*COfA21Alternative(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=0.01) for i in ValuesForgammaLog]
    #print('Corresponding values for COfA21:\n', ValuesForCOfA21)
    LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForElectronLossRate, label='Spectral loss rate (sum)', linewidth=2, color='blue')
    LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForCOfA21, label='Spectral IC down-scattering rate (soft + sync. photons)', linewidth=2, color='red')
    LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForCOfA21Alternative, label='Spectral IC down-scattering rate (soft + sync. photons, alternative times 2)', linewidth=2, color='red', linestyle='--')
    LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForCOfA21WithoutSynchrotron, label='Spectral IC down-scattering rate (soft photons)', linewidth=2, color='orange')
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForCOfA21WithSynchrotron, label='Spectral IC down-scattering rate (sync. photons)', linewidth=2, color='magenta')        
    if IncludeElectronEscape==True:
        ValuesForElectronEscapeRate = [ElectronEscapeRate(i) for i in ValuesForgammaLog] # The escape-probability-rate.
        LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForElectronEscapeRate, label='Spectral escape rate', linewidth=2, color='green')
    if IncludeSynchrotron==True:
        ValuesForSynchrotronSpectralLossRateAna = [SynchrotronSpectralLossRate(i) for i in ValuesForgammaLog] # The probability-rate for losses due to synchrotron-radiation.
        LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForSynchrotronSpectralLossRateAna, label='Spectral synchrotron loss rate (analytically)', linewidth=2, color='cyan')
        if SuperItCounter:
            ValuesForSynchrotronSpectralLossRateNum = [SynchrotronSpectralLossRate(i,ValuesForgaIt) for i in ValuesForgaIt] # The probability-rate for losses due to synchrotron-radiation.
            LeftyAxis.loglog(ValuesForgaIt*x0,ValuesForSynchrotronSpectralLossRateNum, label='Spectral synchrotron loss rate (numerically)', linewidth=2, color='#00FF33')
    LeftyAxis.xaxis.set_tick_params(labelsize=30, direction='inout', width=3, length=10, pad=7)
    LeftyAxis.yaxis.set_tick_params(labelsize=30, direction='inout', width=3, length=10, right=False)
    pl.legend(loc="best", fontsize=18)
    pl.xlabel("Electron energy times SBP energy $\gamma \cdot x_0$", fontsize=35)
    pl.ylabel("Probability rate in s$^{-1}$", fontsize=35)
    Ticks=LeftyAxis.get_yticks() # Store the locations of the ticks.
    LeftyAxis.set_ylim(Ticks[1],Ticks[-2]) # This setting of the limits might change from plot to plot. Usage for pion-decay: Ticks[1],Ticks[7]
    LeftyAxis.set_xlim(ValuesForgammaLog[0]*x0,ValuesForgammaLog[-1]*x0)
    # The right y-axis shows the optical depth, corresponding to the probability-rate.
    RightyAxis = Figure.add_subplot(111, sharex=LeftyAxis, frameon=False)
    RightyAxis.loglog()
    RightyAxis.set_xlim(ValuesForgammaLog[0]*x0,ValuesForgammaLog[-1]*x0)
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
    pl.savefig("Run %s, Electron-loss-rate - Probability-rate versus original electron energy from %s.svg" % (RunIdentifier,GetCurrentDate()))

def EvaluateDotgamma(ValuesForxSyncSuperIt,ValuesFornSyncPs,ValuesForgaIt):
    print('\nCalling EvaluateDotgamma:')
    if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray:
        nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
    else:
        nSyncPs = None
    ValuesForgammaLog = np.logspace(np.log10(ValuesForgaIt[0]),np.log10(ValuesForgaIt[-1]),61)
    pl.rc('font', family='serif')
    Figure = pl.figure(figsize=(18, 10), num="IC-scattering: Dot gamma versus electron energy")
    LeftyAxis = Figure.add_subplot(111)
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        ValuesForDotgamma = [Dotgamma(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=0.01) for i in ValuesForgammaLog]
        ValuesForDotgammaWithSynchrotron = [DotgammaWithSynchrotron(i,nSyncPs,ValuesForxSyncSuperIt,RelativeError=0.01) for i in ValuesForgammaLog]
    else:
        ValuesForDotgamma = [Dotgamma(i,Usedn0,None,None,RelativeError=0.01) for i in ValuesForgammaLog]
    ValuesForDotgammaWithoutSynchrotron = [DotgammaWithoutSynchrotron(i,Usedn0,RelativeError=0.01) for i in ValuesForgammaLog]
    ValuesForDotgammaAlternative = [2*DotgammaAlternative(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=0.01) for i in ValuesForgammaLog]
    LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForDotgamma, label='IC down-scattering (soft + sync. photons)', linewidth=2, color='red')
    LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForDotgammaAlternative, label='IC down-scattering (soft + sync. photons, alternative times 2)', linewidth=2, color='red', linestyle='--')
    LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForDotgammaWithoutSynchrotron, label='IC down-scattering (soft photons)', linewidth=2, color='orange')
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForDotgammaWithSynchrotron, label='IC down-scattering (sync. photons)', linewidth=2, color='magenta')
    if IncludeElectronEscape==True:
        ValuesForElectronEscapeEnergyLossRate = [i*ElectronEscapeRate(i) for i in ValuesForgammaLog]
        LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForElectronEscapeEnergyLossRate, label='Escape of electrons', linewidth=2, color='green')
    if IncludeSynchrotron==True:
        ValuesForSynchrotronEnergyLossRateAbsoluteValue = [SynchrotronEnergyLossRateAbsoluteValue(i) for i in ValuesForgammaLog]
        LeftyAxis.loglog(ValuesForgammaLog*x0,ValuesForSynchrotronEnergyLossRateAbsoluteValue, label='Synchrotron-radiation', linewidth=2, color='cyan')
    LeftyAxis.xaxis.set_tick_params(labelsize=30, direction='inout', width=3, length=10, pad=7)
    LeftyAxis.yaxis.set_tick_params(labelsize=30, direction='inout', width=3, length=10, right=False)
    pl.legend(loc="best", fontsize=18)
    pl.xlabel("Electron energy times SBP energy $\gamma \cdot x_0$", fontsize=35)
    pl.ylabel("Absolute value $|\dot \gamma|$ of\nLorentz-factor's rate-of-change in s$^{-1}$", fontsize=30)
    pl.subplots_adjust(top=0.96, bottom=0.13, left=0.13, right=0.96)
    pl.savefig("Run %s, Electron-energy-loss-rate - Dot gamma versus original electron energy from %s.svg" % (RunIdentifier,GetCurrentDate()))

# Conclusion of version 6: Most things work well. The shape of all functions agrees with Zdziarski's figure 3 for reasonable values of parameters. This agreement tells, that the program works well. The normalisation of the function cOfA19 works for reasonable parameters. However, it is quite obscure, that the normalisation COfA21WithoutSynchrotron with the Dirac-Delta-function adopts negative values.   

# Documentation of version 7: Naming of variables and functions was changed. 

# Documentation of version 9: The basic change was the change of the lower integration border of the integral COfA21WithoutSynchrotron in the case of n0=n0Dirac. Now the former problem with a negative normalisation is extinct. 


## Considering various pair-production formulae according to 1988ApJ...335..786Z

def EOfB2(xga,xOrxSync):
    return xga*xOrxSync      # Equation B2, this is just a definition, it is in units 1. Here, xga is the incoming high-energy photons's energy. This definition is equally valid for xSync instead of x. So, xOrxSync is either x or xSync.

def EAstOfB3(xga,ga):
    return xga**2.0/(4.0*ga*(xga-ga))      # Equation B3, this is just a definition, it is in units 1. ga is the first pair-produced particle's energy. (xga-ga) is (for negligible x) the second pair-produced particle's energy.

def gammaminOf7(x0OrxSync0Used,xga):     # Equation 7. The minimally possible value of gamma for the process of pair-production. This definition is equally valid for xSync0Used instead of x0. This is called gamma_{PP, min} in the PhD thesis.
    return (xga/2.0)*(1.0 - np.sqrt(1.0 - 1.0/(x0OrxSync0Used*xga)))
    
def gammamaxOf7(x0OrxSync0Used,xga):     # Equation 7. The maximally possible value of gamma for the process of pair-production. For all realistic cases it is gammamaxOf7<xga. This definition is equally valid for xSync0Used instead of x0. This is called gamma_{PP, max} in the PhD thesis.
    return (xga/2.0)*(1.0 + np.sqrt(1.0 - 1.0/(x0OrxSync0Used*xga)))

def xgammaLimitOriginal(x0OrxSync0Used,ga): # The limit for the incident hard photon's energy. For ga<1/(4*x0) (which is however a physically forbidden case) it has to be xga<xgammaLimit, while for ga>1/(4*x0) it has to be xga>xgammaLimit.
# As in pair-production ga and gaP are symmetrical, this relation holds for gaP instead of ga, too.
# This definition is equally valid for xSync0Used instead of x0.
    return ga/(1.0-1.0/(4.0*x0OrxSync0Used*ga))

def xgammaLimit(x0OrxSync0Used,ga):
# This definition is equally valid for xSync0Used instead of x0.
# Remark of version 17: It was realised that the evaluation of SecondTermOfNumeratorOfEq8 takes a lot of time sometimes. This might be due to the kink of FourthTermOfEq1 at ga=(1/2*x0). For ga=(1/2*x0), pOfB18 has a singularity at xga=xgammaLimit(x0,(1/2*x0)), which is the lower integration-border of FourthTermOfEq1. This might be a reason for evaluation problems. To overcome this, the lower integration-border is slightly raised now around ga=(1/2*x0). This is done by multiplication of the original xgammaLimit with (1+Gaussian).
# Remark of v2_4: It was realised that since the inclusion of escape, the denominator of NormalisedSpectralPPProbability does not become = 0 anywhere any more. Thus, NormalisedSpectralPPProbability (which was formerly pOfB18) doesn't have a singularity any more. Hence, this inclusion of a Gaussian seems not necessary any more. Thus, in xgammaLimitClean, xgammaLimitOriginal is used instead of xgammaLimit
    xgammaLimitCoefficient = 0.0001 # The maximum relative deviation from xgammaLimitOriginal.
    xgammaLimitLocation = 1.0/(2.0*x0OrxSync0Used) # The location of the Gaussian. This is the site where pOfB18 has critical behaviour.
    xgammaLimitWidth = 0.05/x0OrxSync0Used 
    return xgammaLimitOriginal(x0OrxSync0Used,ga)*(1.0+xgammaLimitCoefficient*np.exp(-1.0/2.0*((ga-xgammaLimitLocation)/xgammaLimitWidth)**2.0))

def xgammaLimitClean(x0OrxSync0Used,ga):
    '''This definition of xgammaLimit is cleaned insofar that it returns infinity for the kinematically non-allowed range. This is called x_{gamma, th} in the PhD thesis.'''
    if ga>1.0/(4*x0OrxSync0Used):
        return xgammaLimitOriginal(x0OrxSync0Used,ga)
    else:
        return np.inf # This is necessary, because in part5 with inclusion of SyncPs and ExPs, a call min(xgammaLimit(x0,ga),xgammaLimit(xSync0Used,ga)) is called for a ga-value below max(1.0/(4*x0),1.0/(4*xSync0Used)). The min()-call mustn't yield negative values however.

def PartOfIntegrandOfB1(xOrxSync,ga,xga):
    '''This is the part behind n_0 of the integrand of eq. B1 of Zdz. This definition is equally valid for xSync instead of x.'''
    EAstOverE = EAstOfB3(xga,ga)/EOfB2(xga,xOrxSync) # Abbreviation for the subsequent equation.
    return IntegrandOfA1Coefficient/(EOfB2(xga,xOrxSync)*xga) * ( rOfA4(ga,(xga-ga)) - (2.0+rOfA4(ga,(xga-ga)))*(EAstOverE) + 2.0*EAstOverE**2.0 + 2.0*EAstOverE*np.log(1.0/EAstOverE) )

def IntegrandOfB1(x,ga,xga,n0):     # The part of equation B1 behind dx. According to a simple dimension analysis, this is in units 1/s. According to my analysis, it describes the double-spectral probability-rate for pair-production of an electron with energy ga by interaction of a high-energy photon with energy xga with a soft photon background of spectral number-density n0(x). Thus it is a probability per unit time dt, per unit photon-energy dx and per unit final electron-energy dga.
    if ga <= 1.0/(x*4.0): # This range is physically not accessible.
        return 0.0
    elif xga >= xgammaLimitOriginal(x,ga):
        return n0(x) * PartOfIntegrandOfB1(x,ga,xga)
    else:
        return 0.0

def IntegrandOfB1WithSynchrotron(xSync,ga,xga,nSyncPs):     # The analogue to the part of equation B1 behind dx, just with the synchrotron-photon-field nSyncPs instead of n0. This is in units 1/s. It describes the double-spectral probability-rate for pair-production of an electron with energy ga by interaction of a high-energy photon with energy xga with a SyncP-background of spectral number-density nSyncPs(xSync). Thus it is a probability per unit time dt, per unit photon-energy dxSync and per unit final electron-energy dga.
    if ga <= 1.0/(xSync*4.0): # This range is physically not accessible.
        return 0.0
    elif xga >= xgammaLimitOriginal(xSync,ga):
        return nSyncPs(xSync) * PartOfIntegrandOfB1(xSync,ga,xga)
    else:
        return 0.0

# Equation B1. According to a simple dimension analysis, this is in units 1/s. According to my analysis, it describes the spectral probability-rate for events of pair-production of an electron with energy ga by interaction of a high-energy photon with energy xga with a soft photon background of spectral number-density n0(x). Thus it is a probability per unit time dt and per unit final electron-energy dga. The kinematically motivated lower integration border is EAstOfB3(xga,ga)/xga.
if Usedn0 == n0Delta:
    # The definition of x0 was already done above.
    def POfB1WithoutSynchrotron(ga,xga,n0,RelativeError=integrate.quad.__defaults__[3]):   # Defining POfB1WithoutSynchrotron separately for the case of n0 being a Dirac-Delta-function is a trick to bypass the integration via making use of the substitution-feature of the Dirac-Delta-function. 
        return IntegrandOfB1(x0,ga,xga,n0) # The case x0<EAstOfB3(xga,ga)/xga needn't be excluded here via an if-test because it is already treated and set =0 in IntegrandOfB1.
elif Usedn0 == n0MultiDelta:
    def POfB1WithoutSynchrotron(ga,xga,n0,RelativeError=integrate.quad.__defaults__[3]):   # Defining POfB1WithoutSynchrotron separately for the case of n0 being a set of Dirac-Delta-functions again makes use of bypassing the integration via making use of the substitution-feature. This is done for each Delta-peak and then the sum of them is determined.
        POfB1WithoutSynchrotronSum = 0.0
        for i in x0MultiDelta:
            POfB1WithoutSynchrotronSum += IntegrandOfB1(i,ga,xga,n0) # The case x0<EAstOfB3(xga,ga)/xga needn't be excluded here via an if-test because it is already treated and set =0 in IntegrandOfB1.
        return POfB1WithoutSynchrotronSum
elif Usedn0==n0Exp or Usedn0==n0PL or Usedn0==n0Planck:  # In the following cases, the integration has to be performed numerically. Again, the definition of x0 was already done above.
    def POfB1WithoutSynchrotron(ga,xga,n0,RelativeError=integrate.quad.__defaults__[3]):
        LowestIntegrationBorder = max(x1,EAstOfB3(xga,ga)/xga)
        if LowestIntegrationBorder >= x0:
            return 0.0
        else:
            POfB1WithoutSynchrotronResult = integrate.quad(IntegrandOfB1, LowestIntegrationBorder, x0, args=(ga,xga,n0), epsrel=RelativeError)[0]
            return POfB1WithoutSynchrotronResult
    # Comment of version 16: The if-test has been changed such that it immediately returns 0 in the case LowestIntegrationBorder >= x0.
elif Usedn0==n0AF: # Again, the definition of x0 was already done above. This case is still experimental.
    def POfB1WithoutSynchrotron(ga,xga,n0,RelativeError=integrate.quad.__defaults__[3]):
        LowestIntegrationBorder = max(ADAFx1,EAstOfB3(xga,ga)/xga)
        BiggestIntegrationBorder = x0
        if LowestIntegrationBorder >= BiggestIntegrationBorder:
            return 0.0
        else:
            NumberOfBorders = max(6,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
            ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
            for PotentialBorder in [ADAFxMin,ADAFxPeak,ADAFxComptonisedCutOff]:
                if LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder and (PotentialBorder not in ListOfIntegrationBorders):
                    ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
            ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
            POfB1WithoutSynchrotronResult = 0.0
            for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                POfB1WithoutSynchrotronResult += integrate.quad(IntegrandOfB1, LeftBorder, RightBorder, args=(ga,xga,n0), epsrel=RelativeError)[0]
            return POfB1WithoutSynchrotronResult
# Comment of v2_3: Was renamed from POfB1 to POfB1WithoutSynchrotron.

def POfB1WithSynchrotron(ga,xga,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the analogue to equation B1 just with SyncPs as soft photons. It is in units 1/s and describes the spectral probability-rate for events of pair-production of an electron with energy ga by interaction of a high-energy photon with energy xga with a synchrotron-photon-background of spectral number-density nSyncPs(xSync). Thus it is a probability per unit time dt and per unit final electron-energy dga. The kinematically motivated lower integration border is EAstOfB3(xga,ga)/xga. However, nSyncPs is non-vanishing only in between of ValuesForxSyncSuperIt, hence one needs to integrate only along this range.'''
    LowestIntegrationBorder = max(ValuesForxSyncSuperIt[0],EAstOfB3(xga,ga)/xga) # Essentially, the first argument is CurrentxSync1.
    BiggestIntegrationBorder = ValuesForxSyncSuperIt[-1] # Essentially, this is CurrentxSync0Used and this is the (physically) correct upper border.
    if LowestIntegrationBorder >= BiggestIntegrationBorder:
        return 0.0
    else:
        if LowestIntegrationBorder<0.01*BiggestIntegrationBorder: # Comment of v2_5: By evaluating EvaluateLookAtPartOfIntegrandOfB1 for reasonable parameters of xga and ga, one can recognise that PartOfIntegrandOfB1 always has its peak slightly right of the lower border (EAstOfB3(xga,ga)/xga) of its non-vanishing range. Above the peak, it is a power-law with slope of approx -1. If PartOfIntegrandOfB1 is multiplied with nSyncPs, the slope of the IntegrandOfB1 is even reduced. Therefore, for those cases where the peak of PartOfIntegrandOfB1 is low enough, it is not necessary to extend the integration range up to ValuesForxSyncSuperIt[-1]. Hence, it is chosen that if the lower border is lower than one tenth of the upper border, then the upper border is reduced by factor 0.5 and if the lower border is lower than one hundredth of the upper border, then the upper border is reduced by factor 0.1.
            NumericalBiggestIntegrationBorder=0.1*BiggestIntegrationBorder
            ListOfIntegrationBorders = np.asarray([LowestIntegrationBorder,NumericalBiggestIntegrationBorder])
        elif LowestIntegrationBorder<0.1*BiggestIntegrationBorder:
            NumericalBiggestIntegrationBorder=0.5*BiggestIntegrationBorder
            ListOfIntegrationBorders = np.asarray([LowestIntegrationBorder,NumericalBiggestIntegrationBorder])
        else:
            NumericalBiggestIntegrationBorder=BiggestIntegrationBorder
            ListOfIntegrationBorders = np.asarray([LowestIntegrationBorder,NumericalBiggestIntegrationBorder])
        for PotentialBorder in ValuesForxSyncSuperIt[1:-1]: # At the elements of ValuesForxSyncSuperIt, nSyncPs has kinks.
            if LowestIntegrationBorder<PotentialBorder<NumericalBiggestIntegrationBorder:
                ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
        ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
        POfB1WithSynchrotronResult = 0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            POfB1WithSynchrotronResult += integrate.quad(IntegrandOfB1, LeftBorder, RightBorder, args=(ga,xga,nSyncPs), epsrel=RelativeError)[0]
        return POfB1WithSynchrotronResult

if IncludeSynchrotron and PerformanceMode!='CommonIteration':
    def POfB1(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
        return POfB1WithoutSynchrotron(ga,xga,n0,RelativeError)+POfB1WithSynchrotron(ga,xga,nSyncPs,ValuesForxSyncSuperIt,RelativeError)
else:
    def POfB1(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]): # For the usage of POfB1 (and all functions that are built on it) in the case without synchrotron-back-reaction, one should just use nSyncPs=None and ValuesForxSyncSuperIt=None.
        return POfB1WithoutSynchrotron(ga,xga,n0,RelativeError)

def POfB11WithoutSynchrotron(xga,n0,RelativeError=integrate.quad.__defaults__[3]):     # Equation B11. According to a simple dimension analysis, this is in units 1/s, too. According to my analysis, it describes the probability-rate for pair-production events, in which a particle with any allowed energy ga is produced via interaction of a high-energy photon with energy xga with a soft photon background of spectral number-density n0. Thus it is a probability per unit time dt.
# Actually the integration range is not given. It would be reasonable to integrate over the allowed range of ga, which is going from gammaminOf7 to gammamaxOf7.
    if xga<1.0/x0: # At this set-in of pair-production, it is gammaminOf7(x0,xga)=gammamaxOf7(x0,xga)=1/(2*x0).
        return 0.0
    else:
        LowestIntegrationBorder=gammaminOf7(x0,xga)
        BiggestIntegrationBorder=gammamaxOf7(x0,xga)
        ListOfIntegrationBorders = SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
        POfB11WithoutSynchrotronResult=0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            POfB11WithoutSynchrotronResult += integrate.quad(POfB1WithoutSynchrotron, LeftBorder, RightBorder, args=(xga,n0,RelativeError), epsrel=RelativeError)[0]
        return POfB11WithoutSynchrotronResult
    # In version 13 the division of the integration-range was added but it does not reduce the error messages.

def POfB11WithSynchrotron(xga,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''Analogue to equation B11. This is in units 1/s, too. It describes the probability-rate for pair-production events, in which a particle with any allowed energy ga is produced via interaction of a high-energy photon with energy xga with a soft SyncP-background of spectral number-density nSyncPs.'''
# Actually the integration range is not given. It would be reasonable to integrate over the allowed range of ga, which is going from gammaminOf7 to gammamaxOf7.
    CurrentxSync0Used = ValuesForxSyncSuperIt[-1] # Recreate the current xSync0Used.
    if xga<1.0/CurrentxSync0Used:
        return 0.0
    else:
        LowestIntegrationBorder=gammaminOf7(CurrentxSync0Used,xga)
        BiggestIntegrationBorder=gammamaxOf7(CurrentxSync0Used,xga)
        ListOfIntegrationBorders = SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
        POfB11WithSynchrotronResult=0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            POfB11WithSynchrotronResult += integrate.quad(POfB1WithSynchrotron, LeftBorder, RightBorder, args=(xga,nSyncPs,ValuesForxSyncSuperIt,RelativeError), epsrel=RelativeError)[0]
        return POfB11WithSynchrotronResult
    # In version 13 the division of the integration-range was added but it does not reduce the error messages.

def POfB11(xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
    '''Analogue to equation B11. This is in units 1/s, too. According to my analysis, it describes the probability-rate for pair-production events, in which a particle with any allowed energy ga is produced via interaction of a high-energy photon with energy xga with a soft photon background of spectral number-density n0 (plus a soft SyncP-background nSyncPs). Thus it is a probability per unit time dt.'''
# Actually the integration range is not given. It would be reasonable to integrate over the allowed range of ga, which is going from gammaminOf7 to gammamaxOf7.
# Again, this is valid only if the two summands of the integrand vanish outside of their definition range.
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        CurrentxTotal0 = DeterminexTotal0(ValuesForxSyncSuperIt[-1]) # Recreate the current xTotal0.
    else:
        CurrentxTotal0 = x0
    if xga<1.0/CurrentxTotal0:
        return 0.0
    else:
        LowestIntegrationBorder=gammaminOf7(CurrentxTotal0,xga)
        BiggestIntegrationBorder=gammamaxOf7(CurrentxTotal0,xga)
        ListOfIntegrationBorders = SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder)[0]
        POfB11Result=0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            POfB11Result += integrate.quad(POfB1, LeftBorder, RightBorder, args=(xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError), epsrel=RelativeError)[0]
        return POfB11Result
    # In version 13 the division of the integration-range was added but it does not reduce the error messages.

# This is an alternative implementation of equation B11 with the used soft photons (hence soft background photons or soft photons plus SyncPs). This is physically the same as POfB11, however, it does not integrate POfB1, but it is the sum of both part-integrals POfB11WithoutSynchrotron and POfB11WithSynchrotron.
if IncludeSynchrotron and PerformanceMode!='CommonIteration':
    def POfB11Alternative(xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
        return POfB11WithoutSynchrotron(xga,n0,RelativeError)+POfB11WithSynchrotron(xga,nSyncPs,ValuesForxSyncSuperIt,RelativeError)
else:
    def POfB11Alternative(xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):
        return POfB11WithoutSynchrotron(xga,n0,RelativeError)

def HEPhotonLossRate(xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]): # This is in units of 1/s and gives the spectral probability-rate of events that effectively destroy/remove photons. It includes destruction via pair-production (included in every case for the case of soft background photons and, depending on the value of IncludeSynchrotron, also for the case of a SyncP-background) as well as via escape from the interaction region.
    return POfB11(xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)+HEPhotonEscapeRate(xga)

def pOfB18(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):        # Equation B18. According to my analysis, this describes a normalised spectral probability for events of pair-production of a particle with energy ga (and an anti-particle with energy gaP) by interaction of a high-energy photon of energy xga with a soft photon background of spectral number-density n0 (plus nSyncPs). Thus it is a probability per unit final particle energy dga.
    return 2.0*POfB1(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)/POfB11(xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def NormalisedSpectralPPProbability(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=integrate.quad.__defaults__[3]):        # This describes a normalised spectral probability for events of pair-production of a particle with energy ga (and an anti-particle with energy gaP) by interaction of a high-energy photon of energy xga with a soft photon background of spectral number-density n0 (plus nSyncPs). Thus it is a probability per unit final particle energy dga. In comparison to pOfB18, it is normalised to the total photon-loss-rate, not only to the pair-production-rate.
    return 2.0*POfB1(ga,xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)/HEPhotonLossRate(xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def OpticalDepthOfHEPhotons(Process,*PositionalArgumentsOfProcess):
    '''This is the optical depth of the high-energetic photons with respect to the process Process, which is either pair-production on the soft background-photons (Process=POfB11, arguments have to be xga and n0), or escape (Process=HEPhotonEscapeRate, arguments have to be xga) or both (Process=HEPhotonLossRate, arguments have to be xga and n0). It is dimensionless.
    It was yielded based on the relation loss-probability-rate-of-process = c*dtau/dlength, where tau is the optical depth and length the path-length. Solving for dtau and integrating yields tau = loss-probability-rate-of-process*MeanEscapeLength/c, where MeanEscapeLength is the mean distance a photon has to pass until it leaves the interaction region. Using c = MeanEscapeLength/MeanEscapeTime gives the following expression.'''
    return Process(*PositionalArgumentsOfProcess)*MeanEscapeTime # This essentially is the equation tau = loss-probability-rate-of-process * MeanEscapeTime.

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.

    def EvaluategammaminmaxOf7():
        ExtraSpace = matplotlib.patches.Rectangle((0, 0), 0.1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
        pl.rc('font', family='serif')
        for i in [0,1]:
            x0=[0.01,0.0025][i]
            Color=['red','blue'][i]
            ColouredSpace = matplotlib.patches.Rectangle((0, 0), 1, 1, fc="%s" % Color, fill=True, edgecolor='none', linewidth=0, alpha=.2)
            #xGammaRange=np.linspace(1.0/x0,3.0/x0,1001)
            xGammaRange=np.linspace(1.0/x0,2000,3001)
            xGammaTotalRange=np.linspace(0,2000,3001)
            ValuesForgammaminOf7=gammaminOf7(x0,xGammaRange)
            ValuesForgammamaxOf7=gammamaxOf7(x0,xGammaRange)
            pl.figure(figsize=(10, 10), num="Pair-production: Limits of gamma versus incident photon energy")
            Aa=pl.plot(xGammaRange,ValuesForgammaminOf7, label='$\gamma = \gamma_{\mathrm{PP,\,min}}(x_\gamma,\,x)$', linewidth=2.2, color='%s' % Color)[0]
            Bb=pl.plot(xGammaRange,ValuesForgammamaxOf7, label='$\gamma = \gamma_{\mathrm{PP,\,max}}(x_\gamma,\,x)$', linewidth=2.5, linestyle='--', color='%s' % Color)[0]
            Cc=pl.plot(xGammaTotalRange,np.ones(len(xGammaTotalRange))*1.0/(4.0*x0),label='$\gamma = 1/(4 \, x)$', linewidth=3, linestyle=':', color='%s' % Color)[0]
            Dd=pl.plot(xGammaTotalRange,xGammaTotalRange-1.0/(4.0*x0),label='$\gamma = x_\gamma - 1/(4 \, x)$', linewidth=2.5, linestyle='-.', color='%s' % Color)[0]
            pl.fill_between(xGammaRange, ValuesForgammaminOf7, ValuesForgammamaxOf7, alpha=.2, color='%s' % Color)
            if i==0:
                Ee=pl.plot(xGammaTotalRange,xGammaTotalRange,label='$\gamma = x_\gamma$', linewidth=1.0, color='gray')[0]
                L1=[Ee, ExtraSpace, ExtraSpace, Aa, Bb, Cc, Dd, ColouredSpace]
                L2=[Ee.get_label(), '   ', '$x=%.4f:$' % x0, Aa.get_label(), Bb.get_label(), Cc.get_label(), Dd.get_label(), 'Allowed Range']
            elif i==1:
                L1.extend([ExtraSpace,ExtraSpace, Aa, Bb, Cc, Dd, ColouredSpace])
                L2.extend(['   ', '$x=%.4f:$' % x0, Aa.get_label(), Bb.get_label(), Cc.get_label(), Dd.get_label(), 'Allowed Range'])
                pl.legend(L1,L2, loc='best', fontsize=20)
            pl.xlabel('Gamma-ray photon energy $x_\gamma$', fontsize=26)
            pl.ylabel('Electron energy $\gamma$', fontsize=26)
            pl.xlim(0,799)
            pl.ylim(0,799)
            ax = pl.gca()
            ax.xaxis.set_tick_params(labelsize=26)
            ax.yaxis.set_tick_params(labelsize=26)
            pl.subplots_adjust(top=0.99, bottom=0.10, left=0.13, right=0.99)
        pl.savefig("Pair-production - gammaPPminmax from %s.pdf" % GetCurrentDate())
    
    def EvaluatexgammaLimit():
        Figure = pl.figure(figsize=(12, 16), num="Pair-production: Minimum allowed value of x_gamma versus final particle energy")
        Plot1 = pl.subplot2grid((4, 1), (0, 0), rowspan=3)
        Plot2 = pl.subplot2grid((4, 1), (3, 0))
        for x0 in [0.02,0.01]:
            print('1/(4*x0) =', 1.0/(4.0*x0))
            GammaRange=np.linspace(-4.0*1.0/(4.0*x0),5.0*1.0/(4.0*x0),1001)
            GammaRange=GammaRange[GammaRange!=1.0/(4.0*x0)]
            ValuesForxgammaLimitOriginal=xgammaLimitOriginal(x0,GammaRange)
            ValuesForxgammaLimit=xgammaLimit(x0,GammaRange)
            ValuesForRelativeDeviation = (ValuesForxgammaLimit-ValuesForxgammaLimitOriginal)/ValuesForxgammaLimitOriginal
            Plot1.plot(GammaRange,ValuesForxgammaLimitOriginal,linestyle='None',marker='.',markersize=1, label='$x_{\gamma,min}(\gamma)$')
            Plot1.plot(GammaRange,ValuesForxgammaLimit,linestyle='None',marker='.',markersize=1, label='$x_{\gamma,min}(\gamma)$ with Gaussian deviation')
            Plot1.plot(GammaRange,GammaRange+1.0/(4.0*x0),label='$\gamma + 1/(4 \cdot x_0)$')
            Plot1.legend()
            pl.legend(loc="best", fontsize=14)
            Plot1.set_ylabel('Original photon energy $x_\gamma$', fontsize=16)
            Plot1.set_ylim(2.0*ValuesForxgammaLimit[0],2.0*ValuesForxgammaLimit[-1])
            Plot1.set_xticklabels([])
            Plot1.set_xlim(GammaRange[0],GammaRange[-1])
            Plot2.plot(GammaRange,ValuesForRelativeDeviation,linestyle='None',marker='.',markersize=1, label='$x_0 = %f$' %x0)
            Plot2.set_xlabel('Final particle energy $\gamma$', fontsize=16)
            Plot2.set_ylabel('Relative deviation\nin photon energy $x_\gamma$', fontsize=16)
            pl.legend(loc="best", fontsize=14)
            Plot2.set_xlim(GammaRange[0],GammaRange[-1])
        pl.subplots_adjust(hspace=0.0)
        
        pl.figure(figsize=(10, 10), num="Pair-production: Minimum allowed value of x_gamma versus final particle energy")
        pl.rc('font', family='serif')
        for i in [0,1]:
            x0=[0.01,0.0025][i]
            Color=['red','blue'][i]
            ColouredSpace = matplotlib.patches.Rectangle((0, 0), 1, 1, fc="%s" % Color, fill=True, edgecolor='none', linewidth=0, alpha=.2)
            print('1/(4*x0) =', 1.0/(4.0*x0))
            #GammaLowerRange=np.linspace(-4.0*1.0/(4.0*x0),1.0/(4.0*x0),1001,endpoint=False)
            GammaLowerRange=np.linspace(-1000,1.0/(4.0*x0),2001,endpoint=False)
            GammaUpperRange=np.linspace(1.0001/(4.0*x0),1000,5001)
            GammaTotalRange=np.append(GammaLowerRange,GammaUpperRange)
            ValuesForxgammaLimitOriginalL=xgammaLimitOriginal(x0,GammaLowerRange)
            ValuesForxgammaLimitOriginalU=xgammaLimitOriginal(x0,GammaUpperRange)
            Aa=pl.plot(GammaLowerRange,ValuesForxgammaLimitOriginalL, label='$x_{\gamma} = x_{\gamma,\,\mathrm{th}}(\gamma,\,x)$', linewidth=1.9, color='%s' % Color)[0]
            pl.plot(GammaUpperRange,ValuesForxgammaLimitOriginalU, color='%s' % Color, linewidth=1.9)
            pl.fill_between(GammaUpperRange, ValuesForxgammaLimitOriginalU, 10000, alpha=.2, color='%s' % Color)
            Bb=pl.plot([1/(4*x0),1/(4*x0)],[-10000,10000],label="$\gamma  = 1/(4 \, x)$", linestyle=':', color='%s' % Color, linewidth=3)[0]
            Cc=pl.plot(GammaTotalRange,GammaTotalRange+1.0/(4.0*x0),label='$x_{\gamma} = \gamma + 1/(4 \, x)$', color='%s' % Color, linestyle='-.', linewidth=2.4)[0]
            pl.scatter(1.0/(2.0*x0),1/x0, color='%s' % Color, s=40)
            if i==0:
                pl.annotate(r'$\left( \frac{1}{2 x},\,\frac{1}{x} \right)$',xy=(1.0/(2.0*x0),1/x0), xycoords='data',xytext=(-120, +10), textcoords='offset points', fontsize=20,arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
                Ee=pl.plot(GammaTotalRange,GammaTotalRange,label='$x_\gamma = \gamma$', linewidth=1.0, color='gray')[0]
                L1=[Ee, ExtraSpace, ExtraSpace, Aa, Bb, Cc, ColouredSpace]
                L2=[Ee.get_label(), '   ', '$x=%.4f:$' % x0, Aa.get_label(), Bb.get_label(), Cc.get_label(), 'Allowed Range']
            elif i==1:
                pl.annotate(r'$\left( \frac{1}{2 x},\,\frac{1}{x} \right)$',xy=(1.0/(2.0*x0),1/x0), xycoords='data',xytext=(-5, +70), textcoords='offset points', fontsize=20,arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
                L1.extend([ExtraSpace,ExtraSpace, Aa, Bb, Cc, ColouredSpace])
                L2.extend(['   ', '$x=%.4f:$' % x0, Aa.get_label(), Bb.get_label(), Cc.get_label(), 'Allowed Range'])
                pl.legend(L1,L2, loc='best', fontsize=20)
        pl.ylabel('Gamma-ray photon energy $x_\gamma$', fontsize=26)
        # pl.ylim(2.0*ValuesForxgammaLimitOriginalL[0],2.0*ValuesForxgammaLimitOriginalU[-1])
        # pl.xlim(GammaTotalRange[0],GammaTotalRange[-1])
        pl.xlim(-201,799)
        pl.ylim(-199,799)
        pl.xlabel('Electron energy $\gamma$', fontsize=26)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=26)
        ax.yaxis.set_tick_params(labelsize=26)
        pl.subplots_adjust(top=0.99, bottom=0.09, left=0.14, right=0.99)
        pl.savefig("Pair-production - xgammath from %s.pdf" % GetCurrentDate())

    def EvaluateLookAtPartOfIntegrandOfB1(Testga,Testxga):
        # Look at the part of the integrand behind n0. If one applies this for example for Testga=10**5,Testxga=2*10**5 in the range 2*10**(-6)<x<10**(-3) one can discern the real course of the function. For increasing Testga-->Testxga the trough moves to the bottom right, but the overall course of the function remains unchanged. For Testga=Testxga the function cannot be evaluated due to zero division.
        ValuesForx=np.logspace(-5.7,-3,100);ValuesForPartOfIntegrandOfB1=PartOfIntegrandOfB1(ValuesForx,Testga,Testxga)
        print(ValuesForPartOfIntegrandOfB1)
        PartOfIntegrandOfB1Max=max(ValuesForPartOfIntegrandOfB1)
        PartOfIntegrandOfB1Min=min(ValuesForPartOfIntegrandOfB1)
        pl.figure(num="Part of integrand of equation B1")
        pl.plot(ValuesForx,ValuesForPartOfIntegrandOfB1, label="$\gamma = %.2f$, $x_{\gamma} = %.2f$" % (Testga, Testxga))
        pl.xscale('log')
        pl.ylim(PartOfIntegrandOfB1Min,PartOfIntegrandOfB1Max)
        pl.xlabel("Soft photon energy $x$")
        pl.legend(loc="lower right")
        # The integrand does behave well for all ga except ga=xga. However one can show, that ga-->xga is possible only for x*xga-->infinity, which is not realised in reality.
    
    # POfB1 can be plotted in dependence on ga for a certain value of xga and for a certain soft photon distribution:
    def EvaluatePOfB1():
        ValuesForgammaLin = np.linspace(gammaminOf7(x0,xgamma),gammamaxOf7(x0,xgamma),1000)  # A sampling-range for ga. Linear division of the interval is used.
        pl.figure(figsize=(12, 9), num="Pair-production: Spectral probability-rate versus electron energy")
        ValuesForPOfB1WithoutSynchrotron = [POfB1WithoutSynchrotron(i,xgamma,Usedn0) for i in ValuesForgammaLin]
        pl.plot(ValuesForgammaLin,ValuesForPOfB1WithoutSynchrotron, label='$x_{\gamma} = %s$\nOnly soft background photons' % xgamma)
        if IncludeSynchrotron:
            ValuesForPOfB1WithSynchrotron = [POfB1WithSynchrotron(i,xgamma,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgammaLin]
            pl.plot(ValuesForgammaLin,ValuesForPOfB1WithSynchrotron, label='$x_{\gamma} = %s$\nOnly synchrotron-photons' % xgamma)
            ValuesForPOfB1 = [POfB1(i,xgamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgammaLin]
            pl.plot(ValuesForgammaLin,ValuesForPOfB1, label='$x_{\gamma} = %s$\nSum' % xgamma)
        #pl.yscale('log')
        pl.legend(loc="upper right")
        pl.xlabel("Electron energy $\gamma$")
        pl.ylabel("Spectral pair-production probability-rate $P(x_{\gamma},\gamma)$ in s$^{-1}$")
    # For all background photon-fields, the function POfB1WithSynchrotron is symmetric (tooth-like) about xgamma/2. Especially one can discern that gammaminOf7 is equal to xgamma-gammamaxOf7. Furthermore, I could confirm that for n0=n0Delta, the function POfB1WithSynchrotron is peaked at 1/(2*x0) and at xga-1/(2*x0).
    # For n0=n0Planck, the graph of POfB1WithSynchrotron changes as follows for increasing xgamma: The peaks become more pronounced. The graph is stretched in direction of the x-axis and compressed along the y-axis. For n0=n0Delta, the graph of POfB1WithSynchrotron changes as follows for increasing xgamma: The peaks become sharper. The graph is stretched in direction of the x-axis and compressed along the y-axis. For n0=n0PL, the graph of POfB1WithSynchrotron is as follows: For low xgamma the graph is just a bump. For increasing xgamma, the following happens: The peaks become sharper and more pronounced while the graph is stretched in direction of the x-axis and compressed along the y-axis. For big xgamma it is again a double peaked tooth-like structure. For n0=n0Exp, the graph of POfB1WithSynchrotron is as follows: For very low xgamma the graph is just a bump. For increasing xgamma, the following happens: The peaks become sharper and more pronounced while the graph is stretched in direction of the x-axis and compressed along the y-axis. For big xgamma it is again a double peaked tooth-like structure.

    def EvaluatePOfB1InDependenceOnxga():
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            MinimumxgaLog=np.log10(min(xgammaLimitClean(x0,gamma),xgammaLimitClean(ValuesForxSyncSuperIt[-1],gamma))) # This is the lower integration border of 4th term.
        else: # This case has to be treated separately, because ValuesForxSyncSuperIt has the value None, which cannot be sliced.
            MinimumxgaLog=np.log10(xgammaLimitClean(x0,gamma))
        MaximumxgaLog = np.log10(DetermineBiggestIntegrationBorderOf4thTermOfEq1(ValuesForxSyncSuperIt,CurrentLastValuesForgaIt[-1]))
        FirstDivisionxgaLog = 1.0050*MinimumxgaLog
        SecondDivisionxgaLog = 0.990*MaximumxgaLog
        ValuesForxgaLog = np.concatenate((np.logspace(MinimumxgaLog,FirstDivisionxgaLog,3000,endpoint=False),np.logspace(FirstDivisionxgaLog,SecondDivisionxgaLog,1000,endpoint=False),np.logspace(SecondDivisionxgaLog,MaximumxgaLog,200)))  # A sampling-range for xga.        
        ValuesForPOfB1WithoutSynchrotron = [POfB1WithoutSynchrotron(gamma,i,Usedn0) for i in ValuesForxgaLog]
        if IncludeSynchrotron:
            ValuesForPOfB1WithSynchrotron = [POfB1WithSynchrotron(gamma,i,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForxgaLog]
            ValuesForPOfB1 = [POfB1(gamma,i,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForxgaLog]
        pl.figure(figsize=(12, 9), num="Pair-production: Spectral probability-rate versus HEP energy")
        pl.loglog(ValuesForxgaLog,ValuesForPOfB1WithoutSynchrotron, label='Only soft background photons')
        if IncludeSynchrotron:
            pl.loglog(ValuesForxgaLog,ValuesForPOfB1WithSynchrotron, label='Only synchrotron-photons')
            pl.loglog(ValuesForxgaLog,ValuesForPOfB1, label='Sum')
        pl.legend(loc="best")
        pl.xlabel("HEP energy $x_{\gamma}$")
        pl.ylabel("Spectral pair-production probability-rate $P(x_{\gamma},\gamma)$ in s$^{-1}$")
        pl.savefig("Pair-production - Spectral probability-rate versus HEP energy from %s.svg" % GetCurrentDate())

    def EvaluatepOfB18AndNormalisedSpectralPPProbability():
        ValuesForgaLin = np.linspace(gammaminOf7(x0,xgamma),gammamaxOf7(x0,xgamma),6000)  # A sampling-range for ga. Actually ga ranges from gammaminOf7 to gammamaxOf7, but here, because of the symmetry, only the left half is of interest. Linear division of the interval is used.
        if not(IncludeSynchrotron):
            ValuesForpOfB18 = np.array([pOfB18(i,xgamma,Usedn0,None,None) for i in ValuesForgaLin])
            #ValuesForNormalisedSpectralPPProbability = np.array([NormalisedSpectralPPProbability(i,xgamma,Usedn0,None,None) for i in ValuesForgaLin])
        elif IncludeSynchrotron:
            ValuesForpOfB18 = np.array([pOfB18(i,xgamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgaLin])
            #ValuesForNormalisedSpectralPPProbability = np.array([NormalisedSpectralPPProbability(i,xgamma,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForgaLin])
        pl.figure(figsize=(12, 9), num="Pair-production: Normalised spectral probability versus final particle energy")
        ScalingFactor = xgamma # This is to rescale the x-axis. However, to conserve the normalisation, one has to rescale the y-axis, too. If one axis is stretched, the other one has to be shrunk, and vice versa. 
        pl.plot(ValuesForgaLin/ScalingFactor,ValuesForpOfB18*ScalingFactor, label='$x_{\gamma} = %3.4g/x_0$' % (xgamma*x0))#Normalised solely to pair-production
        #pl.plot(ValuesForgaLin/ScalingFactor,ValuesForNormalisedSpectralPPProbability*ScalingFactor, label='Normalised to total loss rate')
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=26)
        ax.yaxis.set_tick_params(labelsize=26)
        pl.legend(loc="best", fontsize=24)
        pl.xlabel("Final electron energy divided by gamma ray energy $\gamma/x_{\gamma}$", fontsize=26)
        pl.ylabel("Normalised spectral pair-production probability\n times gamma ray energy $p(x_{\gamma},\gamma) \cdot x_{\gamma}$", fontsize=26)
        #pl.ylim(0,5)
        #pl.savefig("Physik/Astrophysik - AGN/Pair-production - Normalised spectral probability versus final particle energy.png", dpi=72)
    
    def EvaluatepOfB18AndNormalisedSpectralPPProbabilityInDependenceOnxga():
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            MinimumxgaLog=np.log10(min(xgammaLimitClean(x0,gamma),xgammaLimitClean(ValuesForxSyncSuperIt[-1],gamma))) # This is the lower integration border of 4th term.
        else: # This case has to be treated separately, because ValuesForxSyncSuperIt has the value None, which cannot be sliced.
            MinimumxgaLog=np.log10(xgammaLimitClean(x0,gamma))
        MaximumxgaLog = np.log10(DetermineBiggestIntegrationBorderOf4thTermOfEq1(ValuesForxSyncSuperIt,CurrentLastValuesForgaIt[-1]))
        FirstDivisionxgaLog = 1.000500*MinimumxgaLog
        SecondDivisionxgaLog = 0.990*MaximumxgaLog
        ValuesForxgaLog = np.concatenate((np.logspace(MinimumxgaLog,FirstDivisionxgaLog,900,endpoint=False),np.logspace(FirstDivisionxgaLog,SecondDivisionxgaLog,3000,endpoint=False),np.logspace(SecondDivisionxgaLog,MaximumxgaLog,200)))  # A sampling-range for xga.
        if not(IncludeSynchrotron):
            ValuesForpOfB18InDependenceOnxga = np.array([pOfB18(gamma,i,Usedn0,None,None) for i in ValuesForxgaLog])
            ValuesForNormalisedSpectralPPProbabilityInDependenceOnxga = np.array([NormalisedSpectralPPProbability(gamma,i,Usedn0,None,None) for i in ValuesForxgaLog])
        elif IncludeSynchrotron:
            ValuesForpOfB18InDependenceOnxga = np.array([pOfB18(gamma,i,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForxgaLog])
            ValuesForNormalisedSpectralPPProbabilityInDependenceOnxga = np.array([NormalisedSpectralPPProbability(gamma,i,Usedn0,nSyncPs,ValuesForxSyncSuperIt) for i in ValuesForxgaLog])
        pl.figure(figsize=(12, 9), num="Pair-production: Normalised spectral probability versus incident photon energy")
        pl.loglog(ValuesForxgaLog,ValuesForpOfB18InDependenceOnxga, label='Without\nescape')#, label='$p$ for $\gamma = %s \cdot x_0$' % (gamma*x0)
        pl.loglog(ValuesForxgaLog,ValuesForNormalisedSpectralPPProbabilityInDependenceOnxga, label='With\nescape', color='red')#, label='$p$ with esc. for $\gamma = %s \cdot x_0$' % (gamma*x0)
        LineNames={0: r'Helium II Lyman $\alpha$ line', 1: r'Hydrogen Lyman series', 2: r'Hydrogen Lyman $\beta$ line', 3: r'Hydrogen Lyman $\alpha$ line'}
        Linestyles=['-','--','-.',':']
        for i in range(len(x0MultiDelta)):
            x = x0MultiDelta[i]
            pl.loglog([xgammaLimitOriginal(x,gamma),xgammaLimitOriginal(x,gamma)],[10**(-40),10**(40)], linestyle=Linestyles[i%4], color='magenta')#, label='$x_{\gamma ,\mathrm{PP th}}$ for $x_0 = %s$' % i
            pl.loglog([1/x,1/x],[10**(-40),10**(40)], linestyle=Linestyles[i%4], color='olive')#, label='$1/x_i$ for %s' % LineNames[i]
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=26)
        ax.yaxis.set_tick_params(labelsize=26)
        pl.legend(loc="best", fontsize=20)
        pl.xlabel("Incident photon energy $x_{\gamma}$", fontsize=26)
        pl.ylabel("Normalised spectral PP probability     $\,\,\,\,$", fontsize=26)#"Normalised spectral PP probability $p(x_{\gamma},\gamma)$"
        pl.subplots_adjust(top=0.95, bottom=0.12, left=0.17, right=0.97)
        pl.savefig("Pair-production - Normalised probability from %s.pdf" % GetCurrentDate())

def EvaluatePOfB11AndHEPhotonLossRate(ValuesForxSyncSuperIt,ValuesFornSyncPs):
    print('\nCalling EvaluatePOfB11AndHEPhotonLossRate:')
    if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray:
        nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
    else:
        nSyncPs = None
    ValuesForxgammaLog = np.concatenate((np.logspace(np.log10(1.0/x0)-1.4,np.log10(1.0/x0)+1.8,1000,endpoint=False),np.logspace(np.log10(1.0/x0)+1.8,np.log10(1.0/x0)+3.5,100)))
    ValuesForxgammaLogPart = ValuesForxgammaLog[ValuesForxgammaLog>(1.0/x0+0.001)] # Select the values above 1/x0.
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        if ValuesForxSyncSuperIt[-1]>x0:
            ValuesForxgammaLog = np.sort(np.append(ValuesForxgammaLog,np.logspace(np.log10(1.0/ValuesForxSyncSuperIt[-1]),np.log10(1.0/x0)-0.4,50,endpoint=False)))
        ValuesForHEPhotonLossRate = [HEPhotonLossRate(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=0.01) for i in ValuesForxgammaLog]
        ValuesForPOfB11 = [POfB11(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,RelativeError=0.01) for i in ValuesForxgammaLog]
        ValuesForPOfB11WithSynchrotron = [POfB11WithSynchrotron(i,nSyncPs,ValuesForxSyncSuperIt,RelativeError=0.01) for i in ValuesForxgammaLog]
    else:
        ValuesForHEPhotonLossRate = [HEPhotonLossRate(i,Usedn0,None,None,RelativeError=0.01) for i in ValuesForxgammaLog]
        ValuesForPOfB11 = [POfB11(i,Usedn0,None,None,RelativeError=0.01) for i in ValuesForxgammaLogPart]
    ValuesForPOfB11WithoutSynchrotron = [POfB11WithoutSynchrotron(i,Usedn0) for i in ValuesForxgammaLogPart]
    #print('Sample values for xgamma:\n', ValuesForxgammaLog)
    #print('Corresponding values for POfB11:\n', ValuesForPOfB11)
    pl.rc('font', family='serif')
    Figure = pl.figure(figsize=(21, 10), num="Pair-production: Probability-rate versus photon energy")
    # On the left y-axis, the probability-rate of escape- / destruction- / absorption-events is drawn:
    LeftyAxis = Figure.add_subplot(111)
    LeftyAxis.loglog(ValuesForxgammaLog*x0,ValuesForHEPhotonLossRate, label='Spectral loss rate (sum)', linewidth=3, color='blue')#'Total probability-rate of\nescape plus destruction events'
    if IncludeSynchrotron: 
        LeftyAxis.loglog(ValuesForxgammaLogPart*x0,ValuesForPOfB11, color='red', label='Spectral absorption rate\nvia PP (soft + sync. photons)', linewidth=3)
    LeftyAxis.loglog(ValuesForxgammaLogPart*x0,ValuesForPOfB11WithoutSynchrotron, color='purple', label='Spectral absorption rate\nvia PP with contributions\nof the following lines:', linewidth=3)#'Probability-rate $P(x_{\gamma})$ of\ndestruction events due to pair-production\n(with contributions of the following lines:)'
    if Usedn0 == n0MultiDelta:
        Linestyles=['-','--','-.',':']#
        Colors=['orange','orange','orange','orange']#['#FFBF00','#FF4000','#FF00BF','#8000FF']
        LineNames={0: r'Helium II Lyman $\alpha$ line', 1: r'Hydrogen Lyman series', 2: r'Hydrogen Lyman $\beta$ line', 3: r'Hydrogen Lyman $\alpha$ line'}
        #Linewidths={0: 3.5, 1: 5.5, 2: 5, 3: 3}
        for i in range(len(x0MultiDelta)):
            x = x0MultiDelta[i]
            def POfB1OneLineOfMultiDelta(ga,xga,n0,RelativeError=integrate.quad.__defaults__[3]):   # Defining POfB1 separately for the case of n0 being one line-contribution out of n0MultiDelta. 
                return IntegrandOfB1(x,ga,xga,n0)
            def POfB11OneLineOfMultiDelta(xga,n0,RelativeError=integrate.quad.__defaults__[3]):     # Equation B11. According to a simple dimension analysis, this is in units 1/s, too. According to my analysis, it describes the probability-rate for pair-production events, in which a particle with any allowed energy ga is produced via interaction of a high-energy photon with energy xga with a soft photon background of spectral number-density n0. Thus it is a probability per unit time dt.
            # Actually the integration range is not given. It would be reasonable to integrate over the allowed range of ga, which is going from gammaminOf7 to gammamaxOf7.
                if xga<1.0/x0:
                    return 0.0
                else:
                    LowestIntegrationBorder=gammaminOf7(x0,xga)
                    BiggestIntegrationBorder=gammamaxOf7(x0,xga)
                    NumberOfBorders = (max(5,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder))) if xga<=10.0**13.0 else 10.0**4.0) # The else-case is needed since severe integration errors occur for xga>10**13.
                    ListOfIntegrationBorders = np.linspace(LowestIntegrationBorder,BiggestIntegrationBorder,NumberOfBorders)
                    POfB11Result=0.0
                    for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                        POfB11Result += integrate.quad(POfB1OneLineOfMultiDelta, LeftBorder, RightBorder, args=(xga,n0), epsrel=RelativeError)[0]
                    return POfB11Result
            ValuesForxgammaLogLine = np.concatenate((np.logspace(np.log10(1.0/x)+0.001,np.log10(1.0/x)+1.8,500,endpoint=False),np.logspace(np.log10(1.0/x)+1.8,np.log10(1.0/x)+3.5,100)))
            ValuesForPOfB11OneLineOfMultiDelta = [POfB11OneLineOfMultiDelta(i,Usedn0,RelativeError=0.01) for i in ValuesForxgammaLogLine]    # Evaluation for n0=Usedn0 and xga=ValuesForxgammaLog but only one value of x (one single line) is considered.
            LeftyAxis.loglog(ValuesForxgammaLogLine*x0,ValuesForPOfB11OneLineOfMultiDelta, linestyle=Linestyles[i%4], linewidth=2, color=Colors[(i//3)%4], label=LineNames[i])#'Line at $x_{0,i}=%g$'
            pl.loglog([x0*1/x,x0*1/x],[10**(-20),10**(20)], linestyle=Linestyles[i%4], color='olive')#, label='$1/x_i$ for %s' % LineNames[i]
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        LeftyAxis.loglog(ValuesForxgammaLog*x0,ValuesForPOfB11WithSynchrotron, color='magenta', label='Spectral absorption rate\nvia PP (sync. photons)', linewidth=3)
    if IncludeHEPhotonEscape:
        ValuesForHEPhotonEscapeRate = [HEPhotonEscapeRate(i) for i in ValuesForxgammaLog] # The escape-probability-rate.
        LeftyAxis.loglog(ValuesForxgammaLog*x0,ValuesForHEPhotonEscapeRate, label='Spectral escape rate', linewidth=3, color='green')#'Probability-rate of escape events'
    pl.legend(loc="upper left", fontsize=18)
    LeftyAxis.xaxis.set_tick_params(labelsize=35, direction='inout', width=5, length=16, pad=-1)
    LeftyAxis.yaxis.set_tick_params(labelsize=35, direction='inout', width=5, length=16, right=False)
    pl.xlabel("$x_\gamma x_{1}$", fontsize=42)
    pl.ylabel("Probability rate in s$^{-1}$", fontsize=41)
    Ticks=LeftyAxis.get_yticks() # Store the locations of the ticks.
    print(Ticks)
    LeftyAxis.set_ylim(2*10**(-5),2*10**(-2)) # This setting of the limits might change from plot to plot. Usage for pion-decay: Ticks[5],Ticks[8]
    LeftyAxis.set_xlim(ValuesForxgammaLog[0]*x0,ValuesForxgammaLog[-1]*x0)
    # The right y-axis shows the optical depth, corresponding to the probability-rate, and is determined according to OpticalDepthOfHEPhotons.
    RightyAxis = Figure.add_subplot(111, sharex=LeftyAxis, frameon=False)
    RightyAxis.loglog()
    RightyAxis.yaxis.tick_right()
    RightyAxis.yaxis.set_label_position("right")
    RightyAxisTicks=np.asarray([10**(-4),10**(-3),10**(-2),10**(-1),10**(0),10**(1),10**2])
    RightyAxis.yaxis.set_ticks(RightyAxisTicks/MeanEscapeTime) # The values of the spectral loss-rate, that correspond to the ticked optical depth. This essentially is loss-probability-rate-of-process = tau / MeanEscapeTime, where tau represents RightyAxisTicks.
    RightyAxis.yaxis.set_ticklabels(RightyAxisTicks)
    RightyAxis.set_ylim(2*10**(-5),2*10**(-2)) # Usage for pion-decay: Ticks[5],Ticks[8]
    pl.ylabel("Corresponding optical depth", fontsize=41)
    RightyAxis.xaxis.set_tick_params(which='both', bottom=False, top=False, labelbottom=False)
    RightyAxis.yaxis.set_tick_params(labelsize=35, direction='inout', width=5, length=16, pad=9)
    RightyAxis.set_xlim(ValuesForxgammaLog[0]*x0,ValuesForxgammaLog[-1]*x0)
    LeftyAxis.spines['left'].set_linewidth(4)
    LeftyAxis.spines['right'].set_linewidth(4)
    LeftyAxis.spines['top'].set_linewidth(4)
    LeftyAxis.spines['bottom'].set_linewidth(4)
    # The following is to create a different top x-axis:
    UpperAxis=LeftyAxis.twiny()
    UpperAxis.set_xscale('log')
    UpperAxis.set_xlim(LeftyAxis.get_xlim())
    UpperAxisTicks=np.asarray([10**(-3),10**(-2),10**(-1),10**(0),10**(1)]) # The tick-labels in TeV.
    UpperAxis.set_xticks(UpperAxisTicks*x0*e*10**12/(me*c**2)) # Conversion from values of the top axis to the positions on the bottom axis.
    UpperAxis.set_xticklabels(UpperAxisTicks)
    UpperAxis.xaxis.set_tick_params(labelsize=35, direction='inout', width=5, length=16, pad=8)
    UpperAxis.set_xlabel("$x_\gamma m_{\mathrm{e}} c^2 / \, \mathrm{TeV}$", fontsize=42, labelpad=20)
    pl.subplots_adjust(top=0.85, bottom=0.16, left=0.10, right=0.89)
    pl.savefig("Run %s, Photon-loss-rate - Probability-rate versus photon energy from %s.pdf" % (RunIdentifier,GetCurrentDate()))
    # The following conclusions can be drawn: For n0=n0Planck evaluation seems reasonable. The value of POfB11 decreases for increasing xgamma. For n0=n0Delta evaluation is fast and without error messages. Output seems reasonable. The value of POfB11 decreases for increasing xgamma. For n0=n0PL evaluation seems reasonable. The value of POfB11 slowly decreases for increasing xgamma. For n0=n0Exp evaluation seems reasonable. The value of POfB11 slowly decreases for increasing xgamma.


## Considering cut-offs

# For the evaluation of the iteration and super-iteration and for the computation of the integrals in part 5 and 6, it is mandatory to define upper (and lower) cut-offs of the various distributions, expressions, quantities and terms. This is done now:

# The following two cut-offs were obtained by considering the course of SyncPsSpectralNumberDensity and the reasonable integration borders of SyncPsNumberDensity (More precisely: These two coefficients were found by computing SyncPsNumberDensity for different integration borders (cf. FinalCalculationResultQueue). These numbers are those integration borders beyond which SyncPsNumberDensity does not increase any more. Hence, the relative error is kept reasonably small.)
def DeterminexSync0(NElectronsga0):
    '''The physically motivated upper cut-off of the spectral number-density of SyncPs.'''
    return 5.0*SynchrotronCriticalxSync(NElectronsga0)

def DeterminexSync1(Lowestga):
    '''The lower cut-off of the spectral number-density of SyncPs. This might, in any occurring case, be smaller than 0.1.'''
    return 10**(-6)*SynchrotronCriticalxSync(Lowestga)

def DeterminexSync0Used(xSync0):
    '''The used upper cut-off of the spectral number-density of SyncPs, respecting the restriction x<<1.'''
    if xSync0 <= 0.1: # The case satisfying Zdz's assumption x<<1.
        return xSync0
    else: # The case where the SyncP-distribution reaches energies too high with respect to the assumption.
        if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
            print("\nAttention: The SyncP-distribution's upper cut-off had to be reduced to satisfy the assumption x<<1.")
        return 0.1 # Actually, this upper border should be about 0.1. However, until v2_3, it was set to x0 for simplicity.

if IncludeSynchrotron and PerformanceMode!='CommonIteration':
    def DeterminexTotal0(xSync0Used):
        '''The overall upper cut-off of all the low-energy photons.'''
        return max(x0,xSync0Used)
else:
    def DeterminexTotal0(xSync0Used):
        '''The overall upper cut-off of all the low-energy photons.'''
        return x0

# Consider IntegralOfBracketsOfEq1, which is equal to Zdz's equation 5, and define its upper cut-off.
def DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,NElectronsga0):
    '''I could show that IntegralOfBracketsOfEq1 is always vanishing beyond this value of xga.
    Generally, ValuesForxSyncSuperIt can be an array or None. If it is an array, only its last value is valuable because it is the current xSync0Used.'''
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        return max(xmaxOfA20(NElectronsga0,x0),xmaxOfA20(NElectronsga0,ValuesForxSyncSuperIt[-1])) # The first argument was the cut-off without synchrotron. 
    else:
        return xmaxOfA20(NElectronsga0,x0)

# I could show that the integrand of the fourth term has a maximal upper cut-off at a certain value BiggestIntegrationBorderOf4thTermOfEq1. This value is also the upper cut-off of nHEPhotonsSumSpectralInxga, defined in part7. (Obsolete comment: This value could be an input value. The upper integration border of the following integral should be infinity according to Zdz. But firstly, Dotni(xga) is cutted-off at HEPhotonsxga0 and secondly, pOfB18 in dependence on xga is decreasing, so it is cutted-off via the decrease of IntegralOfBracketsOfEq1 (which is an effect of NElectrons decreasing for big enough ga). This is visible in EvaluateIntegrandOfFirstSummandOf4thTerm(). This might be a critical parameter.) (Comment from version 5: I could show that for StartingNElectrons=NElectronsRecommended, the iterated integrand of the fourth term always has an upper cut-off at Integrandxga0. Hence, this can be taken as biggest integration-border.)
if UsedDotni==DotniDelta or UsedDotni==DotniZero: # Note that in this case, the upper cut-off of only the 1st summand of the integrand of the 4th term is meant, hence it is equal to the upper cut-off of equation 5.
    def DetermineBiggestIntegrationBorderOf4thTermOfEq1(ValuesForxSyncSuperIt,NElectronsga0):
        return DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,NElectronsga0)
else: # In this case, the upper cut-off of the entire integrand of the 4th term is meant.
    def DetermineBiggestIntegrationBorderOf4thTermOfEq1(ValuesForxSyncSuperIt,NElectronsga0):
        return max(HEPhotonsxga0,DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,NElectronsga0))

# I could show (cf. Determining an upper cut-off of the fourth term.png) that there is a maximal upper cut-off UpperCutOffOf4thTermOfEq1 in ga, above which the fourth term of equation1 is vanishing for sure. This value is generally determined via DetermineUpperCutOffOf4thTermOfEq1 from input values as well as from xSync0Used and NElectronsga0.
# Furthermore, I could show (cf. Determining an upper cut-off of NElectrons.png), that the electron distribution NElectrons has a maximal upper cut-off at a certain value called NElectronsga0. This value is determined via the function DetermineNElectronsga0, which is generally dependent on input values and on xSync0Used.
# In cases of vanishing UsedDotNi or UsedDotni the values of these cut-offs can be further restricted and thus have to be defined separately, however in all cases, the arguments are kept the same.
if UsedDotNi==DotNiZero:
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        def DetermineUpperCutOffOf4thTermOfEq1(xSync0Used,NElectronsga0=111):
            return max(gammamaxOf7(x0,HEPhotonsxga0),gammamaxOf7(xSync0Used,HEPhotonsxga0)) # The first argument of max(,) is the conventional cut-off. Actually, this is the same function-body as in DetermineNElectronsga0. The keyword argument is only there for syntactical reasons.
        def DetermineNElectronsga0(xSync0Used):
            return DetermineUpperCutOffOf4thTermOfEq1(xSync0Used)
    else:
        def DetermineUpperCutOffOf4thTermOfEq1(xSync0Used,NElectronsga0=111):
            return gammamaxOf7(x0,HEPhotonsxga0) # The conventional cut-off. Actually, this is the same function-body as in DetermineNElectronsga0. The keyword argument is only there for syntactical reasons.
        def DetermineNElectronsga0(xSync0Used):
            return DetermineUpperCutOffOf4thTermOfEq1(xSync0Used) # The conventional cut-off.
    
elif UsedDotni==DotniZero:
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        def DetermineUpperCutOffOf4thTermOfEq1(xSync0Used,NElectronsga0):
            return max(gammamaxOf7(x0,DetermineUpperCutOffOfEq5(xSync0Used,Injectedga0)),gammamaxOf7(xSync0Used,DetermineUpperCutOffOfEq5(xSync0Used,Injectedga0))) # The first argument is the conventional cut-off.
    else:
        def DetermineUpperCutOffOf4thTermOfEq1(xSync0Used,NElectronsga0):
            return gammamaxOf7(x0,DetermineUpperCutOffOfEq5(xSync0Used,Injectedga0)) # The conventional cut-off.   
    def DetermineNElectronsga0(xSync0Used):
        return Injectedga0 # The conventional cut-off.
    
else:
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        def DetermineUpperCutOffOf4thTermOfEq1(xSync0Used,NElectronsga0):
            return max(gammamaxOf7(x0,max(HEPhotonsxga0,DetermineUpperCutOffOfEq5(xSync0Used,NElectronsga0))),gammamaxOf7(xSync0Used,max(HEPhotonsxga0,DetermineUpperCutOffOfEq5(xSync0Used,NElectronsga0)))) # The first argument of max(,) is the conventional cut-off.
        def DetermineNElectronsga0(xSync0Used):
            return max(max(Injectedga0,DetermineUpperCutOffOf4thTermOfEq1(xSync0Used,Injectedga0)),max(Injectedga0,DetermineUpperCutOffOf4thTermOfEq1(xSync0Used,Injectedga0))) # The first argument of max(,) is the conventional cut-off.
    else:
        def DetermineUpperCutOffOf4thTermOfEq1(xSync0Used,NElectronsga0):
            return gammamaxOf7(x0,max(HEPhotonsxga0,DetermineUpperCutOffOfEq5(xSync0Used,NElectronsga0))) # The conventional cut-off.       
        def DetermineNElectronsga0(xSync0Used):
            return max(Injectedga0,DetermineUpperCutOffOf4thTermOfEq1(xSync0Used,Injectedga0)) # The conventional cut-off.

ThomsonKNBoundaryExPs = 1.0/(8.0*x0) # This is the approximate border between the Thomson-regime and the Klein-Nishina-regime with respect to scattering on external photons. Comment of version v1_1: Now, I don't see any more why the third term should vanish at and below this value!

LowestgaParameter1 = 0.01 # Physically, the distribution of the electrons won't have a lower cut-off. Even though, one cannot extend the evaluation range down to 0. This parameter, together with the following function, determines the chosen lower cut-off of the sampling range in units of x0 or xSync0Used, respectively. INPUT VALUE!
LowestgaParameter2 = 30 # This is the minimally possible value that Lowestga can attain. INPUT VALUE!
def DetermineLowestga(xSync0Used):
    '''The chosen lower border of the evaluation- and plotting-range of ga.
    The max(LowestgaParameter2,...) ensures that Lowestga is minimally equal to LowestgaParameter2. This implicates that for xSync0Used=0.1 and for interactions with SyncPs, the Thomson-regime is never reached.'''
    if IncludeSynchrotron and PerformanceMode!='CommonIteration':
        return max(LowestgaParameter2,min(LowestgaParameter1/x0,LowestgaParameter1/xSync0Used))
    else:
        return max(LowestgaParameter2,LowestgaParameter1/x0) # The conventional cut-off.

CurrentNElectronsga0 = DetermineNElectronsga0(x0) # Actually, one had to insert CurrentxSync0Used here as argument. However this would cause a bootstrap problem, because CurrentxSync0Used is in the end depending on CurrentNElectronsga0. As the first iteration of the super-iteration will be performed without synchrotron-photons as targets, it is valid to use x0 here instead of CurrentxSync0Used.
CurrentLowestga = DetermineLowestga(x0) # Again, one can use x0 here instead of CurrentxSync0Used.
CurrentxSync0 = DeterminexSync0(CurrentNElectronsga0)    
CurrentxSync1 = DeterminexSync1(CurrentLowestga)
CurrentxSync0Used = DeterminexSync0Used(CurrentxSync0)
CurrentxTotal0 = DeterminexTotal0(CurrentxSync0Used)


## Initialisation of the iteration

# Preparation of the iteration (sampling of points of gamma):

if UsedIterationScheme == 'PointsFromRToLIterationPointwise': # Create the array of sampling points. In the case 'IterationStepByStepPointsFromLToR', this preparation is performed in the iteration-loop EvaluateIteration because a changing array of sampling points is possible there:

    Range5Parameter = 40.0 # This affects the compilation of the sampling range. INPUT VALUE!
    Range6Parameter = 4.0 # This affects the compilation of the sampling range. INPUT VALUE!
    Range7Parameter = 22.0 # This affects the compilation of the sampling range. INPUT VALUE!
    NumberOfSamplingPointsForRange1 = 32 # Number of sampling points for the lowest energetic range below the Thomson-KN-boundary. INPUT VALUE!
    NumberOfSamplingPointsForRange2 = 6 # Number of sampling points for the range 2, immediately above the Thomson-KN-boundary. INPUT VALUE!
    NumberOfSamplingPointsForRange3 = 6 # Number of sampling points for the range 3. INPUT VALUE!
    NumberOfSamplingPointsForRange4 = 30 # Number of sampling points for the range 4. INPUT VALUE!
    NumberOfSamplingPointsForRange5 = 30 # Number of sampling points for the range 5. INPUT VALUE!
    NumberOfSamplingPointsForRange6 = 25 # Number of sampling points for the range 6. INPUT VALUE!
    NumberOfSamplingPointsForRange7 = 30 # Number of sampling points for the range 7, the highest-energetic range. INPUT VALUE!
    
    def CompileSamplingRange(Lowestga,NElectronsga0,xSync0Used):
        Range1Beginning = Lowestga # The beginning of the sampling range.
        ThomsonKNBoundarySyncPs = 1.0/(8.0*xSync0Used) # This is the approximate border between the Thomson-regime and the Klein-Nishina-regime with respect to scattering on synchrotron-photons
        ThomsonKNBoundaryLower = min(ThomsonKNBoundaryExPs,ThomsonKNBoundarySyncPs) # The lower one of the two boundaries.
        ThomsonKNBoundaryHigher = max(ThomsonKNBoundaryExPs,ThomsonKNBoundarySyncPs) # The higher one.
        Range2Beginning = ThomsonKNBoundaryExPs
        Range3Beginning = 2.0*ThomsonKNBoundaryExPs # 1/(4*x0), which is the set-in of pair-production on ExPs. In other words, this is the value of ga, below and at which the contribution due to external photons to the fourth term of the kinetic equation is vanishing. 1/(4*xSync0Used) would be the value of ga, below and at which the contribution due to SyncPs to the fourth term of the kinetic equation is vanishing.
        Range4Beginning = 4.0*ThomsonKNBoundaryExPs # At 1/(2*x0) the contribution due to ExPs to FourthTermOfEq1 has a kink.
        Range5Beginning = Range5Parameter*ThomsonKNBoundaryExPs
        Range6BeginningLog = (1.0*np.log10(ThomsonKNBoundaryExPs)+Range6Parameter*np.log10(NElectronsga0))/(Range6Parameter+1) # This is a point, which is situated between ThomsonKNBoundary and NElectronsga0 in the ratio Range6Parameter to 1 in logarithmic space.
        Range7BeginningLog = (1.0*np.log10(ThomsonKNBoundaryExPs)+Range7Parameter*np.log10(NElectronsga0))/(Range7Parameter+1) # This is a point, which is situated between ThomsonKNBoundary and NElectronsga0 in the ratio Range7Parameter to 1 in logarithmic space.

        ValuesForgaLogRange1 = np.logspace(np.log10(Range1Beginning),np.log10(Range2Beginning),NumberOfSamplingPointsForRange1,endpoint=False) # The sampling-range for ga, at very low values of ga.
        ValuesForgaLogRange2 = np.logspace(np.log10(Range2Beginning),np.log10(Range3Beginning),NumberOfSamplingPointsForRange2,endpoint=False) # Sampling-range 2.
        ValuesForgaLogRange3 = np.logspace(np.log10(Range3Beginning),np.log10(Range4Beginning),NumberOfSamplingPointsForRange3,endpoint=False) # Sampling-range 3.
        ValuesForgaLogRange4 = np.logspace(np.log10(Range4Beginning),np.log10(Range5Beginning),NumberOfSamplingPointsForRange4,endpoint=False) # Sampling-range 4.
        ValuesForgaLogRange5 = np.logspace(np.log10(Range5Beginning),Range6BeginningLog,NumberOfSamplingPointsForRange5,endpoint=False) # Range 5.
        ValuesForgaLogRange6 = np.logspace(Range6BeginningLog,Range7BeginningLog,NumberOfSamplingPointsForRange6,endpoint=False) # The sampling-range 6.
        ValuesForgaLogRange7 = np.logspace(Range7BeginningLog,np.log10(NElectronsga0),NumberOfSamplingPointsForRange7,endpoint=False) # The highest-energetic sampling-range.
    
        ValuesForgaLogTotalRange=np.append(np.append(np.append(np.append(np.append(np.append(np.append(np.asarray([NElectronsga0]),ValuesForgaLogRange1),ValuesForgaLogRange2),ValuesForgaLogRange3),ValuesForgaLogRange4),ValuesForgaLogRange5),ValuesForgaLogRange6),ValuesForgaLogRange7) # Initialise the complete list of sampling points.
        ValuesForgaLogTotalRange = np.sort(ValuesForgaLogTotalRange)
        return ValuesForgaLogTotalRange

    CurrentValuesForgaIt = CompileSamplingRange(CurrentLowestga,CurrentNElectronsga0,x0) # Initialise the current sampling range, that is used during the iteration and perhaps changed during the super-iteration. Again x0 is used instead of xSync0Used, because the later is unknown until now.

## Electron distribution

# The spectral number-density N(ga) of electrons is now considered. It is given in 1/m^3. This quantity is to be determined, however an initial function has to be specified.

def NElectronsZero(ga):
    '''Provisionally it can be set equal to an easy quantity. Attention, this is not physically reasonable.'''
    return 0.0

def NElectronsUnity(ga):
    '''Provisionally it can be set equal to an easy quantity. Attention, this is not physically reasonable.'''
    return 1.0

def NElectronsRecommendedZdzExact(ga):
    '''The distribution, that is recommended by Zdziarski to be used as starting distribution for the iteration. In this case, the exact recommended value is used.'''
    return UsedDotNi(ga)/COfA21WithoutSynchrotron(ga,Usedn0)

AveragedValueOfCOfA21WithoutSynchrotron = (COfA21WithoutSynchrotron(Injectedga1,Usedn0)+COfA21WithoutSynchrotron(Injectedga0,Usedn0))/2.0
def NElectronsRecommendedZdzApprox(ga):
    '''However the repeated evaluation of COfA21WithoutSynchrotron can take much time. By looking at EvaluateLowerBorderOfNonVanRangeOfCOfA1AndCOfA21 one can realise, that COfA21WithoutSynchrotron does not vary strongly over the range of ga, in which UsedDotNi is non-vanishing. Hence one can use one single averaged value of COfA21WithoutSynchrotron for all ga.'''
    return UsedDotNi(ga)/AveragedValueOfCOfA21WithoutSynchrotron
# Plotting-test:    X=np.linspace(Injectedga1,Injectedga0,301);Y=np.array([NElectronsRecommendedZdzApprox(i) for i in X]);pl.plot(X*x0,X*Y)

def NElectronsRecommendedAllProcesses(ga):
    '''All included processes (IC-scattering on the soft background photons + escape + synchrotron-radiation) are considered in the denominator.'''
    return UsedDotNi(ga)/ElectronLossRate(ga,Usedn0,lambda xSync : 0.0,[0.001, 1000],None,0.01) # The arguments-choice nSyncPs= lambda xSync : 0.0 (a lambda function should be used here, because the argument has to be a function) and ValuesForxSyncSuperIt=[0.001, 1000] guarantees that the IC-scattering term on the synchrotron-background is set =0.
# Plotting-test:   X=np.logspace(np.log10(Injectedga1),np.log10(Injectedga0),301);Y=np.array([NElectronsRecommendedAllProcesses(i) for i in X]);pl.loglog(X*x0,X*Y)

def NElectronsPremonition(ga):
    return 7.0*10**(5)*ga**(-1.0)#0.0015*ga**(-1.4)*Heavi(1.0,ga,1010.0) # This is a smart guess, which comes from previous findings.
# Plotting-test:    X=np.linspace(0.005/(x0),Injectedga0,301);Y=np.array([NElectronsPremonition(i) for i in X]);pl.loglog(X*x0,Y*X);pl.ylim(10.0**(-3),10.0**(3.0))

if InjectionType=='ResultsOfOldIt':
    StartingNElectrons = OldCurrentNElectrons
    StartingNElectrons.__name__ = 'OldCurrentNElectrons'
elif InjectionType=='Chosen':
    StartingNElectrons = NElectronsRecommendedZdzApprox # The electron distribution, that is used at the beginning of the iteration. INPUT VALUE!
CurrentNElectrons = StartingNElectrons        # The actually used electron distribution. It will be changed later during the iteration.

ArtificialPushUpOfInitialisation = 2.0 # This is used in determining the initialisation of NElectrons in the case of super-iterations. To attain the initialisation of NElectrons in one SuperIt-step, CurrentNElectrons (in other words the result of the last SuperIt-step) is multiplied with ArtificialPushUpOfInitialisation. If ArtificialPushUpOfInitialisation is equal to 1.0, CurrentNElectrons is directly used as initialisation. The need for this is to impede the sticking of NElectrons at the initialisation. INPUT VALUE!

DefaultTake1stIntRangeOf4thTermIntoAccount = True # If True, the 1st integration range is taken into account in every iteration step. If False, the 1st range is taken into account only after a certain, required precision was achieved during an iteration (cf. MoreIterationsNecessary).

if UsedIterationScheme == 'PointsFromRToLIterationPointwise': 
    def CreateANewNElectronsIterated(CurrentNElectrons,ValuesForgaIt,SuperItCounter):
        '''This creates, initialises and returns a new dictionary to contain the iterated values of NElectrons.'''
        NElectronsIterated = {} # A data-structure to save the results of the iteration. In the case 'PointsFromRToLIterationPointwise', it will be a dictionary of dictionaries. The outer dictionary has the elements of ValuesForgaIt as keys. The value of each key is a dictionary, which will contain the keys 'StackOfIteratedValues', 'IterationFinished', 'OneStepDoneAfterTaking1stIntRangeOf4thTermIntoAccount', 'Take1stIntRangeOf4thTermIntoAccount', 'RelativeChangeInFinalIterationStep' and 'IterationCounter'. 'StackOfIteratedValues' contains an array with the iterated values of NElectrons that correspond to the respective key.
        # Initialise NElectronsIterated with the starting function of NElectrons:
        for ga in ValuesForgaIt:
            NElectronsIterated[ga] = {}
            if SuperItCounter==1:
                NElectronsIterated[ga]['StackOfIteratedValues'] = np.asarray([0.0,CurrentNElectrons(ga)]) # The initial filling with two elements is necessary, for the convergence checking function MoreIterationsNecessary being applicable.
            elif SuperItCounter>1:
                NElectronsIterated[ga]['StackOfIteratedValues'] = np.asarray([0.0,ArtificialPushUpOfInitialisation*CurrentNElectrons(ga)]) # In these cases, CurrentNElectrons is not directly used as initialisation.
            NElectronsIterated[ga]['IterationFinished'] = False # This flags, whether the iteration at this ga was already finished (True) or whether it is still running or still has to be started. For this flag to be set True, all criteria in the testing function MoreIterationsNecessary have to be satisfied.
            NElectronsIterated[ga]['OneStepDoneAfterTaking1stIntRangeOf4thTermIntoAccount'] = False # This flags, whether one additional iteration step at this ga was already performed after the inclusion of the 1st integration range of the 4th term.
            NElectronsIterated[ga]['Take1stIntRangeOf4thTermIntoAccount'] = DefaultTake1stIntRangeOf4thTermIntoAccount # This flags, whether the 1st integration range is supposed to be taken into account in the computation of the 4th term during the beginning of an iteration.
            NElectronsIterated[ga]['RelativeChangeInFinalIterationStep'] = None # This will, after the iteration was finished at this point of ga, store the relative change that occurred from the second last to the first last iterated value, in other words, it stores the actually achieved accuracy.
            NElectronsIterated[ga]['IterationCounter'] = None # This will, after the iteration was finished at this point of ga, store the value of IterationCounter, in other words, it stores the number of iteration steps, that have been necessary.
        return NElectronsIterated

    if PerformanceMode!='SuperIteration':
        NElectronsIterated = CreateANewNElectronsIterated(StartingNElectrons,CurrentValuesForgaIt,1) # NElectronsIterated is a global object in the case of a common iteration without synchrotron-back-reaction. 

    CurrentValuesForNElectrons = np.asarray([StartingNElectrons(ga) for ga in CurrentValuesForgaIt]) # Recreate the list of latest values of NElectrons.
    CurrentLastValuesForgaIt = CurrentValuesForgaIt # This is necessary for some functions to be callable also in this 'PointsFromRToLIterationPointwise' case.

elif UsedIterationScheme == 'IterationStepByStepPointsFromLToR':
    NElectronsIterated = {} # A data-structure to save the findings of the iteration. In the case 'IterationStepByStepPointsFromLToR', it will be a dictionary of dictionaries. The outer dictionary has keys 'i. Iteration' where i is the number of the respective iteration. Each value of the outer dictionary is again a dictionary. These inner dictionaries have the keys 'Values for gamma', 'Interpolated object' and 'Values for next NElectrons' with the corresponding values ValuesForgaIt, InterpolatedObjectForValuesForNextNElectrons and ValuesForNextNElectrons.
    CurrentIterationCounter=0 # This has a meaning and is used only in the case UsedIterationScheme = IterationStepByStepPointsFromLToR. It is an integer, that gives the number of iteration steps. If it is equal to zero, the first step was not performed yet or is just performing.
    CurrentValuesForgaIt = np.logspace(np.log10(10**5),np.log10(10**7),50) # The existence of ValuesForgaIt is necessary for the following functions to be callable even though the iteration was not performed yet.
    CurrentLastValuesForgaIt = CurrentValuesForgaIt # The existence of LastValuesForgaIt is necessary for the following functions to be callable even though the iteration was not performed yet.

def CompareNElectronsSuperIt(NElectronsSuperIt1,NElectronsSuperIt2):
    '''If no output is printed, the dictionaries are equal.'''
    def CompareNElectronsSuperItAuxiliary(NESuperIt1,NESuperIt2):
        for SuperItCounter in NESuperIt1.keys():
            for gamma in NESuperIt1[SuperItCounter].keys():
                if NESuperIt1[SuperItCounter][gamma]['IterationCounter']!=NESuperIt2[SuperItCounter][gamma]['IterationCounter']:
                    print('False')
                if len(NESuperIt1[SuperItCounter][gamma]['StackOfIteratedValues'])!=len(NESuperIt2[SuperItCounter][gamma]['StackOfIteratedValues']):
                    print('False')
                for k in range(len(NESuperIt1[SuperItCounter][gamma]['StackOfIteratedValues'])):
                    if NESuperIt1[SuperItCounter][gamma]['StackOfIteratedValues'][k]!=NESuperIt2[SuperItCounter][gamma]['StackOfIteratedValues'][k]:
                        print('False')
                if NESuperIt1[SuperItCounter][gamma]['OneStepDoneAfterTaking1stIntRangeOf4thTermIntoAccount']!=NESuperIt2[SuperItCounter][gamma]['OneStepDoneAfterTaking1stIntRangeOf4thTermIntoAccount']:
                    print('False')
                if NESuperIt1[SuperItCounter][gamma]['RelativeChangeInFinalIterationStep']!=NESuperIt2[SuperItCounter][gamma]['RelativeChangeInFinalIterationStep']:
                    print('False')
                if NESuperIt1[SuperItCounter][gamma]['IterationFinished']!=NESuperIt2[SuperItCounter][gamma]['IterationFinished']:
                    print('False')
                if NESuperIt1[SuperItCounter][gamma]['Take1stIntRangeOf4thTermIntoAccount']!=NESuperIt2[SuperItCounter][gamma]['Take1stIntRangeOf4thTermIntoAccount']:
                    print('False')
    CompareNElectronsSuperItAuxiliary(NElectronsSuperIt1,NElectronsSuperIt2) # All items of NElectronsSuperIt1 are checked with NElectronsSuperIt2.
    CompareNElectronsSuperItAuxiliary(NElectronsSuperIt2,NElectronsSuperIt1) # All items of NElectronsSuperIt2 are checked with NElectronsSuperIt1.

# Determine the power-density of the injected electrons, which is in units of 1/(m^3*s):
def InjectedPowerDensitySpectralInga(ga,DotNi):
    return ga*DotNi(ga)

def InjectedPowerDensity(DotNi,RelativeError=integrate.quad.__defaults__[3]):
    if isinstance(DotNi, types.FunctionType):
        LowestIntegrationBorder = Injectedga1
        BiggestIntegrationBorder = Injectedga0
        if LowestIntegrationBorder >= BiggestIntegrationBorder:
            return 0.0
        else:
            NumberOfBorders = max(6,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
            ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
            InjectedPowerDensityResult = 0.0
            for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                InjectedPowerDensityResult += integrate.quad(InjectedPowerDensitySpectralInga, LeftBorder, RightBorder, args=(DotNi), epsrel=RelativeError)[0]
            return InjectedPowerDensityResult
    elif isinstance(DotNi, interp1d):
        ListOfIntegrationBorders = OldValuesForgaIt
        InjectedPowerDensityResult = 0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            InjectedPowerDensityResult += integrate.quad(InjectedPowerDensitySpectralInga, LeftBorder, RightBorder, args=(DotNi), epsrel=RelativeError)[0]
        return InjectedPowerDensityResult

InjectedElectronPowerDensity = InjectedPowerDensity(UsedDotNi,RelativeError=0.01) # In units of 1/(m^3*s).

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n\n--------------------------- Consider the injected electrons -----------------------------\n')
    print('Injected electron power-density:', InjectedElectronPowerDensity, "1/(m^3*s)")
    print('    Corresponding injected electron power:', SphericalVolume('InteractionRadius')[1]*InjectedElectronPowerDensity*me*c**2, "W (using R=%s)" % SphericalVolume('InteractionRadius')[0])
    print('    Corresponding injected electron power:', SphericalVolume('RClouds')[1]*InjectedElectronPowerDensity*me*c**2, "W (using R=%s)" % SphericalVolume('RClouds')[0])


# Documentation of version 8: Pair-production section was developed. Minor changes in previous sections were performed. Concerning the IC-section, I did realise certain connections concerning the minimally possible value of gaP, LowerBorderOfNonVanRangeOfCOfA1, the integration range of A21 and the so far unsolved problem with the Dirac-Delta-function. This will be treated in version 9.

# Documentation of version 10: Beside minor changes, the pair-production section was completed up to here. 

# Documentation of version 11: Structuring of loose stetements into _functions was done and test-blocks have been included. Minor changes were performed.

# Documentation of version 12: GammaLimit, GammaMinOnxga and xmaxOfA20 and their evaluations have been included.

# Documentation of version 13: The following was added: _EvaluategammaminmaxOf7(), xgammaLimit(x0,ga), _EvaluatexgammaLimit(), _EvaluatepOfB18InDependenceOnxga()

# Documentation of version 14: The order of structuring was reorganised. Now the definitions always appear at first. Afterwards the _Evaluate-definitions and test blocks appear. Furthermore, the names of the labels were changed. The definition of the Heaviside-function was included. GetCurrentDate() was added. Most test-blocks have been deleted.

# Documentation of version 15: The subsections "Injection of primary particles" and "Injection of primary HE-photons" were moved to part 1 from part 2. The normalisation of DeltaCoefficient was introduced. The natural constants were precised. Auxiliary functions were interchanged with part 2 and 3. Minor changes.

# Documentation of version 16: COfA1, COfA21 and POfB1 have been altered. For all functions, beginning with underscore, the underscore was deleted and the function-name was prefixed with "Evaluate", if not done yet. Minor changes.

# Documentation of version 17: xgammaLimit and xgammaLimitOriginal were distinguished. DotNiRelMax was introduced.

# Documentation of version 18: IntelligentStorage and ApplyIntelligentStorage have been added. Instead, the complicated and intricate procedure within SecondTermOfNumeratorOfEq8 to use AuxiliaryQuantityOfSecondTermOfNumeratorOfEq8 could be erased. The if-test in EvaluateExportInputValueshas been deleted as it is not needed any more.

# Documentation of version 19: Up to now, there occurred differences in the results between python2- and python3-evaluations. These were due to python2 treating 1/2 as integer-division. To find remedy / make the code also valid for python2, all numbers that are definitively no integer are rendered with a decimal point. Furthermore, quantities that can be simplified, like e. g. 1/2, are simplified, e. g. to 0.5.
# The Heaviside-function was implemented another way. Now, at the steps, it adopts the value 1.0 instead of 0.5.
# DeleteFile was added.
# NElectronsGaussian was removed.
# mToMpc and MpcTom was added.
# The conversion from HEPhotonsGaussianCoefficient to HEPhotonsGaussianNorm and vice versa and its manual input was added.
# The integration of COfA1 was divided into a number of integrations along part-intervals.
# The epsrel-keyword-argument has been introduced in some integrals.
# Minor changes.

# Documentation from v1_0 on will proceed in part 3...
