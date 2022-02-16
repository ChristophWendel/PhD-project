## ------------------------------------------------------------------------------------     o   o
## Electromagnetic cascades                                                                   I
##                                                                                          \.../
## Numerical solution of equation 8; HE-photons' spectrum
## ------------------------------------------------------------------------------------

if __name__ == "__main__": # This prevents evaluation in child-processes.
    print('\n\n_________________________________________________________________________________________')
    print('_________________________________________________________________________________________\n')
    print('-----------------------------    AGN- and cascade-code    -------------------------------')
    print('_________________________________________________________________________________________')
    print('_________________________________________________________________________________________\n\n')

import os
import multiprocessing

def SearchAFile0rDirectory(Name,WhereToSearch):
    '''Search for files or directories called Name in the directory WhereToSearch and in all its subfolders. An example for Name would be "Lost File.txt" or "Any Dir" and an example for WhereToSearch would be "C:\\Messy Dir\\". "." denotes the current directory. The names of all paths containing such a file or directory are printed and if exactly one such file or one directory is found, its path is returned.'''
    if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
        print('Searching:', Name)
    Filecounter=0
    Directorycounter=0
    ListOfPathsContainingThisName=[]
    for dirpath, subdirnames, filenames in os.walk(WhereToSearch):
        if Name in filenames:
            Filecounter+=1
            ListOfPathsContainingThisName.append(dirpath)
        if Name in subdirnames:
            Directorycounter+=1
            ListOfPathsContainingThisName.append(dirpath)
    if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
        print("    Finished with finding", Filecounter+Directorycounter ,"result(s) in the following folder(s):")
    for Pathname in ListOfPathsContainingThisName:
        if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
            print('   ', os.path.abspath(Pathname))
    if (Filecounter+Directorycounter)==1:
        return os.path.abspath(ListOfPathsContainingThisName[0])
    elif (Filecounter+Directorycounter)==0:
        if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
            print('    No such file with this filename.')
        return None
    elif (Filecounter+Directorycounter)>1:
        if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
            print('    Ambiguous occurrence of this filename.')
        return 'Ambiguous occurrence of this filename'

def SearchAFile0rDirectoryEverywhere(Filename,GoThere=True):
    '''This searches the directory that is containing the file or directory Filename. Filename has to be a string including the filetype-suffix if it is a file. Filename needn't be located in subsidiary directories it can also be located somewhere in superordinate directories. If GoThere is True, the working directory is moved to the found path, else this path is just returned.'''
    PathContainingFilename = SearchAFile0rDirectory(Filename,".") # This is None, if the searched file wasn't found among the current working directory. It contains the path of the file if it was found.
    while PathContainingFilename==None: # This loop breaks as soon as the file was found.
        if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
            print('    Searching in the superordinate directory...')
        os.chdir('..') # Go one directory upwards...
        PathContainingFilename = SearchAFile0rDirectory(Filename,".") # ...and search again.
    if PathContainingFilename =='Ambiguous occurrence of this filename':
        raise ImportFailedError
    if GoThere==True:
        os.chdir(PathContainingFilename)
        if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
            print('    Changed to this directory.')
    elif GoThere==False:
        pass
    return PathContainingFilename

SearchAFile0rDirectoryEverywhere("part7.py",GoThere=True)
from part6 import * # Import of the file.


## Test the terms of equation 1 

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n\n----------------------------- Consider terms of equation 1 ------------------------------\n')
    if IncludeSynchrotron==False:
        print('Test for ThirdTermOfEq1:', ThirdTermOfEq1(0.03/x0,Usedn0,None,None,CurrentNElectrons,1,CurrentLastValuesForgaIt),'\n')
        print('Test for IntegralOfBracketsOfEq1:', IntegralOfBracketsOfEq1(10000.0,Usedn0,None,None,CurrentNElectrons,0,CurrentLastValuesForgaIt,RelativeError=0.01),'\n')
        print('Test for AuxiliaryIntegrate:', AuxiliaryIntegrateIntegrandOf4thTerm(IntegrandOfFirstSummandOf4thTerm, 700000, 800000, 400000, Usedn0, None, None, UsedDotni, 0, CurrentLastValuesForgaIt, 0.01, NElectrons=CurrentNElectrons),'\n')
        if UseMultiprocessingIn4thTermOfEq1==False:
            #print('\nTest for FourthTermOfEq1 with function-object:', FourthTermOfEq1(400000.0,Usedn0,None,None,CurrentNElectrons,UsedDotni,0,CurrentLastValuesForgaIt,RelativeError=0.01),'\n')
            print('Test for FourthTermOfEq1 with interp1d-object:', FourthTermOfEq1(400000.0,Usedn0,None,None,interp1d(CurrentLastValuesForgaIt,[StartingNElectrons(ga) for ga in CurrentLastValuesForgaIt], kind='linear', bounds_error=False,fill_value=0.0),UsedDotni,1,CurrentLastValuesForgaIt,ValuesForNextNElectrons=[StartingNElectrons(ga) for ga in CurrentLastValuesForgaIt],RelativeError=0.01),'\n')
    elif IncludeSynchrotron==True:
        print('Test for ThirdTermOfEq1:', ThirdTermOfEq1(0.03/x0,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons,1,CurrentLastValuesForgaIt),'\n')
        print('Test for IntegralOfBracketsOfEq1:', IntegralOfBracketsOfEq1(10000.0,Usedn0,nSyncPs,ValuesForxSyncSuperIt,CurrentNElectrons,0,CurrentLastValuesForgaIt),'\n')
        print('Test for AuxiliaryIntegrate:', AuxiliaryIntegrateIntegrandOf4thTerm(IntegrandOfFirstSummandOf4thTerm, 700000, 800000, 400000, Usedn0, ValuesForxSyncSuperIt, ValuesFornSyncPs, UsedDotni, 0, CurrentLastValuesForgaIt, 0.01, NElectrons=CurrentNElectrons),'\n')
        if UseMultiprocessingIn4thTermOfEq1==False:
            print('\nTest for FourthTermOfEq1 with function-object:', FourthTermOfEq1(400000.0,Usedn0,ValuesForxSyncSuperIt,ValuesFornSyncPs,CurrentNElectrons,UsedDotni,0,CurrentLastValuesForgaIt),'\n')
            print('Test for FourthTermOfEq1 with interp1d-object:', FourthTermOfEq1(400000.0,Usedn0,ValuesForxSyncSuperIt,ValuesFornSyncPs,interp1d(CurrentLastValuesForgaIt,[StartingNElectrons(ga) for ga in CurrentLastValuesForgaIt], kind='linear', bounds_error=False,fill_value=0.0),UsedDotni,1,CurrentLastValuesForgaIt,ValuesForNextNElectrons=StartingNElectrons(CurrentLastValuesForgaIt)),'\n')
    print('Test for SynchrotronSpectralProductionRate:', SynchrotronSpectralProductionRate(gamma*100,CurrentNElectrons),'1/(m^3*s)\n')
    print('Test for SynchrotronSpectralLossRate:', SynchrotronSpectralLossRate(gamma),'1/s')


## Considering the first term of the numerator of equation 8

# The integrand of the first term of the numerator of equation 8 is just DotNi(gaP). It is in units of 1/(s*m^3).
# Now, it is to be integrated over gaP, from ga to infinity. According to my analysis, the resulting quantity is the injection-rate of particles with energy above ga per unit volume. Its units are still 1/(s*m^3).
def FirstTermOfNumeratorOfEq8(ga,DotNi):
    if ga < Injectedga0:
        if ga <= Injectedga1:
            LeftBorder = Injectedga1
        else:
            LeftBorder = ga
        FirstTermOfNumeratorOfEq8Result = integrate.quad(DotNi, LeftBorder, Injectedga0)[0]
        return FirstTermOfNumeratorOfEq8Result
    else:
        return 0.0

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    def EvaluateFirstTermOfNumeratorOfEq8():
        ValuesForgammaLin = np.linspace(0.8*Injectedga1,1.2*Injectedga0,1000)
        ValuesForFirstTermOfNumeratorOfEq8 = [FirstTermOfNumeratorOfEq8(i,UsedDotNi) for i in ValuesForgammaLin]
        pl.figure(figsize=(12, 9), num="Thomson-solution: First term of the numerator of equation 8 versus final electron energy in the case Dot N_i=%s" % UsedDotNi.__name__)
        pl.plot(ValuesForgammaLin,ValuesForFirstTermOfNumeratorOfEq8)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)


## Considering the second term of the numerator of equation 8

# Comment of v2_4: This computation of the second term of the numerator of eq. 8 was not used since the introduction of "PointsFromRToLIterationPointwise" any more and was not maintained as it seems outdated now.

def AuxiliaryFunctionForSecondTermOfNumeratorOfEq8(n0, Dotni, IterationCounter, LastValuesForgaIt, ValuesForNextNElectrons, LeftBorder, RightBorder,RelativeError=integrate.quad.__defaults__[3]):
    # This function will be executed in the child-processes. This function is necessary, because the interpolation and the integration have to be done inside the child-processes.
    StartingTime = time.time() # To measure the spent time of the computation, determine the point of time before the start of the computation.
    print('\nDoing internal interpolation.')
    InternalInterpolatedObjectForValuesForNextNElectrons = interp1d(LastValuesForgaIt, ValuesForNextNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
    print('Doing internal Assignment.')
    NElectrons = InternalInterpolatedObjectForValuesForNextNElectrons
    print('Integrating from gamma = %s to %s.' % (LeftBorder,RightBorder))
    AuxiliaryFunctionForSecondTermOfNumeratorOfEq8Result = integrate.quad(FourthTermOfEq1, LeftBorder, RightBorder, args=(n0,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,ValuesForNextNElectrons,RelativeError), epsrel=RelativeError)
    EndingTime = time.time() # Now, determine the point of time after the end of the integration.
    TimeInterval = (EndingTime-StartingTime)/3600.0 # Determine the time in hours, that was spent by the computation. 
    return (AuxiliaryFunctionForSecondTermOfNumeratorOfEq8Result, (LeftBorder, RightBorder, TimeInterval))

StorageOfSecondTermOfNumeratorOfEq8 = None # This is a storage for the value of the following SecondTermOfNumeratorOfEq8. As soon as SecondTermOfNumeratorOfEq8 was computed once, the result is saved in this quantity.

# The integrand of the second term of the numerator of equation 8 is FourthTermOfEq1(ga,n0,NElectrons,Dotni,IterationCounter,LastValuesForgaIt). It is in units of 1/(s*m^3). It gives the number-density of particles, that are pair-produced per unit time and per unit energy-interval. 
# Now, it is to be integrated over gaP, from 1/(4*x0) to infinity. According to my analysis, the resulting quantity is the injection-rate of pair-produced particles per unit volume. Its unit is still 1/(s*m^3).
def SecondTermOfNumeratorOfEq8(n0,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,ValuesForNextNElectrons,RelativeError=integrate.quad.__defaults__[3]):
    global StorageOfSecondTermOfNumeratorOfEq8 # Get access to the global quantity.
    if StorageOfSecondTermOfNumeratorOfEq8 == None: # Compute the integral, if it was not computed before.
        LowestIntegrationBorder=1.0/(4.0*x0) # The lower border of the integration range actually is 1/(4*x0).
        BiggestIntegrationBorder=UpperCutOffOf4thTermOfEq1
        NumberOfBorders = max(5,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
        if NumberOfIteratingProcesses>1: # If multiprocessing is applied, use more integration ranges.
            NumberOfBorders *= 5
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        ListOfIntegrationBorders = np.sort(np.append(ListOfIntegrationBorders,(ListOfIntegrationBorders[0]+ListOfIntegrationBorders[1])/2.0)) # Slightly above LowestIntegrationBorder, the integrand FourthTermOfEq1 rises very fast and therefore the integration is difficult here. Thus, a point between the first two integration borders is added to the list.
        ListOfIntegrationBorders = np.sort(np.append(ListOfIntegrationBorders,1.0/(2.0*x0))) # The point 1/(2*x0) is included into the range. This makes the integration a lot faster, as the integrand FourthTermOfEq1 has a kink here.
        if NumberOfIteratingProcesses==1: # Do not use multiprocessing:
            SecondTermOfNumeratorOfEq8Result = 0.0
            SecondTermOfNumeratorOfEq8BordersAndTiming = []
            for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
                StartingTime = time.time() # To measure the spent time of the computation, determine the point of time before the start of the computation.
                print('Integrating from gamma = %s to %s.' % (LeftBorder,RightBorder))
                SecondTermOfNumeratorOfEq8Result += integrate.quad(FourthTermOfEq1, LeftBorder, RightBorder, args=(n0,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,ValuesForNextNElectrons,RelativeError), epsrel=RelativeError)[0]
                EndingTime = time.time() # Now, determine the point of time after the end of the integration.
                TimeInterval = (EndingTime-StartingTime)/3600.0 # Determine the time in hours, that was spent by the computation. 
                SecondTermOfNumeratorOfEq8BordersAndTiming.append((LeftBorder, RightBorder, TimeInterval))
        else: # Use multiprocessing:
            PoolOfWorkers = multiprocessing.Pool(processes=NumberOfIteratingProcesses) # Initialise the multiple processes.
            QueueOfPartresultsOfSecondTermOfNumeratorOfEq8 = [PoolOfWorkers.apply_async(AuxiliaryFunctionForSecondTermOfNumeratorOfEq8, args=[n0, Dotni, IterationCounter, LastValuesForgaIt, ValuesForNextNElectrons, LeftBorder, RightBorder,RelativeError]) for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:])] # Distribute the tasks and start the work.
            ListOfPartresultsOfSecondTermOfNumeratorOfEq8 = [r.get() for r in QueueOfPartresultsOfSecondTermOfNumeratorOfEq8] # Collect the results into a list.
            PoolOfWorkers.close() # Clear the Pool-object to exit the processes and to...
            PoolOfWorkers.join() # ...retrieve the used RAM.
            SecondTermOfNumeratorOfEq8Result = sum([i[0][0] for i in ListOfPartresultsOfSecondTermOfNumeratorOfEq8]) # Take only the integration-results (and not the error) and add up all the part-results.
            SecondTermOfNumeratorOfEq8BordersAndTiming = [i[1] for i in ListOfPartresultsOfSecondTermOfNumeratorOfEq8] # Take the duration measurements and the corresponding borders of the integration-range.
        print("Time in hours spent on the integration of partresult of SecondTermOfNumeratorOfEq8 and corresponding integration-borders:")
        print("         Left border                  Right border             Time in hours")
        for i in SecondTermOfNumeratorOfEq8BordersAndTiming:
            print("%20.5f          %20.5f          %16.3f" % (i[0], i[1], i[2]))
        StorageOfSecondTermOfNumeratorOfEq8 = SecondTermOfNumeratorOfEq8Result # To pretend repetitive computation of the same quantity, save it in the global constant.
        return SecondTermOfNumeratorOfEq8Result
    else: # In case StorageOfSecondTermOfNumeratorOfEq8 already has a value:
        print("    SecondTermOfNumeratorOfEq8 was already determined. It is:", StorageOfSecondTermOfNumeratorOfEq8)
        return StorageOfSecondTermOfNumeratorOfEq8 # ... return the value that was determined before.
    # Comment of version 2: IterationCounter was added to the function arguments.
    # Comment of version 3: ValuesForgaLogComputation was added to the function arguments.
    # Comment of version 7: LowestIntegrationBorder was lowered to exactly 1.0/(4*x0).


## Considering the denominator of equation 8

# Now, compute the energy-loss-rate, which is given in units of 1/s.

CoefficientOfDenominatorOfEq8 = (4.0/3.0)*sigmaT*c # Abbreviation for the following function.

if Usedn0==n0Exp or Usedn0==n0PL or Usedn0==n0Planck or Usedn0 == n0Delta or Usedn0 == n0MultiDelta:
    def DenominatorOfEq8(ga,n0):
        DenominatorOfEq8Result = CoefficientOfDenominatorOfEq8*(ga**2.0)*IntegralOfDenominatorOfEq8(n0)
        return DenominatorOfEq8Result
elif Usedn0 == n0AF:
    def DenominatorOfEq8(ga,n0):
        DenominatorOfEq8Result = CoefficientOfDenominatorOfEq8*(ga**2.0)*(AFEnergyDensity(ADAFx1,ADAFxBremsstrahlungCutOff)/(me*c**2)) # AFEnergyDensity is in units of J/m^3, thus it has to be divided by me*c^2 again.
        return DenominatorOfEq8Result


## Considering all terms of equation 8

# Comment of v2_4: The computations in this chapter, too, were not used since the introduction of "PointsFromRToLIterationPointwise" any more and were not maintained as they seem outdated now.

def NElectronsOfEq8(ga,n0,NElectrons,Dotni,DotNi,IterationCounter,LastValuesForgaIt,ValuesForNextNElectrons):
    '''Determining N from equation 8:'''
    Term1 = FirstTermOfNumeratorOfEq8(ga,DotNi)
    Term2 = SecondTermOfNumeratorOfEq8(n0,NElectrons,Dotni,IterationCounter,LastValuesForgaIt,ValuesForNextNElectrons,0.01) # Pay attention to the choice of RelativeError!
    print('    Numerator of Eq8:', Term1+Term2)
    Term3 = DenominatorOfEq8(ga,n0)
    return (Term1+Term2)/Term3

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    def EvaluateNElectronsOfEq8(CurrentNElectrons,CurrentIterationCounter,LastValuesForgaIt,ValuesForNextNElectrons):
        '''Evaluate equation 8 based on Usedn0, CurrentNElectrons, UsedDotni and UsedDotNi.'''
        HighestgaOfEq8 = 1.0/(8.0*x0)
        LowestgaOfEq8 = 0.001*HighestgaOfEq8
        NumberOfSamplingPointsOfEq8 = 200.0
        ValuesForgaLogOfEq8 = np.logspace(np.log10(LowestgaOfEq8),np.log10(HighestgaOfEq8),NumberOfSamplingPointsOfEq8)
        ValuesForNElectronsOfEq8 = np.asarray([NElectronsOfEq8(i,Usedn0,CurrentNElectrons,UsedDotni,UsedDotNi,CurrentIterationCounter,LastValuesForgaIt,ValuesForNextNElectrons) for i in ValuesForgaLogOfEq8])
        # Append the Thomson-regime-values to the KN-regime-values...:
        ValuesForgaLogComplete = np.append(ValuesForgaLogOfEq8,LastValuesForgaIt)
        ValuesForNElectronsComplete = np.append(ValuesForNElectronsOfEq8,ValuesForNextNElectrons)
        # ... and again, create an interpolated object, that is assigned to CurrentNElectrons in order to being able to evaluate equation 5 over the Thomson- and KN-regime.
        InterpolatedObjectForValuesForNElectronsComplete = interp1d(ValuesForgaLogComplete, ValuesForNElectronsComplete, kind='linear', bounds_error=False, fill_value=0.0)   # If fill_value=0.0 is included, the value of the interpolating object outside the interpolation range is 0.0.
        CurrentNElectrons = InterpolatedObjectForValuesForNElectronsComplete
        return CurrentNElectrons, ValuesForgaLogComplete, ValuesForNElectronsComplete
    
    def EvaluateExportDataComplete(ValuesForgaLogComplete, ValuesForNElectronsComplete, StorageOfSecondTermOfNumeratorOfEq8):
        '''The final data points and the ValueForSecondTermOfNumeratorOfEq8 are written as a list into a .dat-file.'''
        Outputfile = open("Run %s, N(gamma) versus particle energy from %s (complete, data).dat" % (RunIdentifier,GetCurrentDate()), 'w')
        Outputfile.write('%s\n' % list(ValuesForgaLogComplete))
        Outputfile.write('%s\n' % list(ValuesForNElectronsComplete))
        Outputfile.write('%s' % StorageOfSecondTermOfNumeratorOfEq8)
        Outputfile.close()
    
    def EvaluateImportDataComplete():
        '''The string of a data-file is imported here to extract ValuesForgaLogComplete, ValuesForNElectronsComplete and the ValueForSecondTermOfNumeratorOfEq8. ValuesForNElectronsComplete is interpolated as above to restore CurrentNElectrons. The name of the file has to be inserted.'''
        OpenedLines = ImportFileToString('N(gamma) versus particle energy from 2018-06-15 20-03 (complete, data).dat')
        ValuesForgaLogCompleteString = OpenedLines[0].replace("\n", "") # Get the list ValuesForgaLogComplete.
        ValuesForgaLogComplete = np.asarray(eval(ValuesForgaLogCompleteString))
        ValuesForNElectronsCompleteString = OpenedLines[1].replace("\n", "") # Get the list ValuesForNElectronsComplete.
        ValuesForNElectronsComplete = np.asarray(eval(ValuesForNElectronsCompleteString))
        StorageOfSecondTermOfNumeratorOfEq8 = eval(OpenedLines[2]) # Get the current StorageOfSecondTermOfNumeratorOfEq8.
        InterpolatedObjectForValuesForNElectronsComplete = interp1d(ValuesForgaLogComplete, ValuesForNElectronsComplete, kind='linear', bounds_error=False, fill_value=0.0)   # The interpolation as above.
        CurrentNElectrons = InterpolatedObjectForValuesForNElectronsComplete
        return CurrentNElectrons, ValuesForgaLogComplete, ValuesForNElectronsComplete, StorageOfSecondTermOfNumeratorOfEq8
    
    def EvaluatePlotForNElectronsComplete(ValuesForgaLogComplete, ValuesForNElectronsComplete):
        '''Plotting of NElectrons in both the Thomson- and the KN-range.'''
        pl.figure(figsize=(18, 13), num="Cascade-equation: N(gamma) versus particle energy over complete range")
        # In figure 1 of Zdz., gamma*N(gamma) is plottet versus gamma*x0:
        pl.loglog(ValuesForgaLogComplete*x0,ValuesForgaLogComplete*ValuesForNElectronsComplete)
        ax = pl.gca()
        #pl.ylim(10**(1),10**(4))
        #pl.xlim(10**(-2),10**(0))
        ax.xaxis.set_tick_params(labelsize=26)
        ax.yaxis.set_tick_params(labelsize=26)
        pl.xlabel("Lorentz-factor times dim.less background photon energy $\gamma \cdot x_0$", family='serif', fontsize=30)
        pl.ylabel("Lorentz-factor times spectral\nnumber-density $\gamma \cdot N(\gamma)$ in $\mathrm{m}^{-3}$", family='serif', fontsize=32)
        #pl.savefig("Run %s, N(gamma) versus particle energy from %s (complete, plot).svg" % (RunIdentifier,GetCurrentDate()))
    
    def EvaluateSmearingOutTheStepViaFitting(ValuesForgaLogComplete, ValuesForNElectronsComplete, CurrentNElectrons, CurrentIterationCounter, LastValuesForgaIt):
        '''The step at gamma=1/(8*x0) is present, because neither the Thomson-limit nor the Klein-Nishina-limit is valid here. Hence, the true solution might be a mixture of both limiting cases. Thus, the discontinuity can be overcome by finding a function, that bridges the part around 1/(8*x0). This is done by fitting a power-law function with non-constant power-law-index to some data points. The data points, that are fitted, are taken partly from the Thomson-part and partly from the left edge of the Klein-Nishina-part.'''
        # First, select the range, out of which the fitting points are taken.
        LeftEdgeOfFittingRange = 0.04*1.0/(4.0*x0) # This might be varied.
        RightEdgeOfFittingRange = 1.9*1.0/(4.0*x0) # This might be a delicate parameter and might be varied!
        BooleanMaskForFittedRange = operator.and_(LeftEdgeOfFittingRange<=ValuesForgaLogComplete,ValuesForgaLogComplete<=RightEdgeOfFittingRange) # This array is True at those items, that correspond to data points, that lie within the range, that is to be fitted.
        SetOfIndicesOfPointsForFitting = [0,2,-6,-4,-3] # This list specifies the indices of the points, that are picked up (from the set of points, that are allowed by BooleanMaskForFittedRange) for being fitting. These might be delicate parameters and might be varied!
        gaPointsForFitting=ValuesForgaLogComplete[BooleanMaskForFittedRange][SetOfIndicesOfPointsForFitting] # Select the points out of ValuesForgaLogComplete ...
        NElectronsPointsForFitting=ValuesForNElectronsComplete[BooleanMaskForFittedRange][SetOfIndicesOfPointsForFitting] # ... and the corresponding points of ValuesForNElectronsComplete.
        # For the function, that is to be fitted to the data points, a smoothly broken power-law (sometimes called LogParabola) FittingParameter0 * x**(FittingParameter1+FittingParameter2*log10(x/LeftEdgeOfFittingRange)) is used.
        def ResidualsFunctionForFitting(FittingParameters, gaDataVector, NElectronsDataVector): # The residual function, whose sum of squares should be minimized. It takes as arguments a vector FittingParameters, which contains the fitting parameters, the gamma-data (x-data) and the NElectrons-data (y-data) points. It returns the residual vector.
            FittingParameter0, FittingParameter1, FittingParameter2 = FittingParameters
            ErrorVector = NElectronsDataVector - FittingParameter0 * gaDataVector**(FittingParameter1+FittingParameter2*np.log10(gaDataVector/LeftEdgeOfFittingRange))
            return ErrorVector
        # Define the starting positions for the fitting parameters:
        FittingParameter0Starting = NElectronsOfEq8(LeftEdgeOfFittingRange,Usedn0,CurrentNElectrons,UsedDotni,UsedDotNi,CurrentIterationCounter,LastValuesForgaIt,ValuesForNextNElectrons)/LeftEdgeOfFittingRange**(-2.0)
        FittingParameter1Starting = -2.1
        FittingParameter2Starting = -0.1
        FittingParametersStarting = [FittingParameter0Starting, FittingParameter1Starting, FittingParameter2Starting] # Collect them in a list.
        FittedParameters = optimize.leastsq(ResidualsFunctionForFitting, FittingParametersStarting, args=(gaPointsForFitting, NElectronsPointsForFitting))[0] # Fit it.
        print("    The parameters found via a least-squares fitting algorithm are:", FittedParameters)
        def FittedFunction(ga): # This gives the fitted function.
            return FittedParameters[0] * ga**(FittedParameters[1]+FittedParameters[2]*np.log10(ga/LeftEdgeOfFittingRange))
        # Plot the fitted function together with the original discontinuous function and the fitted data points:
        pl.figure(figsize=(15, 12), num="Cascade-equation: N(gamma) versus particle energy with fitted function (smoothing details)")
        pl.loglog(ValuesForgaLogComplete*x0,ValuesForNElectronsComplete)
        pl.plot(gaPointsForFitting*x0,NElectronsPointsForFitting,'o')
        gaFittedRange = np.linspace(gaPointsForFitting[0],gaPointsForFitting[-1],100)
        pl.plot(gaFittedRange*x0, FittedFunction(gaFittedRange))
        pl.xlim(0.9*LeftEdgeOfFittingRange*x0,1.5*RightEdgeOfFittingRange*x0)
        pl.ylim(0.1*NElectronsPointsForFitting[-1],2.0*NElectronsPointsForFitting[0])
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=26)
        ax.yaxis.set_tick_params(labelsize=26)
        pl.xlabel("Lorentz-factor times dim.less background photon energy $\gamma \cdot x_0$", family='serif', fontsize=28)
        pl.ylabel("Spectral number-density $N(\gamma)$", family='serif', fontsize=28)
        # At the Thomson-end of the fitted range the fitted function very precisely merges with the original ValuesForNElectronsComplete. However, at the Klein-Nishina-side the fitted function does not exactly merge with ValuesForNElectronsComplete but intersects it once or twice. The first intersection-point is searched now, by finding the first zero of the difference between FittedFunction and CurrentNElectrons right of the discontinuity.
        def AuxiliaryFunctionForFindingIntersection(ga):
            return FittedFunction(ga)-CurrentNElectrons(ga)
        gaIntersection = optimize.fsolve(AuxiliaryFunctionForFindingIntersection, 1.2/(4.0*x0))[0]
        print("    Intersection of FittedFunction and ValuesForNElectronsComplete is at gamma =", gaIntersection)
        pl.plot(gaIntersection*x0,FittedFunction(gaIntersection),'o')
        # Now, ValuesForNElectronsComplete is replaced by FittedFunction between LeftEdgeOfFittingRange and gaIntersection
        ValuesForgaLogReplacedRange = np.logspace(np.log10(LeftEdgeOfFittingRange),np.log10(gaIntersection),200)
        ValuesForNElectronsReplacedRange = FittedFunction(ValuesForgaLogReplacedRange)
        BooleanMaskForLeftNonfittedRange = ValuesForgaLogComplete<LeftEdgeOfFittingRange # This array is True at those items, that correspond to data points, that lie left outside of the fitted range.
        BooleanMaskForRightNonfittedRange = gaIntersection<ValuesForgaLogComplete # This array is True at those items, that correspond to data points, that lie right of gaIntersection.
        ValuesForgaLogSmoothed = np.concatenate((ValuesForgaLogComplete[BooleanMaskForLeftNonfittedRange],ValuesForgaLogReplacedRange,ValuesForgaLogComplete[BooleanMaskForRightNonfittedRange]))
        ValuesForNElectronsSmoothed = np.concatenate((ValuesForNElectronsComplete[BooleanMaskForLeftNonfittedRange],ValuesForNElectronsReplacedRange,ValuesForNElectronsComplete[BooleanMaskForRightNonfittedRange]))
        pl.plot(ValuesForgaLogSmoothed*x0,ValuesForNElectronsSmoothed)
        # Interpolate again:
        InterpolatedObjectForValuesForNElectronsSmoothed = interp1d(ValuesForgaLogSmoothed, ValuesForNElectronsSmoothed, kind='linear', bounds_error=False, fill_value=0.0)   # If fill_value=0.0 is included, the value of the interpolating object outside the interpolation range is 0.0.
        CurrentNElectrons = InterpolatedObjectForValuesForNElectronsSmoothed
        pl.savefig("Run %s, N(gamma) versus particle energy from %s (smoothing details).svg" % (RunIdentifier,GetCurrentDate()))
        return ValuesForgaLogSmoothed, ValuesForNElectronsSmoothed, CurrentNElectrons
    
    def EvaluatePlotForNElectronsSmoothed(ValuesForgaLogSmoothed, ValuesForNElectronsSmoothed):
        '''Plotting of the smoothed NElectrons in both the Thomson- and the KN-range.'''
        pl.figure(figsize=(18, 13), num="Cascade-equation: N(gamma) versus particle energy over complete range (smoothed)")
        # In figure 1 of Zdz., gamma*N(gamma) is plottet versus gamma*x0:
        pl.loglog(ValuesForgaLogSmoothed*x0,ValuesForgaLogSmoothed*ValuesForNElectronsSmoothed)
        ax = pl.gca()
        #pl.ylim(10.0**(1.0),10.0**(4.0))
        ax.xaxis.set_tick_params(labelsize=26)
        ax.yaxis.set_tick_params(labelsize=26)
        pl.xlabel("Lorentz-factor times dim.less background photon energy $\gamma \cdot x_0$", family='serif', fontsize=30)
        pl.ylabel("Lorentz-factor times spectral\nnumber-density $\gamma \cdot N(\gamma)$ in $\mathrm{m}^{-3}$", family='serif', fontsize=32)
        pl.savefig("Run %s, N(gamma) versus particle energy from %s (smoothed).svg" % (RunIdentifier,GetCurrentDate()))


## The energy-density of the electrons

def ElectronsSpectralEnergyDensity(ga,CurrentNElectrons):
    '''The spectral energy-density of the electrons in units of 1/m^3.'''
    return ga*CurrentNElectrons(ga)

def ElectronsEnergyDensity(LastValuesForgaIt,ValuesForNElectrons):
    '''The total energy-density of the electrons in units of 1/m^3.'''
    CurrentNElectrons = interp1d(LastValuesForgaIt, ValuesForNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
    ElectronsEnergyDensityResult = 0.0
    for LeftBorder, RightBorder in zip(LastValuesForgaIt[:-1], LastValuesForgaIt[1:]):
        ElectronsEnergyDensityResult += integrate.quad(ElectronsSpectralEnergyDensity, LeftBorder, RightBorder, args=(CurrentNElectrons))[0]
    return ElectronsEnergyDensityResult


## Considering equation 5 and the stationary HE-photon-spectrum:

# In all the following functions of this section, in the case IncludeSynchrotron==False, and one has to set nSyncPs=None and ValuesForxSyncSuperIt=None, while in the case IncludeSynchrotron==True, one has to set in the corresponding arrays.

# Equation 5 is equal to IntegralOfBracketsOfEq1. So, certain implementation-steps are not necessary, here.

# Determining the (spectral and total) power-density of the IC-upscattered photons:
def ICCreatedPowerDensitySpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the spectral energy-density injection-rate (= spectral power-density) of the inverse-Compton-scattered photons in units of 1/(m^3*s).'''
    return xga*IntegralOfBracketsOfEq1(xga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError)

def ICCreatedPowerDensity(n0,IterationCounter,LowestIntegrationBorder,LastValuesForgaIt,ValuesForNElectrons,ValuesForxSyncSuperIt,ValuesFornSyncPs,RelativeError=integrate.quad.__defaults__[3]):
    '''This is the total energy-density injection-rate (= total power-density) of the inverse-Compton-scattered photons in units of 1/(m^3*s).'''
    NElectrons = interp1d(LastValuesForgaIt, ValuesForNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
    if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray:
        nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
    else:
        nSyncPs = None
    BiggestIntegrationBorder=DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,LastValuesForgaIt[-1])
    if LowestIntegrationBorder>=BiggestIntegrationBorder:
        ICCreatedPowerDensityResult = 0.0 
    else:
        NumberOfBorders = max(10,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        ICCreatedPowerDensityResult = 0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            ICCreatedPowerDensityResult += integrate.quad(ICCreatedPowerDensitySpectralInxga, LeftBorder, RightBorder, args=(n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError), epsrel=RelativeError)[0]
    return ICCreatedPowerDensityResult

# The stationary HE-photon-spectrum can be determined via solving a corresponding kinetic equation in steady-state. The result is nHEPhotonsSumSpectralInxga = (Equation5+Dotni)/HEPhotonLossRate = Equation5/HEPhotonLossRate + Dotni/HEPhotonLossRate := nICCreatedSpectralInxga + nHEPhotonsSpectralInxga. Here, nHEPhotonsSumSpectralInxga is the spectral number-density of high-energetic photons in units of 1/m^3.
# Attention: This determination of nHEPhotonsSumSpectralInxga (and of all the quantities that are derived on it) makes sense, only if HE photon escape is included. Otherwise, HEPhotonLossRate would have to be substituted for by POfB11, which vanishes below the pair-production threshold. This again would cause nHEPhotonsSumSpectralInxga to be infinite below the threshold.
def nICCreatedSpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,IterationCounter,LastValuesForgaIt,NElectrons,RelativeError=integrate.quad.__defaults__[3]):
    '''The spectral number-density of high-energetic photons that are injected via IC-upscattering.'''
    return IntegralOfBracketsOfEq1(xga,n0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,LastValuesForgaIt,RelativeError)/HEPhotonLossRate(xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)

def ICCreatedEnergyDensitySpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,IterationCounter,LastValuesForgaIt,NElectrons,RelativeError=integrate.quad.__defaults__[3]):
    '''The corresponding spectral energy-density of high-energetic photons that are injected via IC-upscattering in units of 1/m^3.'''
    return xga*nICCreatedSpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,IterationCounter,LastValuesForgaIt,NElectrons,RelativeError)

def ICCreatedEnergyDensity(n0,IterationCounter,LastValuesForgaIt,ValuesForNElectrons,ValuesForxSyncSuperIt,ValuesFornSyncPs,RelativeError=integrate.quad.__defaults__[3]):
    '''The corresponding total (i. e. integrated) energy-density of high-energetic photons that are injected via IC-upscattering in units of 1/m^3.'''
    NElectrons = interp1d(LastValuesForgaIt, ValuesForNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
    if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray:
        nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
    else:
        nSyncPs = None
    # The integration essentially goes back to an integration of IntegralOfBracketsOfEq1/HEPhotonLossRate.
    LowestIntegrationBorder=10#0.00001 # This might be a delicate parameter.
    BiggestIntegrationBorder=DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,LastValuesForgaIt[-1]) # This is the upper cut-off of the integral in brackets of equation 1. 
    if LowestIntegrationBorder>=BiggestIntegrationBorder:
        ICCreatedEnergyDensityResult = 0.0
    # Now the integration: 
    else:
        NumberOfBorders = max(10,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
        ListOfIntegrationBorders = np.linspace(LowestIntegrationBorder,BiggestIntegrationBorder,NumberOfBorders)
        # There are jumps in the integrand due to the jumps in HEPhotonLossRate, so insert additional borders there.
        if Usedn0==n0MultiDelta:
            PotentialAdditionalBorders=1/x0MultiDelta
        else:
            PotentialAdditionalBorders=np.asarray([1/x0])
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,ValuesForxSyncSuperIt[-1]) # ValuesForxSyncSuperIt[-1] is xSync0Used.
        for PotentialBorder in PotentialAdditionalBorders:
            if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
                ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
        ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
        #print('    Integration borders of ICCreatedEnergyDensity:', ListOfIntegrationBorders)
        ICCreatedEnergyDensityResult = 0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            ICCreatedEnergyDensityResult += integrate.quad(ICCreatedEnergyDensitySpectralInxga, LeftBorder, RightBorder, args=(n0,nSyncPs,ValuesForxSyncSuperIt,IterationCounter,LastValuesForgaIt,NElectrons,RelativeError), epsrel=RelativeError)[0]
    return ICCreatedEnergyDensityResult

if UsedDotni==DotniDelta or UsedDotni==DotniZero: # Here, Dotni is just set =0 in the numerator due to its Delta-function-like nature or due to its vanishing.
    def nHEPhotonsSpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError=integrate.quad.__defaults__[3]):
        '''The spectral number-density of high-energetic photons that are directly injected.'''
        return 0.0
    def HEPhotonsEnergyDensitySpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError=integrate.quad.__defaults__[3]):
        '''The corresponding spectral energy-density of high-energetic photons that are directly injected in units of 1/m^3.'''
        return 0.0
    if UsedDotni==DotniDelta:
        def HEPhotonsEnergyDensity(n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError=integrate.quad.__defaults__[3]):
            '''The corresponding total (i. e. integrated) energy-density of high-energetic photons that are directly injected in units of 1/m^3. The expression in the second factor is the total number-density of the injected HE-photons.'''
            return HEPhotonsxga0Delta*DotniDelta/HEPhotonLossRate(HEPhotonsxga0Delta,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)
    elif UsedDotni==DotniZero:
        def HEPhotonsEnergyDensity(n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError=integrate.quad.__defaults__[3]):
            '''The corresponding total (i. e. integrated) energy-density of high-energetic photons that are directly injected in units of 1/m^3.'''
            return 0.0
else:
    def nHEPhotonsSpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError=integrate.quad.__defaults__[3]):
        '''The spectral number-density of high-energetic photons that are directly injected.'''
        return Dotni(xga)/HEPhotonLossRate(xga,n0,nSyncPs,ValuesForxSyncSuperIt,RelativeError)
    def HEPhotonsEnergyDensitySpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError=integrate.quad.__defaults__[3]):
        '''The corresponding spectral energy-density of high-energetic photons that are directly injected in units of 1/m^3.'''
        return xga*nHEPhotonsSpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError)
    def HEPhotonsEnergyDensity(n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError=integrate.quad.__defaults__[3]):
        '''The corresponding total (i. e. integrated) energy-density of high-energetic photons that are directly injected in units of 1/m^3.'''
        # The integration essentially goes back to an integration of Dotni/HEPhotonLossRate.
        LowestIntegrationBorder=HEPhotonsxga1NonDelta
        BiggestIntegrationBorder=HEPhotonsxga0NonDelta
        NumberOfBorders = max(10,round(np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)))
        ListOfIntegrationBorders = np.linspace(LowestIntegrationBorder,BiggestIntegrationBorder,NumberOfBorders)
        # There are jumps in the integrand due to the jumps in HEPhotonLossRate, so insert additional borders there.
        if Usedn0==n0MultiDelta:
            PotentialAdditionalBorders=1/x0MultiDelta
        else:
            PotentialAdditionalBorders=np.asarray([1/x0])
        if IncludeSynchrotron and PerformanceMode!='CommonIteration':
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,ValuesForxSyncSuperIt[-1]) # ValuesForxSyncSuperIt[-1] is xSync0Used.
        if InjectionType=='ResultsOfOldIt': # Add the kinks of Dotni.
            PotentialAdditionalBorders = np.append(PotentialAdditionalBorders,OldValuesForxgammanHEPhotons)
        for PotentialBorder in PotentialAdditionalBorders:
            if (LowestIntegrationBorder<PotentialBorder<BiggestIntegrationBorder) and (PotentialBorder not in ListOfIntegrationBorders):
                ListOfIntegrationBorders=np.append(ListOfIntegrationBorders,PotentialBorder)
        ListOfIntegrationBorders=np.sort(ListOfIntegrationBorders)
        #print('    Integration borders of HEPhotonsEnergyDensity:', ListOfIntegrationBorders)
        HEPhotonsEnergyDensityResult = 0.0
        for LeftBorder, RightBorder in zip(ListOfIntegrationBorders[:-1], ListOfIntegrationBorders[1:]):
            HEPhotonsEnergyDensityResult += integrate.quad(HEPhotonsEnergyDensitySpectralInxga, LeftBorder, RightBorder, args=(n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError), epsrel=RelativeError)[0]
        return HEPhotonsEnergyDensityResult

def nHEPhotonsSumSpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,IterationCounter,LastValuesForgaIt,NElectrons,RelativeError=integrate.quad.__defaults__[3]):
    '''The spectral number-density of all the high-energetic photons.'''
    return nICCreatedSpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,IterationCounter,LastValuesForgaIt,NElectrons,RelativeError)+nHEPhotonsSpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError)

def HEPhotonsSumEnergyDensitySpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,IterationCounter,LastValuesForgaIt,NElectrons,RelativeError=integrate.quad.__defaults__[3]):
    '''The corresponding spectral energy-density of all high-energetic photons in units of 1/m^3.'''
    return ICCreatedEnergyDensitySpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,IterationCounter,LastValuesForgaIt,NElectrons,RelativeError)+HEPhotonsEnergyDensitySpectralInxga(xga,n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError)

def HEPhotonsSumEnergyDensity(n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,IterationCounter,LastValuesForgaIt,NElectrons,RelativeError=integrate.quad.__defaults__[3]):
    '''The corresponding total (i. e. integrated) energy-density of all the high-energetic photons in units of 1/m^3.'''
    return ICCreatedEnergyDensity(n0,nSyncPs,ValuesForxSyncSuperIt,IterationCounter,LastValuesForgaIt,NElectrons,RelativeError)+HEPhotonsEnergyDensity(n0,nSyncPs,ValuesForxSyncSuperIt,Dotni,RelativeError)

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n\n-------------------------- Consider the injected HE-photons -----------------------------\n')
    if not(SuperItCounter): # The case where no synchrotron-radiation field has been determined yet via a super-iteration.
        CurrentHEPhotonsEnergyDensity = HEPhotonsEnergyDensity(Usedn0,None,None,UsedDotni,0.01)*me*c**2
        print("Injected HE-photons' energy-density (using no pair-absorption via SyncPs):", CurrentHEPhotonsEnergyDensity, "J/m^3")
        print("    Corresponding injected HE-photons' luminosity (using no pair-absorption via SyncPs):", CurrentHEPhotonsEnergyDensity*FactorEnergyDensityToLuminosity('InteractionRadius')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('InteractionRadius')[0])
        print("    Corresponding injected HE-photons' luminosity (using no pair-absorption via SyncPs):", CurrentHEPhotonsEnergyDensity*FactorEnergyDensityToLuminosity('RClouds')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('RClouds')[0])
    else: # The case where the synchrotron-radiation field has already been determined.
        CurrentHEPhotonsEnergyDensity = HEPhotonsEnergyDensity(Usedn0,nSyncPs,ValuesForxSyncSuperIt,UsedDotni,0.01)*me*c**2
        print("Injected HE-photons' energy-density (including pair-absorption via SyncPs):", CurrentHEPhotonsEnergyDensity, "J/m^3")
        print("    Corresponding injected HE-photons' luminosity (including pair-absorption via SyncPs):", CurrentHEPhotonsEnergyDensity*FactorEnergyDensityToLuminosity('InteractionRadius')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('InteractionRadius')[0])
        print("    Corresponding injected HE-photons' luminosity (including pair-absorption via SyncPs):", CurrentHEPhotonsEnergyDensity*FactorEnergyDensityToLuminosity('RClouds')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('RClouds')[0])
    CurrentHEPhotonsPowerDensity = HEPhotonsPowerDensity(UsedDotni)
    print("Injected HE-photons' power-density:", CurrentHEPhotonsPowerDensity, "1/(m^3*s)")
    print('    Corresponding injected HE-photon power:', SphericalVolume('InteractionRadius')[1]*CurrentHEPhotonsPowerDensity*me*c**2, "W (using R=%s)" % SphericalVolume('InteractionRadius')[0])
    print('    Corresponding injected HE-photon power:', SphericalVolume('RClouds')[1]*CurrentHEPhotonsPowerDensity*me*c**2, "W (using R=%s)" % SphericalVolume('RClouds')[0])

def EvaluateEquation5AndnHEPhotons(IterationCounter,ValuesForgaIt,ValuesForNElectrons,ValuesForxSyncSuperIt,ValuesFornSyncPs):
    '''Sample and plot x_gamma^2 * dot n_gamma and x_gamma^2 * n_gamma versus x_gamma. In the case IncludeSynchrotron==False, and one has to set nSyncPs=None and ValuesForxSyncSuperIt=None, while in the case IncludeSynchrotron==True, one has to set in the corresponding arrays.'''
    print('\nCalling EvaluateEquation5AndnHEPhotons:')
    NElectrons = interp1d(ValuesForgaIt, ValuesForNElectrons, kind='linear', bounds_error=False, fill_value=0.0)
    if type(ValuesForxSyncSuperIt)==np.ndarray and type(ValuesFornSyncPs)==np.ndarray:
        nSyncPs = interp1d(ValuesForxSyncSuperIt, ValuesFornSyncPs, kind='linear', bounds_error=False, fill_value=0.0)
    else:
        nSyncPs = None
    
    # One sampling-range of xga for the evaluation of functions based on nHEPhotons and one for the evaluation of Eq5:
    ValuesForxgammanHEPhotons = np.logspace(np.log10(ValuesForgaIt[0]),np.log10(DetermineBiggestIntegrationBorderOf4thTermOfEq1(ValuesForxSyncSuperIt,ValuesForgaIt[-1])),200)
    ValuesForxgammaEq5 = np.logspace(np.log10(ValuesForgaIt[0]),np.log10(DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,ValuesForgaIt[-1])),200)
    
    # Sample nICCreatedSpectralInxga, nHEPhotonsSumSpectralInxga and nHEPhotonsSpectralInxga: (ValuesForgaIt and LastValuesForgaIt is again set equal here)
    if UsedIterationScheme == 'IterationStepByStepPointsFromLToR':
        ValuesFornICCreatedSpectralInxga = np.asarray([nICCreatedSpectralInxga(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,IterationCounter,ValuesForgaIt,NElectrons,0.01) for i in ValuesForxgammanHEPhotons])
        ValuesFornHEPhotonsSumSpectralInxga = np.asarray([nHEPhotonsSumSpectralInxga(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,UsedDotni,IterationCounter,ValuesForgaIt,NElectrons,0.01) for i in ValuesForxgammanHEPhotons])
    elif UsedIterationScheme == 'PointsFromRToLIterationPointwise':
        ValuesFornICCreatedSpectralInxga = np.asarray([nICCreatedSpectralInxga(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,1,ValuesForgaIt,NElectrons,0.01) for i in ValuesForxgammanHEPhotons])
        ValuesFornHEPhotonsSumSpectralInxga = np.asarray([nHEPhotonsSumSpectralInxga(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,UsedDotni,1,ValuesForgaIt,NElectrons,0.01) for i in ValuesForxgammanHEPhotons])
    ValuesFornHEPhotonsSpectralInxga = np.asarray([nHEPhotonsSpectralInxga(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,UsedDotni,0.01) for i in ValuesForxgammanHEPhotons])
    
    # Sample IntegralOfBracketsOfEq1, i.e. the IC-created production-rate: (ValuesForgaIt and LastValuesForgaIt is again set equal here)
    if UsedIterationScheme == 'IterationStepByStepPointsFromLToR':
        ValuesForEq5 = np.asarray([IntegralOfBracketsOfEq1(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,IterationCounter,ValuesForgaIt,0.01) for i in ValuesForxgammaEq5])
    elif UsedIterationScheme == 'PointsFromRToLIterationPointwise':
        ValuesForEq5 = np.asarray([IntegralOfBracketsOfEq1(i,Usedn0,nSyncPs,ValuesForxSyncSuperIt,NElectrons,1,ValuesForgaIt,0.01) for i in ValuesForxgammaEq5])
    
    # Sample the injection-rate:
    if UsedDotni!=DotniZero and UsedDotni!=DotniDelta:
        ValuesForDotni = np.asarray([UsedDotni(i) for i in ValuesForxgammaEq5])
    
    # Prepare a plot that draws all lines in one figure with two y-axes:
    pl.rc('font', family='serif')
    Fig, LeftAxis = pl.subplots(figsize=(14, 12), num="Cascade-equation: Dotn(xgamma)*xgamma^2 and n(xgamma)*xgamma^2 versus final photon energy xgamma")
    LeftAxis.set_xlabel("$x_\gamma x_{1}$", fontsize=50)
    if UsedDotni==DotniZero:
        UpperLimitOfPlot = 1.2*DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,ValuesForgaIt[-1])*x0
    else:
        UpperLimitOfPlot = 1.2*DetermineBiggestIntegrationBorderOf4thTermOfEq1(ValuesForxSyncSuperIt,ValuesForgaIt[-1])*x0
    LeftAxis.set_xlim((1*ValuesForgaIt[0]*x0,UpperLimitOfPlot))
    LeftAxis.set_ylabel("$x_\gamma^2 \dot n_{\gamma, \, \mathrm{IC}} \, / \, \mathrm{m}^{-3} \mathrm{s}^{-1}$", fontsize=50, color='purple')
    
    # Plot x_gamma^2 * dot n_gamma:
    LeftAxis.loglog(ValuesForxgammaEq5*x0,ValuesForxgammaEq5**2.0*ValuesForEq5, label='$\dot n_{\gamma, \, \mathrm{IC}}$ (Spectral IC production rate)', color='purple', linewidth=3, linestyle='-') # Plot the IC-created production-rate.
    if UsedDotni!=DotniZero:
        if UsedDotni==DotniDelta:
            LeftAxis.loglog(np.asarray([HEPhotonsxga0Delta,HEPhotonsxga0Delta])*x0,np.asarray([10**(-30),10**(30)]), label='$\dot n_{\gamma, \, \mathrm{i}}$ (Spectral injection rate)', color='purple', linewidth=2, linestyle='--') # Plot the injection-rate.
            LeftAxis.set_ylim((0.5*min((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()]),2.0*max((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()]))) # Automated scaling.
        else:
            LeftAxis.loglog(ValuesForxgammaEq5*x0,ValuesForxgammaEq5**2.0*ValuesForDotni, label='$\dot n_{\gamma, \, \mathrm{i}}$ (Spectral injection rate)', color='purple', linewidth=2, linestyle='--') # Plot the injection-rate.
            LeftAxis.loglog(ValuesForxgammaEq5*x0,ValuesForxgammaEq5**2.0*(ValuesForEq5+ValuesForDotni), label='$\dot n_{\gamma, \, \mathrm{i}} + \dot n_{\gamma, \, \mathrm{IC}}$ (Sum)', color='purple', linewidth=2, linestyle='-') # Plot the summed production-rate.
            #LeftAxis.set_ylim((0.5*min([min((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()]),min((ValuesForxgammaEq5[ValuesForDotni.nonzero()])**2.0*ValuesForDotni[ValuesForDotni.nonzero()]),min((ValuesForxgammaEq5[(ValuesForEq5+ValuesForDotni).nonzero()])**2.0*(ValuesForEq5+ValuesForDotni)[(ValuesForEq5+ValuesForDotni).nonzero()])]),2.0*max([max((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()]),max((ValuesForxgammaEq5[ValuesForDotni.nonzero()])**2.0*ValuesForDotni[ValuesForDotni.nonzero()]),max((ValuesForxgammaEq5[(ValuesForEq5+ValuesForDotni).nonzero()])**2.0*(ValuesForEq5+ValuesForDotni)[(ValuesForEq5+ValuesForDotni).nonzero()])]))) # Automated scaling.
            LeftAxis.set_ylim(2.0*10**2,2.0*10**11) # Special for Mrk501-bump.
    else:
        LeftAxis.set_ylim((0.5*min((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()]),2.0*max((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()]))) # Automated scaling.
    LeftAxis.tick_params(axis='y', labelsize=40, direction='inout', width=3, length=14, labelcolor='purple', pad=-8)
    LeftAxis.tick_params(axis='y', which='minor', labelsize=22, direction='inout', width=2, length=10, labelcolor='purple')
    LeftAxis.tick_params(axis='x', labelsize=40, direction='inout', width=3, length=14, pad=-5)
    LeftAxis.tick_params(axis='x', which='minor', labelsize=22, direction='inout', width=2, length=10, pad=7)
    #LeftAxis.legend(loc="upper left", fontsize=18)
    LeftAxis.spines['left'].set_linewidth(3)
    LeftAxis.spines['right'].set_linewidth(3)
    LeftAxis.spines['top'].set_linewidth(3)
    LeftAxis.spines['bottom'].set_linewidth(3)
    
    # Instantiate a second axes that shares the same x-axis:
    RightAxis = LeftAxis.twinx()
    RightAxis.set_ylabel("$x_\gamma^2 n_{\gamma} \, / \, \mathrm{m}^{-3}$", fontsize=50, color='blue')
    RightAxis.set_xlim((1*ValuesForgaIt[0]*x0,UpperLimitOfPlot))
    
    # Plot x_gamma^2 * n_gamma:
    RightAxis.loglog(ValuesForxgammanHEPhotons*x0,ValuesForxgammanHEPhotons**2.0*ValuesFornICCreatedSpectralInxga, label='$n_{\gamma, \, \mathrm{IC}}$ (IC contribution)', color='blue', linewidth=3, linestyle='-') # Plot the IC-created number-density.
    if UsedDotni!=DotniZero and UsedDotni!=DotniDelta:
        RightAxis.loglog(ValuesForxgammanHEPhotons*x0,ValuesForxgammanHEPhotons**2.0*ValuesFornHEPhotonsSpectralInxga, label='$n_{\gamma, \, \mathrm{i}}$ (Injection contribution)', color='blue', linewidth=2, linestyle='--') # Plot the injected number-density.
        RightAxis.loglog(ValuesForxgammanHEPhotons*x0,ValuesForxgammanHEPhotons**2.0*ValuesFornHEPhotonsSumSpectralInxga, label='$n_\gamma$ (Sum)', color='blue', linewidth=2, linestyle='-') # Plot the summed number-density.
    RightAxis.set_ylim(3.0*10**11,1.0*10**14) # Special for Mrk501-bump.
    RightAxis.tick_params(axis='y', labelsize=40, direction='inout', width=3, length=14, labelcolor='blue')
    RightAxis.tick_params(axis='y', which='minor', labelsize=22, direction='inout', width=2, length=10, labelcolor='blue')
    #RightAxis.legend(loc="lower right", fontsize=18)
    # The following is to create a different top x-axis:
    UpperAxis=LeftAxis.twiny()
    UpperAxis.set_xscale('log')
    UpperAxis.set_xlim(LeftAxis.get_xlim())
    UpperAxisMajorTicks=np.asarray([10**(-4),10**(-3),10**(-2),10**(-1),10**(0)]) # The tick-labels in TeV.
    UpperAxisMinorTicks=np.asarray([8*10**(-5),9*10**(-5),2*10**(-4),3*10**(-4),4*10**(-4),5*10**(-4),6*10**(-4),7*10**(-4),8*10**(-4),9*10**(-4),2*10**(-3),3*10**(-3),4*10**(-3),5*10**(-3),6*10**(-3),7*10**(-3),8*10**(-3),9*10**(-3),2*10**(-2),3*10**(-2),4*10**(-2),5*10**(-2),6*10**(-2),7*10**(-2),8*10**(-2),9*10**(-2),2*10**(-1),3*10**(-1),4*10**(-1),5*10**(-1),6*10**(-1),7*10**(-1),8*10**(-1),9*10**(-1),2*10**(0),3*10**(0),4*10**(0),5*10**(0),6*10**(0)]) # The minor tick-labels in TeV.
    UpperAxis.set_xticks(UpperAxisMajorTicks*x0*e*10**12/(me*c**2)) # Conversion from values of the top axis to the positions on the bottom axis.
    UpperAxis.set_xticks(UpperAxisMinorTicks*x0*e*10**12/(me*c**2), minor=True) # Conversion from values of the top axis to the positions on the bottom axis.
    UpperAxis.set_xticklabels(np.asarray(['$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$']))
    UpperAxis.set_xticklabels([], minor=True)
    UpperAxis.xaxis.set_tick_params(labelsize=40, direction='inout', width=5, length=16, pad=8)
    UpperAxis.xaxis.set_tick_params(which='minor', labelsize=22, direction='inout', width=2, length=10, pad=8)
    UpperAxis.set_xlabel("$x_\gamma m_{\mathrm{e}} c^2 / \, \mathrm{TeV}$", fontsize=50, labelpad=20)
    #Fig.tight_layout()  # otherwise the right y-label is slightly clipped
    Fig.subplots_adjust(top=0.86, bottom=0.14, left=0.18, right=0.83)
    pl.savefig("Run %s, Dotn and n times xgamma^2 versus xgamma from %s.svg" % (RunIdentifier,GetCurrentDate()))
    pl.savefig("Run %s, Dotn and n times xgamma^2 versus xgamma from %s.pdf" % (RunIdentifier,GetCurrentDate()))

    # Prepare a plot that draws one figure for the rates:
    pl.rc('font', family='serif')
    Fig2, LeftAxis2 = pl.subplots(figsize=(13, 10), num="Cascade-equation: Dotn(xgamma)*xgamma^2 versus final photon energy xgamma")
    LeftAxis2.set_xlabel("$x_\gamma \cdot \, x_{0,1}$", fontsize=44)
    if UsedDotni==DotniZero:
        UpperLimitOfPlot = 1.2*DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,ValuesForgaIt[-1])*x0
    else:
        UpperLimitOfPlot = 1.2*DetermineBiggestIntegrationBorderOf4thTermOfEq1(ValuesForxSyncSuperIt,ValuesForgaIt[-1])*x0
    LeftAxis2.set_xlim((ValuesForgaIt[0]*x0,UpperLimitOfPlot))
    LeftAxis2.set_ylabel("$x_\gamma^2 \cdot \, \dot n_{\ldots}(x_\gamma) \, \, [\mathrm{m}^{-3} \mathrm{s}^{-1}]$", fontsize=44)
    
    # Plot x_gamma^2 * dot n_gamma:
    LeftAxis2.loglog(ValuesForxgammaEq5*x0,ValuesForxgammaEq5**2.0*ValuesForEq5, label='$x_\gamma^2 \dot n_{\gamma,\mathrm{IC}}$ (Inverse-Compton)', color='purple', linewidth=6, linestyle=':') # Plot the IC-created production-rate.
    if UsedDotni!=DotniZero:
        if UsedDotni==DotniDelta:
            LeftAxis2.loglog(np.asarray([HEPhotonsxga0Delta,HEPhotonsxga0Delta])*x0,np.asarray([10**(-40),10**(40)]), label='$x_\gamma^2 \dot n_{\gamma,\mathrm{i}}$ (Injection)', color='purple', linewidth=4, linestyle='--') # Plot the injection-rate.
            LeftAxis2.set_ylim((0.5*min((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()]),2.0*max((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()])))
        else:
            LeftAxis2.loglog(ValuesForxgammaEq5*x0,ValuesForxgammaEq5**2.0*ValuesForDotni, label='$x_\gamma^2 \dot n_{\gamma,\mathrm{i}}$ (Injection)', color='purple', linewidth=4, linestyle='--') # Plot the injection-rate.
            LeftAxis2.loglog(ValuesForxgammaEq5*x0,ValuesForxgammaEq5**2.0*(ValuesForEq5+ValuesForDotni), label='$x_\gamma^2 \dot n_\gamma$ (Sum)', color='purple', linewidth=3, linestyle='-') # Plot the summed production-rate.
            #LeftAxis2.set_ylim((0.5*min([min((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()]),min((ValuesForxgammaEq5[ValuesForDotni.nonzero()])**2.0*ValuesForDotni[ValuesForDotni.nonzero()]),min((ValuesForxgammaEq5[(ValuesForEq5+ValuesForDotni).nonzero()])**2.0*(ValuesForEq5+ValuesForDotni)[(ValuesForEq5+ValuesForDotni).nonzero()])]),2.0*max([max((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()]),max((ValuesForxgammaEq5[ValuesForDotni.nonzero()])**2.0*ValuesForDotni[ValuesForDotni.nonzero()]),max((ValuesForxgammaEq5[(ValuesForEq5+ValuesForDotni).nonzero()])**2.0*(ValuesForEq5+ValuesForDotni)[(ValuesForEq5+ValuesForDotni).nonzero()])])))
            LeftAxis2.set_ylim(2.0*10**8,3.0*10**11)
    else:
        LeftAxis2.set_ylim(2.0*10**8,3.0*10**11)#LeftAxis2.set_ylim((0.5*min((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()]),2.0*max((ValuesForxgammaEq5[ValuesForEq5.nonzero()])**2.0*ValuesForEq5[ValuesForEq5.nonzero()])))
    LeftAxis2.tick_params(axis='y', labelsize=38, direction='inout', width=5, length=20, pad=10)
    LeftAxis2.tick_params(axis='y', which='minor', labelsize=32, direction='inout', width=3, length=10)
    LeftAxis2.tick_params(axis='x', labelsize=38, direction='inout', width=5, length=20, pad=0, top=False)
    LeftAxis2.tick_params(axis='x', which='minor', labelsize=32, direction='inout', width=3, length=10, top=False)
    
    # Formatting:
    LeftAxis2.spines['left'].set_linewidth(4)
    LeftAxis2.spines['right'].set_linewidth(4)
    LeftAxis2.spines['top'].set_linewidth(4)
    LeftAxis2.spines['bottom'].set_linewidth(4)
    LeftAxis2.legend(loc="lower left", fontsize=30)
    
    # The following is to create a different top x-axis:
    UpperAxis2=LeftAxis2.twiny()
    UpperAxis2.set_xscale('log')
    UpperAxis2.set_xlim(LeftAxis2.get_xlim())
    UpperAxisTicks2=np.asarray([10**(-1),10**(0),10**(1),10**(2)]) # The tick-labels.
    UpperAxis2.set_xticks(UpperAxisTicks2*x0*e*10**9/(me*c**2)) # Conversion from values of the top axis to the positions on the bottom axis).
    UpperAxis2.set_xticklabels(UpperAxisTicks2)
    UpperAxis2.xaxis.set_tick_params(labelsize=35, direction='inout', width=5, length=16, pad=8)
    UpperAxis2.set_xlabel("$x_\gamma \cdot \, m_{\mathrm{e}} c^2 \, \, [{\mathrm{GeV}}]$", fontsize=42, labelpad=20)
    #Fig.tight_layout()  # otherwise the right y-label is slightly clipped
    Fig2.subplots_adjust(top=0.85, bottom=0.16, left=0.18, right=0.97)
    #pl.savefig("Run %s, Dotn times xgamma^2 versus final photon energy from %s.svg" % (RunIdentifier,GetCurrentDate()))
    pl.savefig("Run %s, Dotn times xgamma^2 versus final photon energy from %s.pdf" % (RunIdentifier,GetCurrentDate()))

    # Prepare a plot that draws one figure for the densities:
    pl.rc('font', family='serif')
    Fig3, LeftAxis3 = pl.subplots(figsize=(13, 10), num="Cascade-equation: n(xgamma)*xgamma^2 versus final photon energy xgamma")
    LeftAxis3.set_xlabel("$x_\gamma \cdot \, x_{0,1}$", fontsize=44)
    if UsedDotni==DotniZero:
        UpperLimitOfPlot = 1.2*DetermineUpperCutOffOfEq5(ValuesForxSyncSuperIt,ValuesForgaIt[-1])*x0
    else:
        UpperLimitOfPlot = 1.2*DetermineBiggestIntegrationBorderOf4thTermOfEq1(ValuesForxSyncSuperIt,ValuesForgaIt[-1])*x0
    LeftAxis3.set_xlim((ValuesForgaIt[0]*x0,UpperLimitOfPlot))
    LeftAxis3.set_ylabel("$x_\gamma^2 \cdot \, n_{\ldots}(x_\gamma) \, \, [\mathrm{m}^{-3}]$", fontsize=44)
    
    # Plot x_gamma^2 * n_gamma:
    LeftAxis3.loglog(ValuesForxgammanHEPhotons*x0,ValuesForxgammanHEPhotons**2.0*ValuesFornICCreatedSpectralInxga, label='$x_\gamma^2 n_{\gamma,\mathrm{IC}}$ (Inverse-Compton)', color='blue', linewidth=6, linestyle=':') # Plot the IC-created number-density.
    if UsedDotni!=DotniZero and UsedDotni!=DotniDelta:
        LeftAxis3.loglog(ValuesForxgammanHEPhotons*x0,ValuesForxgammanHEPhotons**2.0*ValuesFornHEPhotonsSpectralInxga, label='$x_\gamma^2 n_{\gamma,\mathrm{i}}$ (Injection)', color='blue', linewidth=4, linestyle='--') # Plot the injected number-density.
        LeftAxis3.loglog(ValuesForxgammanHEPhotons*x0,ValuesForxgammanHEPhotons**2.0*ValuesFornHEPhotonsSumSpectralInxga, label='$x_\gamma^2 n_\gamma$ (Sum)', color='blue', linewidth=3, linestyle='-') # Plot the summed number-density.
    LeftAxis3.set_ylim(3.0*10**11,1.0*10**14)#LeftAxis3.set_ylim(10**6,10**11)#LeftAxis3.set_ylim(3.0*10**7,5.0*10**9)
    LeftAxis3.tick_params(axis='y', labelsize=38, direction='inout', width=5, length=20, pad=10)
    LeftAxis3.tick_params(axis='y', which='minor', labelsize=32, direction='inout', width=3, length=10)
    LeftAxis3.tick_params(axis='x', labelsize=38, direction='inout', width=5, length=20, pad=0, top=False)
    LeftAxis3.tick_params(axis='x', which='minor', labelsize=32, direction='inout', width=3, length=10, top=False)
    
    # Formatting:
    LeftAxis3.spines['left'].set_linewidth(4)
    LeftAxis3.spines['right'].set_linewidth(4)
    LeftAxis3.spines['top'].set_linewidth(4)
    LeftAxis3.spines['bottom'].set_linewidth(4)
    LeftAxis3.legend(loc="lower left", fontsize=30)
    
    # The following is to create a different top x-axis:
    UpperAxis3=LeftAxis3.twiny()
    UpperAxis3.set_xscale('log')
    UpperAxis3.set_xlim(LeftAxis3.get_xlim())
    UpperAxisTicks3=np.asarray([10**(-1),10**(0),10**(1),10**(2)]) # The tick-labels.
    UpperAxis3.set_xticks(UpperAxisTicks3*x0*e*10**9/(me*c**2)) # Conversion from values of the top axis. to the positions )on the bottom axis).
    UpperAxis3.set_xticklabels(UpperAxisTicks3)
    UpperAxis3.xaxis.set_tick_params(labelsize=35, direction='inout', width=5, length=16, pad=8)
    UpperAxis3.set_xlabel("$x_\gamma \cdot \, m_{\mathrm{e}} c^2 \, \, [{\mathrm{GeV}}]$", fontsize=42, labelpad=20)
    #Fig.tight_layout()  # otherwise the right y-label is slightly clipped
    Fig3.subplots_adjust(top=0.85, bottom=0.16, left=0.18, right=0.97)
    #pl.savefig("Run %s, n times xgamma^2 versus final photon energy from %s.svg" % (RunIdentifier,GetCurrentDate()))
    pl.savefig("Run %s, n times xgamma^2 versus final photon energy from %s.pdf" % (RunIdentifier,GetCurrentDate()))
    
    # Save the number-density:
    Outputfile = open("Run %s, nHEP versus energy from %s.dat" % (RunIdentifier,GetCurrentDate()), 'w')
    Outputfile.write('Values for the energy of the high-energetic photons in units of m_e*c^2, Corresponding values for the sum of the number-density in m^(-3), Corresponding values for the number-density of the injected HEPs in m^(-3)\n')
    for i in range(len(ValuesForxgammanHEPhotons)):
        Outputfile.write('%s %s %s\n' % (ValuesForxgammanHEPhotons[i],ValuesFornHEPhotonsSumSpectralInxga[i],ValuesFornHEPhotonsSpectralInxga[i]))
    Outputfile.close()

    # Plot dot n_gamma versus x_gamma:
    pl.figure(figsize=(12, 11), num="Cascade-equation: Dotn(xgamma) versus final photon energy xgamma")
    pl.loglog(ValuesForxgammaEq5*x0,ValuesForEq5, label='$\dot n_{\mathrm{IC}}$ (Inverse-Compton)')
    if UsedDotni!=DotniZero and UsedDotni!=DotniDelta:
        pl.loglog(ValuesForxgammaEq5*x0,ValuesForDotni, label='$\dot n_{\mathrm{i}}$ (Injection)')
        pl.loglog(ValuesForxgammaEq5*x0,ValuesForEq5+ValuesForDotni, label='$\dot n_\gamma$ (Sum)')
    ax = pl.gca()
    ax.xaxis.set_tick_params(labelsize=26)
    ax.yaxis.set_tick_params(labelsize=26)
    pl.legend(loc="best", fontsize=22)
    pl.xlabel("Energy $x \cdot x_0$", fontsize=32)
    pl.ylabel("HEP spectral production rate $\dot n_\gamma(x_\gamma)$", fontsize=32)
    pl.subplots_adjust(top=0.97, bottom=0.11, left=0.15, right=0.96)
    pl.savefig("Run %s, Dotn versus xgamma from %s.svg" % (RunIdentifier,GetCurrentDate()))
    
    return ValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, ValuesFornHEPhotonsSpectralInxga


## Obtaining observable quantities

def EvaluateConversionOfData(ValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, ValuesFornHEPhotonsSpectralInxga):
    ValuesForEnergy = ValuesForxgammanHEPhotons*me*c**2 # Conversion from dimensionless energy x_gamma to energy, here in units of J.
    ValuesForEnergyGeV = ValuesForEnergy/(10**9*e) # Conversion of energy from J to GeV.
    # In the following "flux" is short for "spectral flux" and denotes a number of photons that streams through a surface per unit time, per unit area and per unit energy-interval.
    ValuesForDirectedFluxEarth = ValuesFornHEPhotonsSumSpectralInxga*FactornHEPsToFlux # This is the flux and has the unit 1/(J*s*m^2). Directed means, that it is assumed, that that the entire radiation is considered.
    ValuesForDirectedFluxEarthPerGeVcm2s = ValuesForDirectedFluxEarth*e*10**9/(10**4) # Conversion of the flux from 1/(J*s*m^2) to 1/(GeV*s*cm^2).
    ValuesForObsFluxGeVPercm2s = ValuesForEnergyGeV**2*ValuesForDirectedFluxEarthPerGeVcm2s # Conversion to the observational quantity which is in units of GeV/(cm^2*s).
    InterpolatedObjectForValuesForObsFluxGeVPercm2s = interp1d(ValuesForEnergyGeV, ValuesForObsFluxGeVPercm2s, kind='linear', bounds_error=False, fill_value=0.0) # This is an interpolation of the theoretical values.
    ValuesForDirectedTransmittedFluxEarth = ValuesFornHEPhotonsSpectralInxga*FactornHEPsToFlux # This is the transmitted flux and has the unit 1/(J*s*m^2). Transmitted means those HE photons that have been directly injected and traversed the interaction region without having had an interaction.
    ValuesForDirectedTransmittedFluxEarthPerGeVcm2s = ValuesForDirectedTransmittedFluxEarth*e*10**9/(10**4) # Conversion of the transmitted flux from 1/(J*s*m^2) to 1/(GeV*s*cm^2).
    ValuesForTransmittedObsFluxGeVPercm2s = ValuesForEnergyGeV**2*ValuesForDirectedTransmittedFluxEarthPerGeVcm2s # Conversion to the observational quantity which is in units of GeV/(cm^2*s).
    return ValuesForEnergyGeV, ValuesForObsFluxGeVPercm2s, InterpolatedObjectForValuesForObsFluxGeVPercm2s, ValuesForTransmittedObsFluxGeVPercm2s

def EvaluateExportEnergyFluxDensity(ValuesForEnergyGeV,ValuesForObsFluxGeVPercm2s,ValuesForTransmittedObsFluxGeVPercm2s):
    '''The final data points are saved into a .dat-file.'''
    # First, the total emitted flux.
    Outputfile = open("Run %s, Energy-flux-density versus photon energy from %s (total flux).dat" % (RunIdentifier,GetCurrentDate()), 'w')
    Outputfile.write('Values for the energy of the high-energetic photons in GeV, Corresponding values for the flux (nu*F_nu) in GeV/(s*cm^2)\n')
    for i in range(len(ValuesForEnergyGeV)):
        Outputfile.write('%s %s\n' % (ValuesForEnergyGeV[i],ValuesForObsFluxGeVPercm2s[i]))
    Outputfile.close()
    # Second, the flux, that was transmitted without escaping or getting absorbed.
    Outputfile2 = open("Run %s, Energy-flux-density versus photon energy from %s (transmitted flux).dat" % (RunIdentifier,GetCurrentDate()), 'w')
    Outputfile2.write('Values for the energy of the high-energetic photons in GeV, Corresponding values for the transmitted flux (nu*F_nu) in GeV/(s*cm^2)\n')
    for i in range(len(ValuesForEnergyGeV)):
        Outputfile2.write('%s %s\n' % (ValuesForEnergyGeV[i],ValuesForTransmittedObsFluxGeVPercm2s[i]))
    Outputfile2.close()


## Comparison to the Mrk 501 data

def ImportDataPoints(NameOfFileToImport):
    '''The string of a data-file is imported here to extract the frequency and the flux. The name of the file, e. g. 'Mrk501_MWL_MJD_56858.98.txt', has to be inserted.'''
    OpenedLines = ImportFileToString(NameOfFileToImport)
    DataFrequencyHz = np.asarray([]) # An array containing the frequency values from the observational data. This is the frequency in units of Hz.
    DataObsFluxergPercm2s = np.asarray([]) # An array containing the flux- (proportional to nu*F_nu) values of the observational data. This is in units of erg*cm^-2*s^-1.
    DataObsFluxErrorergPercm2s = np.asarray([]) # An array containing the error of the flux-values of the observational data. This is in units of erg*cm^-2*s^-1.
    for Line in OpenedLines:
        if '#' not in Line:
            DataQuarted = Line.replace("\n", "").split() # A DataQuarted contains the number of the data-point, the frequency, the corresponding nu*F_nu value and two additional items.
            DataFrequencyHz = np.append(DataFrequencyHz,eval(DataQuarted[0]))
            DataObsFluxergPercm2s = np.append(DataObsFluxergPercm2s,eval(DataQuarted[1]))
            DataObsFluxErrorergPercm2s = np.append(DataObsFluxErrorergPercm2s,eval(DataQuarted[2]))
    pl.figure(figsize=(12, 9), num="Imported observational data")
    pl.loglog(DataFrequencyHz, DataObsFluxergPercm2s, marker='o', linestyle='None')
    return DataFrequencyHz, DataObsFluxergPercm2s, DataObsFluxErrorergPercm2s

def EvaluateMrk501_35Bins(ValuesForEnergyGeV, ValuesForObsFluxGeVPercm2s, InterpolatedObjectForValuesForObsFluxGeVPercm2s):
    # External observational data for the major flaring day MJD 56857.98 with 35 bins:
    DataFrequencyHz_35Bins, DataObsFluxergPercm2s_35Bins, DataObsFluxErrorergPercm2s_35Bins = ImportDataPoints('Mrk 501 observational data.txt') # The following is imported and assigned: Raw data for the frequency in units of Hz of Mrk501 from MJD 56857.98 with 35 bins, Raw data for the flux in units of erg*cm^-2*s^-1 of Mrk501 from MJD 56857.98 with 35 bins, Raw data for the flux-error in units of erg*cm^-2*s^-1 of Mrk501 from MJD 56857.98 with 35 bins
    
    # Conversion to the plotted units:
    DataEnergyGeV_35Bins = DataFrequencyHz_35Bins*h/(10**9*e) # Data for the energy in units of GeV.
    DataObsFluxGeVPercm2s_35Bins = DataObsFluxergPercm2s_35Bins/(10**16*e) # Data for the flux in units of GeV*cm^-2*s^-1.
    DataObsFluxErrorGeVPercm2s_35Bins = DataObsFluxErrorergPercm2s_35Bins/(10**16*e) # Data for the error of the flux in units of GeV*cm^-2*s^-1.
    
    # Plot:
    pl.figure(figsize=(16, 12), num="Mrk501 with 35 bins: Energy-flux versus photon energy")
    pl.loglog(ValuesForEnergyGeV,ValuesForObsFluxGeVPercm2s, label='Theory (electr. cascade)')
    pl.errorbar(DataEnergyGeV_35Bins,DataObsFluxGeVPercm2s_35Bins, yerr=DataObsFluxErrorGeVPercm2s_35Bins, label='Mrk 501, MJD 56857.98, 35 bins', marker='.', linestyle='None')
    pl.legend(loc="best", fontsize=20)
    ax = pl.gca()
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    pl.xlim(3*10**(-1.0),10**4.0)
    pl.ylim(1*10**(-9),2.0*10**(-7))
    pl.xlabel("High-energy photon energy $\epsilon_\gamma$ in GeV", family='serif', fontsize=26)
    pl.ylabel("Energy flux density $\epsilon_\gamma^2 F_{\epsilon_\gamma}$ in $\mathrm{GeV}/(\mathrm{cm}^2\mathrm{s})$", family='serif', fontsize=26)
    #pl.savefig("Run %s, Energy-flux-density versus photon energy from %s (35 bins).svg" % (RunIdentifier,GetCurrentDate()))
    
    # Chi-square:
    DataEnergyGeV_35Bins_bump = DataEnergyGeV_35Bins[-6:-1] # Only the bump data points.
    DataObsFluxGeVPercm2s_35Bins_bump = DataObsFluxGeVPercm2s_35Bins[-6:-1] # Only the bump data points.
    DataObsFluxErrorGeVPercm2s_35Bins_bump = DataObsFluxErrorGeVPercm2s_35Bins[-6:-1] # Only the bump data points.
    DataObsFluxGeVPercm2s_35Bins_bump_Modelled = InterpolatedObjectForValuesForObsFluxGeVPercm2s(DataEnergyGeV_35Bins_bump) # The theoretical values that correspond to the observed bump data points.
    ChiSquare_35Bins_bump = np.sum(((DataObsFluxGeVPercm2s_35Bins_bump_Modelled-DataObsFluxGeVPercm2s_35Bins_bump)/DataObsFluxErrorGeVPercm2s_35Bins_bump)**2.0)/(len(DataObsFluxGeVPercm2s_35Bins_bump_Modelled)-1.0)
    print('\nConsider the bump data points of the 35 binned data:')
    print("chi-square:", ChiSquare_35Bins_bump)
    
    return DataEnergyGeV_35Bins, DataObsFluxGeVPercm2s_35Bins, DataObsFluxErrorGeVPercm2s_35Bins

def EvaluateImportSSCModel():
    '''The string of the data-file from Pepa's modelling is imported here to extract the frequency and the flux. The name of the file has to be inserted.'''
    OpenedLines = ImportFileToString('SSC model for Mrk 501.dat')
    SSCFrequencyHz = np.asarray([]) # An array containing the frequency values from the SSC-fit by Pepa to the observational data. This is the frequency in units of Hz.
    SSCFluxergPercm2s = np.asarray([]) # An array containing the energy-flux values from the SSC-fit by Pepa to the observational data. This is the energy-flux in units of erg*cm^-2*s^-1..
    for Line in OpenedLines:
        DataPair = Line.replace("\n", "").split() # A DataPair contains the exponent of the frequency and the corresponding exponent of the flux.
        SSCFrequencyHz = np.append(SSCFrequencyHz,10**(eval(DataPair[0])))
        SSCFluxergPercm2s = np.append(SSCFluxergPercm2s,10**(eval(DataPair[1])))
    pl.figure(figsize=(12, 9), num="Imported SSC Data")
    pl.loglog(SSCFrequencyHz, SSCFluxergPercm2s)
    return SSCFrequencyHz, SSCFluxergPercm2s

def EvaluateMrk501_35Bins_IncludeSSC(CurrentNElectrons, DataEnergyGeV_35Bins, DataObsFluxGeVPercm2s_35Bins, DataObsFluxErrorGeVPercm2s_35Bins, ValuesForEnergyGeV, ValuesForObsFluxGeVPercm2s, SSCFrequencyHz, SSCFluxergPercm2s, InterpolatedObjectForValuesForObsFluxGeVPercm2s):
    # External SSC-data:
    SSCFrequencyHz_35Bins = SSCFrequencyHz # These are the SSC-fits by Pepa to the observational data. This is the frequency in units of Hz.
    SSCObsFluxergPercm2s_35Bins = SSCFluxergPercm2s # These are the SSC-fits by Pepa to the observational data. This is the energy-flux in units of erg*cm^-2*s^-1.
    
    # Conversion to the plotted units:
    SSCEnergyGeV_35Bins = SSCFrequencyHz_35Bins*h/(10.0**9*e) # Conversion from frequency to energy in units of GeV.
    SSCObsFluxGeVPercm2s_35Bins = SSCObsFluxergPercm2s_35Bins/(10.0**16*e) # Conversion into units of GeV*cm^-2*s^-1.
    
    # Determining my model at the energy points of the SSC-data:
    ValuesForObsFluxGeVPercm2s_35Bins_SSCEnergies = np.asarray(InterpolatedObjectForValuesForObsFluxGeVPercm2s(SSCEnergyGeV_35Bins)) # Evaluation of the interpolation at those points where the SSC-model is sampled.
    
    # Interpolating the SSC-fit:
    InterpolatedObjectForSSCObsFluxGeVPercm2s_35Bins = interp1d(SSCEnergyGeV_35Bins, SSCObsFluxGeVPercm2s_35Bins, kind='linear', bounds_error=False, fill_value=0.0) # This is an interpolation of the SSC-fit.
    SSCObsFluxGeVPercm2s_35Bins_ValuesForEnergyGeV = np.asarray(InterpolatedObjectForSSCObsFluxGeVPercm2s_35Bins(ValuesForEnergyGeV)) # Sample the SSC-fit at each single point of ValuesForEnergyGeV.
    
    # Addition of SSC-data with cascade-data:
    SumObsFluxGeVPercm2s_35Bins = SSCObsFluxGeVPercm2s_35Bins_ValuesForEnergyGeV+ValuesForObsFluxGeVPercm2s # This is the addition of the SSC-fit with this cascade model evaluated at ValuesForEnergyGeV.
    InterpolatedObjectForMrk501SumObsFluxGeVPercm2s_35Bins = interp1d(ValuesForEnergyGeV, SumObsFluxGeVPercm2s_35Bins, kind='linear', bounds_error=False, fill_value=0.0) # This is an interpolation of the complete theoretical model.
    
    # Plot:
    pl.figure(figsize=(16, 12), num="Mrk501 with 35 bins: Energy-flux-density versus photon-energy, including SSC")
    pl.errorbar(DataEnergyGeV_35Bins,DataObsFluxGeVPercm2s_35Bins, yerr=DataObsFluxErrorGeVPercm2s_35Bins, label='Mrk 501, MJD 56857.98, 35 bins', marker='.', linestyle='None')
    pl.loglog(ValuesForEnergyGeV,ValuesForObsFluxGeVPercm2s, label='Theory (electr. cascade)')
    pl.loglog(SSCEnergyGeV_35Bins,SSCObsFluxGeVPercm2s_35Bins, label='Theory (SSC)')
    pl.loglog(ValuesForEnergyGeV,SumObsFluxGeVPercm2s_35Bins, label='Theory (sum)')
    pl.rc('font',family='serif')
    pl.legend(loc="best", fontsize=22)
    ax = pl.gca()
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    ax.xaxis.set_tick_params(labelsize=26, direction='inout', width=3, length=10, pad=7)
    ax.yaxis.set_tick_params(labelsize=26, direction='inout', width=3, length=10)
    ax.xaxis.set_tick_params(which='minor', labelsize=20, direction='inout', width=2, length=6, pad=7)
    ax.yaxis.set_tick_params(which='minor', labelsize=20, direction='inout', width=2, length=6)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    pl.xlim(3*10**(-1.0),10**4.0)
    pl.ylim(1*10**(-9),2.0*10**(-7))
    pl.xlabel("High-energy photon energy $\epsilon_\gamma$ in GeV", family='serif', fontsize=30)
    pl.ylabel("Energy flux density $\epsilon_\gamma^2 F_{\epsilon_\gamma}$ in $\mathrm{GeV}/(\mathrm{cm}^2\mathrm{s})$", family="serif", fontsize=32)
    pl.savefig("Run %s, Energy-flux-density versus photon energy from %s (35 bins, SSC and cascade).svg" % (RunIdentifier,GetCurrentDate()))
    
    # Chi-square:
    DataEnergyGeV_35Bins_bump = DataEnergyGeV_35Bins[-6:-1] # Only the bump data points.
    DataObsFluxGeVPercm2s_35Bins_bump = DataObsFluxGeVPercm2s_35Bins[-6:-1] # Only the bump data points.
    DataObsFluxErrorGeVPercm2s_35Bins_bump = DataObsFluxErrorGeVPercm2s_35Bins[-6:-1] # Only the bump data points.
    SumObsFluxGeVPercm2s_35Bins_bump_Modelled = np.asarray(InterpolatedObjectForMrk501SumObsFluxGeVPercm2s_35Bins(DataEnergyGeV_35Bins_bump)) # The theoretical values (sum) that correspond to the observed points.
    ChiSquare_35Bins = np.sum(((SumObsFluxGeVPercm2s_35Bins_bump_Modelled-DataObsFluxGeVPercm2s_35Bins_bump)/DataObsFluxErrorGeVPercm2s_35Bins_bump)**2.0)/(len(SumObsFluxGeVPercm2s_35Bins_bump_Modelled)-1.0)
    print('\nConsider the bump data points of the 35 binned data:')
    print("chi-square:", ChiSquare_35Bins)


## Comparison to the 3C 279 data

def Evaluate3C279(ValuesForEnergyGeV,ValuesForObsFluxGeVPercm2s,ValuesForTransmittedObsFluxGeVPercm2s):
    # Read and plot observational data:
    MJD58129FileName = '3C 279 observational data.txt'
    MJD58129FilePath = SearchAFile0rDirectoryEverywhere(MJD58129FileName,False)
    MJD58129Data = np.genfromtxt(os.path.join(MJD58129FilePath,MJD58129FileName)) # This reads all the data (usual points as well as upper limits) in the file into an array.
    MJD58129EnergyGeV = MJD58129Data[:,1]*(1.0+Redshift3C279) # Select the second column, which is the energy in units of GeV and in the observer's frame. The multiplication with (1.0+Redshift3C279) corrects for cosmological redshift and, thus converts into the source's frame.
    MJD58129FluxGeVPercm2s = MJD58129Data[:,4]*10**(-16)/e # This selects the flux-density in units of erg/(cm*s) and converts it to GeV/(cm*s).
    MJD58129FluxErrorGeVPercm2s = MJD58129Data[:,5]*10**(-16)/e # This selects the error (upper and lower error) of the flux-density in units of erg/(cm*s) and converts it to GeV/(cm*s).
    MJD58129EnergyLowerErrorGeV = np.abs(MJD58129Data[:,6]*(1.0+Redshift3C279)-MJD58129EnergyGeV) # This selects the lower fringe of the lower energy error bar and determines the lower error of the energy in units of GeV.
    MJD58129EnergyUpperErrorGeV = np.abs(MJD58129Data[:,7]*(1.0+Redshift3C279)-MJD58129EnergyGeV) # This selects the upper fringe of the upper energy error bar and determines the upper error of the energy in units of GeV.    
    pl.figure(figsize=(16, 12), num="3C279: Energy-flux-density versus photon energy")
    #fig, ax = pl.subplots(figsize=(14, 10))
    pl.rc('font', family='serif')
    pl.errorbar(MJD58129EnergyGeV[:-2], MJD58129FluxGeVPercm2s[:-2], xerr=[(MJD58129EnergyLowerErrorGeV[:-2]),(MJD58129EnergyUpperErrorGeV[:-2])], yerr=MJD58129FluxErrorGeVPercm2s[:-2], marker=".", color="g", elinewidth=2.5, linestyle='None', capsize=0.0, label="MJD 58129-58150", zorder=2) # First, plot the usual points with errors.
    pl.errorbar(MJD58129EnergyGeV[-2:], MJD58129FluxGeVPercm2s[-2:], xerr=[(MJD58129EnergyLowerErrorGeV[-2:]),(MJD58129EnergyUpperErrorGeV[-2:])], yerr=(0.2*MJD58129FluxGeVPercm2s[-2:]), uplims=np.array([1, 1], dtype=bool), elinewidth=2.5, color="g", linestyle='None', capsize=0.0, zorder=2) # Second, plot the upper limits.
    
    # Read and plot a theoretical fit (seems obsolete):
    # FitFileName = 'Energy-flux versus photon energy from 2018-06-06 22-51 (total flux).dat'
    # FitFilePath = SearchAFile0rDirectoryEverywhere(FitFileName,False)
    # FitData = np.genfromtxt(os.path.join(FitFilePath,FitFileName), skip_header=1) # This reads all the data (usual points as well as upper limits) in the file into an array.
    # FitEnergyGeV = FitData[:,0]
    # FitFluxGeVPercm2s = FitData[:,1]
    # pl.plot(FitEnergyGeV, FitFluxGeVPercm2s, label='Fit %s' % FitFileName[38:-4])
    
    # Plot theoretical fit of the current run:
    #pl.plot(ValuesForEnergyGeV, ValuesForObsFluxGeVPercm2s, label='Total energy flux-density')
    #pl.plot(ValuesForEnergyGeV, ValuesForTransmittedObsFluxGeVPercm2s, label='Injected energy-flux-density, diminished by pair-absorption')
    
    # Import and plot old fits
    Theo1=np.genfromtxt(os.path.join(SearchAFile0rDirectoryEverywhere("3C 279 modelling case 1a.dat",False),"3C 279 modelling case 1a.dat"), skip_header=1) # Import theoretical scenario 1, also called case 1a in the thesis (run 188). Electron + photon injection.
    Theo2=np.genfromtxt(os.path.join(SearchAFile0rDirectoryEverywhere("3C 279 modelling case 1b.dat",False),"3C 279 modelling case 1b.dat"), skip_header=1) # Import theoretical scenario 2, also called case 1b in the thesis (run 189). Only photon injection.
    Theo3=np.genfromtxt(os.path.join(SearchAFile0rDirectoryEverywhere("3C 279 modelling case 1c.dat",False),"3C 279 modelling case 1c.dat"), skip_header=1) # Import theoretical scenario 3, also called case 1c in the thesis (run 193). Only electron injection.
    Theo4=np.genfromtxt(os.path.join(SearchAFile0rDirectoryEverywhere("3C 279 modelling case 2'.dat",False),"3C 279 modelling case 2'.dat"), skip_header=1) # Import theoretical scenario 1, also called case 2' in the thesis (run 155). Has reduced ExPs.
    Theo5=np.genfromtxt(os.path.join(SearchAFile0rDirectoryEverywhere("3C 279 modelling case 2''.dat",False),"3C 279 modelling case 2''.dat"), skip_header=1) # Import theoretical scenario 1, also called case 2'' in the thesis (run 158). Has reduced ExPs and compensation in HEP normalisation.
    #Theo6=np.genfromtxt(os.path.join(SearchAFile0rDirectoryEverywhere("Energy-flux-density versus photon energy from 2020-03-29 14-02 (total flux).dat",False),"Energy-flux-density versus photon energy from 2020-03-29 14-02 (total flux).dat"), skip_header=1) # Import theoretical scenario 1, (run 160). Has reduced ExPs and compensation in HEP normalisation and electron injection.
    #Theo7=np.genfromtxt(os.path.join(SearchAFile0rDirectoryEverywhere("Energy-flux-density versus photon energy from 2021-03-18 13-14 (total flux).dat",False),"Energy-flux-density versus photon energy from 2021-03-18 13-14 (total flux).dat"), skip_header=1) # Import theoretical scenario 1, also called case 2'' in the thesis (run 195, equivalent to run 158). Has reduced ExPs and compensation in HEP normalisation and electron injection.
    Theo8=np.genfromtxt(os.path.join(SearchAFile0rDirectoryEverywhere("3C 279 modelling case 2'''.dat",False),"3C 279 modelling case 2'''.dat"), skip_header=1) # Import theoretical scenario 1, also called case 2''' in the thesis (run 196). Has reduced ExPs and compensation in HEP normalisation and electron injection.
    Theo9=np.genfromtxt(os.path.join(SearchAFile0rDirectoryEverywhere("3C 279 modelling case 2''''.dat",False),"3C 279 modelling case 2''''.dat"), skip_header=1) # Import theoretical scenario 1, also called case 2'''' in the thesis (run 198). Has reduced ExPs and compensation in HEP normalisation and electron injection.
    Theo1x=Theo1[:,0]
    Theo2x=Theo2[:,0]
    Theo3x=Theo3[:,0]
    Theo4x=Theo4[:,0]
    Theo5x=Theo5[:,0]
    #Theo6x=Theo6[:,0]
    #Theo7x=Theo7[:,0]
    Theo8x=Theo8[:,0]
    Theo9x=Theo9[:,0]
    Theo1y=Theo1[:,1]
    Theo2y=Theo2[:,1]
    Theo3y=Theo3[:,1]
    Theo4y=Theo4[:,1]
    Theo5y=Theo5[:,1]
    #Theo6y=Theo6[:,1]
    #Theo7y=Theo7[:,1]
    Theo8y=Theo8[:,1]
    Theo9y=Theo9[:,1]
    pl.plot(Theo1x, Theo1y, label='Case 1a and 1c', color="b", linewidth=3.2, zorder=1)
    pl.plot(Theo2x, Theo2y, label='Case 1b', color="magenta", linewidth=3.2, zorder=1)
    #pl.plot(Theo3x, Theo3y, label='Case 1c', color="c", linewidth=2.5, zorder=1)
    pl.plot(Theo4x, Theo4y, label="Case 2'", color='#444444', linewidth=2.9, linestyle=':', zorder=1)
    pl.plot(Theo5x, Theo5y, label="Case 2''", color='#333333', linewidth=2.9, linestyle='-.', zorder=1)
    #pl.plot(Theo6x, Theo6y, label="Run160", color='#111111', linewidth=2.8, linestyle=':', zorder=1)
    #pl.plot(Theo7x, Theo7y, label="Run195", color='#111111', linewidth=2.8, linestyle='-', zorder=1)
    pl.plot(Theo8x, Theo8y, label="Case 2'''", color='#111111', linewidth=2.8, linestyle='--', zorder=1)
    pl.plot(Theo9x, Theo9y, label="Case 2''''", color='#111111', linewidth=2.3, linestyle='-', zorder=1)
    
    ax = pl.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_tick_params(labelsize=34, direction='inout', width=5, length=20, pad=17)
    ax.yaxis.set_tick_params(labelsize=34, direction='inout', width=5, length=20, pad=10)
    ax.xaxis.set_tick_params(which='minor', labelsize=26, direction='inout', width=3, length=10, pad=7)
    ax.yaxis.set_tick_params(which='minor', labelsize=26, direction='inout', width=3, length=10)
    ax.spines['left'].set_linewidth(4)
    ax.spines['right'].set_linewidth(4)
    ax.spines['top'].set_linewidth(4)
    ax.spines['bottom'].set_linewidth(4)
    pl.xlabel('$x_{\gamma} \, m_{\mathrm{e}} c^2 / \, \mathrm{GeV}$', fontsize=36)
    pl.ylabel('$(x_{\gamma} \, m_{\mathrm{e}} c^2)^2 \, F_{\mathrm{casc}} \, / \, \mathrm{GeV} \, \mathrm{cm}^{-2} \mathrm{s}^{-1}$', fontsize=36)
    pl.axis([0.08,600,2.7*10**(-8),1.2*10**(-6)])#pl.axis([0.02,800,1.0*10**(-8),2.0*10**(-6)])
    pl.legend(loc='lower left', fontsize=29)
    pl.subplots_adjust(top=0.97, bottom=0.16, left=0.14, right=0.98)
    pl.savefig("Run %s, Energy-flux-density versus photon energy from %s.pdf" % (RunIdentifier,GetCurrentDate()))
    pl.savefig("Run %s, Energy-flux-density versus photon energy from %s.svg" % (RunIdentifier,GetCurrentDate()))
    
    # Chi-square of the current fit:
    # InterpolatedObjectCurrent = interp1d(ValuesForEnergyGeV, ValuesForObsFluxGeVPercm2s, kind='linear', bounds_error=False, fill_value=0.0) # Interpolation of the current fit. Energy in GeV is the independent variable and energy-flux-density in units of GeV*cm^{-2}*s^{-1} is returned.
    # CurrentAtDataPoints = InterpolatedObjectCurrent(MJD58129EnergyGeV[:-2]) # The values that correspond to the energies of the observed data points, without the upper limits.
    # CurrentChiSquare = np.sum(((CurrentAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0) # Determine chi-square for the current fit. The upper limits are excluded.
    # CurrentChiSquareRed = CurrentChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-7) # Determine chi-square for the current fit. The upper limits are excluded.
    # print("Current fit: chi-square =", CurrentChiSquare)
    # print("Current fit: Reduced chi-square =", CurrentChiSquareRed)
    
    # Chi-square of old fits:
    InterpolatedObjectTheo1 = interp1d(Theo1x, Theo1y, kind='linear', bounds_error=False, fill_value=0.0) # Interpolation of the theoretical model case 1a. Energy in GeV is the independent variable and energy-flux-density in units of GeV*cm^{-2}*s^{-1} is returned.
    InterpolatedObjectTheo2 = interp1d(Theo2x, Theo2y, kind='linear', bounds_error=False, fill_value=0.0) # Same for model case 1b.
    InterpolatedObjectTheo3 = interp1d(Theo3x, Theo3y, kind='linear', bounds_error=False, fill_value=0.0) # Same for model case 1c.
    InterpolatedObjectTheo4 = interp1d(Theo4x, Theo4y, kind='linear', bounds_error=False, fill_value=0.0) # Same for model case 2'.
    InterpolatedObjectTheo5 = interp1d(Theo5x, Theo5y, kind='linear', bounds_error=False, fill_value=0.0) # Same for model case 2''.
    #InterpolatedObjectTheo6 = interp1d(Theo6x, Theo6y, kind='linear', bounds_error=False, fill_value=0.0) # Same for model 6.
    #InterpolatedObjectTheo7 = interp1d(Theo7x, Theo7y, kind='linear', bounds_error=False, fill_value=0.0) # Same for model case 2''.
    InterpolatedObjectTheo8 = interp1d(Theo8x, Theo8y, kind='linear', bounds_error=False, fill_value=0.0) # Same for model case 2'''.
    InterpolatedObjectTheo9 = interp1d(Theo9x, Theo9y, kind='linear', bounds_error=False, fill_value=0.0) # Same for model case 2''''.
    Theo1yAtDataPoints = InterpolatedObjectTheo1(MJD58129EnergyGeV[:-2]) # The theoretical values that correspond to the energies of the observed data points, without the upper limits.
    Theo2yAtDataPoints = InterpolatedObjectTheo2(MJD58129EnergyGeV[:-2])
    Theo3yAtDataPoints = InterpolatedObjectTheo3(MJD58129EnergyGeV[:-2])
    Theo4yAtDataPoints = InterpolatedObjectTheo4(MJD58129EnergyGeV[:-2])
    Theo5yAtDataPoints = InterpolatedObjectTheo5(MJD58129EnergyGeV[:-2])
    #Theo6yAtDataPoints = InterpolatedObjectTheo6(MJD58129EnergyGeV[:-2])
    #Theo7yAtDataPoints = InterpolatedObjectTheo7(MJD58129EnergyGeV[:-2])
    Theo8yAtDataPoints = InterpolatedObjectTheo8(MJD58129EnergyGeV[:-2])
    Theo9yAtDataPoints = InterpolatedObjectTheo9(MJD58129EnergyGeV[:-2])
    Sce1ChiSquare = np.sum(((Theo1yAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0) # Determine chi-square. The upper limits are excluded.
    Sce2ChiSquare = np.sum(((Theo2yAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0)
    Sce3ChiSquare = np.sum(((Theo3yAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0)
    Sce4ChiSquare = np.sum(((Theo4yAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0)
    Sce5ChiSquare = np.sum(((Theo5yAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0)
    #Sce6ChiSquare = np.sum(((Theo6yAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0)
    #Sce7ChiSquare = np.sum(((Theo7yAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0)
    Sce8ChiSquare = np.sum(((Theo8yAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0)
    Sce9ChiSquare = np.sum(((Theo9yAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0)
    Sce1ChiSquareRed = Sce1ChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-7) # Determine reduced chi-square. The upper limits are excluded.
    Sce2ChiSquareRed = Sce2ChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-5)
    Sce3ChiSquareRed = Sce3ChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-7)
    Sce4ChiSquareRed = Sce4ChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-7)
    Sce5ChiSquareRed = Sce5ChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-7)
    #Sce6ChiSquareRed = Sce6ChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-7)
    #Sce7ChiSquareRed = Sce7ChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-7)
    Sce8ChiSquareRed = Sce8ChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-7)
    Sce9ChiSquareRed = Sce9ChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-7)
    print("Case 1a:    chi-square =", Sce1ChiSquare)
    print("Case 1b:    chi-square =", Sce2ChiSquare)
    print("Case 1c:    chi-square =", Sce3ChiSquare)
    print("Case 2':    chi-square =", Sce4ChiSquare)
    print("Case 2'':   chi-square =", Sce5ChiSquare)
    #print("Case Theo6: chi-square =", Sce6ChiSquare)
    #print("Case 2'':   chi-square =", Sce7ChiSquare)
    print("Case 2''':  chi-square =", Sce8ChiSquare)
    print("Case 2'''': chi-square =", Sce9ChiSquare)
    print("Case 1a:    Reduced chi-square =", Sce1ChiSquareRed)
    print("Case 1b:    Reduced chi-square =", Sce2ChiSquareRed)
    print("Case 1c:    Reduced chi-square =", Sce3ChiSquareRed)
    print("Case 2':    Reduced chi-square =", Sce4ChiSquareRed)
    print("Case 2'':   Reduced chi-square =", Sce5ChiSquareRed)
    #print("Case Theo6: Reduced chi-square =", Sce6ChiSquareRed)
    #print("Case 2'':   Reduced chi-square =", Sce7ChiSquareRed)
    print("Case 2''':  Reduced chi-square =", Sce8ChiSquareRed)
    print("Case 2'''': Reduced chi-square =", Sce9ChiSquareRed)
    
    # Conventional fits with chi-square:
    def FitFunctionLogparabolic(xga,P): # This is a logparabola and it can be used as the fit function.
        Norm=P[0]
        Slope=P[1]
        Curvature=P[2]
        return Norm*xga**(Slope+Curvature*np.log10(xga))
    
    def FitFunctionPLExpCutOff(xga,P): # A power-law with an exponential cut-off.
        Norm=P[0]
        Slope=P[1]
        CutOff=P[2]
        return Norm*xga**Slope*np.exp(-xga/CutOff)
    
    def ResidualsFunctionForFit(FitParameters, xgaData, yData, FitFunction): # The residual function, whose sum of squares should be minimized. It takes as arguments a vector FitParameters, which contains the fitting parameters, the gamma-data (XData) and the flux-density-data (YData). It returns the residual vector.
        ErrorVector = yData - FitFunction(xgaData,FitParameters)
        return ErrorVector
    
    StartingParametersFitLogparabolic = [10**(-6), 0, 0] # Collect the starting positions for the fitting parameters in a list.
    StartingParametersFitPLExpCutOff = [10**(-6), 0, 100] # Collect the starting positions for the fitting parameters in a list.
    
    ParametersFitLogparabolic = optimize.leastsq(ResidualsFunctionForFit, StartingParametersFitLogparabolic, args=(MJD58129EnergyGeV[:-2], MJD58129FluxGeVPercm2s[:-2], FitFunctionLogparabolic))[0]
    print("The parameters found via a least-squares fitting algorithm for a logparabola are:", ParametersFitLogparabolic)
    ParametersFitPLExpCutOff = optimize.leastsq(ResidualsFunctionForFit, StartingParametersFitPLExpCutOff, args=(MJD58129EnergyGeV[:-2], MJD58129FluxGeVPercm2s[:-2], FitFunctionPLExpCutOff))[0]
    print("The parameters found via a least-squares fitting algorithm for a PL + exp cut-off are:", ParametersFitPLExpCutOff)
    
    # Plot the fitted function together with the data points:
    pl.figure(figsize=(15, 12), num="Flux-density versus photon energy with fitted function")
    pl.errorbar(MJD58129EnergyGeV[:-2], MJD58129FluxGeVPercm2s[:-2], xerr=[(MJD58129EnergyLowerErrorGeV[:-2]),(MJD58129EnergyUpperErrorGeV[:-2])], yerr=MJD58129FluxErrorGeVPercm2s[:-2], marker=".", color="b", elinewidth=3.5, linestyle='None', capsize=0.0, label="MJD 58129-58150") # First, plot the usual points with errors.
    PlotRange = np.logspace(np.log10(0.1),np.log10(300),100)
    pl.loglog(PlotRange, FitFunctionLogparabolic(PlotRange,ParametersFitLogparabolic), label='Logparabolic fit', color='r')
    pl.loglog(PlotRange, FitFunctionPLExpCutOff(PlotRange,ParametersFitPLExpCutOff), label='PL + exp. cut-off fit', color='g')
    ax = pl.gca()
    ax = pl.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_tick_params(labelsize=40, direction='inout', width=5, length=20, pad=17)
    ax.yaxis.set_tick_params(labelsize=40, direction='inout', width=5, length=20, pad=10)
    ax.xaxis.set_tick_params(which='minor', labelsize=32, direction='inout', width=3, length=10, pad=7)
    ax.yaxis.set_tick_params(which='minor', labelsize=32, direction='inout', width=3, length=10)
    ax.spines['left'].set_linewidth(4)
    ax.spines['right'].set_linewidth(4)
    ax.spines['top'].set_linewidth(4)
    ax.spines['bottom'].set_linewidth(4)
    pl.xlabel('$x_{\gamma} \cdot \, m_{\mathrm{e}} c^2 \, \, \, [\mathrm{GeV}]$', fontsize=42)
    pl.ylabel('$F_{\mathrm{casc}}(x_{\gamma}) \, \, \, [\mathrm{GeV} \cdot \, \mathrm{cm}^{-2} \mathrm{s}^{-1}]$', fontsize=42)
    pl.axis([0.08,600,2.7*10**(-8),1.2*10**(-6)])#pl.axis([0.02,800,1.0*10**(-8),2.0*10**(-6)])
    pl.legend(loc='lower left', fontsize=26)
    pl.subplots_adjust(top=0.97, bottom=0.18, left=0.17, right=0.97)
    
    # Chi-square:
    FitFunctionLogparabolicAtDataPoints = FitFunctionLogparabolic(MJD58129EnergyGeV[:-2],ParametersFitLogparabolic) # The values of the fit function that correspond to the energies of the observed data points, without the upper limits.
    FitFunctionPLExpCutOffAtDataPoints = FitFunctionPLExpCutOff(MJD58129EnergyGeV[:-2],ParametersFitPLExpCutOff) # The values of the fit function that correspond to the energies of the observed data points, without the upper limits.
    LogparabolicChiSquare = np.sum(((FitFunctionLogparabolicAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0) # Determine chi-square for the logparabolic fit.
    PLExpCutOffChiSquare = np.sum(((FitFunctionPLExpCutOffAtDataPoints-MJD58129FluxGeVPercm2s[:-2])/MJD58129FluxErrorGeVPercm2s[:-2])**2.0) # Same for the PL + exp cut-off fit.
    LogparabolicChiSquareRed = LogparabolicChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-3) # Determine reduced chi-square for the logparabolic fit.
    PLExpCutOffChiSquareRed = PLExpCutOffChiSquare/(len(MJD58129FluxGeVPercm2s[:-2])-3) # Same for the PL + exp cut-off fit.
    print("Logparabolic fit: Reduced chi-square =", LogparabolicChiSquareRed)
    print("PL + exp. cut-off fit: Reduced chi-square =", PLExpCutOffChiSquareRed)


## Evaluation via a shell:

if __name__ == "__main__" and PerformanceMode == 'CommonIteration':      # Test block

    print('\n\n-------------------------- Main evaluation with common iteration --------------------------')
    time.sleep(3)
    StartingTimeOfEntireEvaluation = time.time() # Measure the spent time of the entire evaluation.
    
    print('\nCalling EvaluateExportInputValues:');EvaluateExportInputValues()
    # print('\nCalling EvaluateImportDataOfIteration:')
    # if UsedIterationScheme == 'IterationStepByStepPointsFromLToR':
    #     CurrentNElectrons, CurrentValuesForgaIt, CurrentLastValuesForgaIt, ValuesForNextNElectrons, NElectronsIterated, CurrentIterationCounter = EvaluateImportDataOfIteration()
    #     LastValuesForNextNElectrons=ValuesForNextNElectrons # This has to be performed if an old result was imported. In this case, LastValuesForNextNElectrons has to be inserted as sixth argument in EvaluateIteration.
    # elif UsedIterationScheme == 'PointsFromRToLIterationPointwise':
    #     CurrentNElectrons, CurrentValuesForgaIt, CurrentLastValuesForgaIt, CurrentValuesForNElectrons, NElectronsIterated = EvaluateImportDataOfIteration('N(gamma) versus gamma from 2020-07-30 11-51 (data).dat')
    #     ValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, ValuesFornHEPhotonsSpectralInxga, nHEPsInterpolated = EvaluateImportDataOfHEPs('nHEP versus energy from 2020-08-06 12-15.dat')
    # print('\nCalling EvaluatePreprocessOldIterationForNewOne:');CurrentIterationCounter, CurrentNElectrons, CurrentLastValuesForgaIt, NElectronsIterated, CurrentValuesForgaIt, ValuesForNextNElectrons = EvaluatePreprocessOldIterationForNewOne(CurrentValuesForgaIt,ValuesForNextNElectrons,NElectronsIterated)
    
    print('\nCalling EvaluateIteration and EvaluatePlotForIteration:')
    if UsedIterationScheme == 'IterationStepByStepPointsFromLToR':
        CurrentIterationCounter, CurrentNElectrons, CurrentLastValuesForgaIt, NElectronsIterated, CurrentValuesForgaIt, ValuesForNextNElectrons = EvaluateIteration(CurrentIterationCounter, CurrentNElectrons, CurrentLastValuesForgaIt, NElectronsIterated, CurrentNElectronsga0)
        EvaluatePlotForIteration(CurrentIterationCounter, NElectronsIterated, CurrentNElectrons)
    elif UsedIterationScheme == 'PointsFromRToLIterationPointwise':
        CurrentNElectrons, CurrentValuesForNElectrons, NElectronsIterated = EvaluateIteration(CurrentNElectrons, CurrentValuesForgaIt, CurrentValuesForNElectrons, NElectronsIterated, None, None)
        EvaluatePlotForIteration(CurrentValuesForgaIt, NElectronsIterated, CurrentNElectrons)
    
    # if UsedIterationScheme == 'IterationStepByStepPointsFromLToR':
    #     print('\nCalling EvaluateNElectronsOfEq8:');CurrentNElectrons, ValuesForgaLogComplete, ValuesForNElectronsComplete = EvaluateNElectronsOfEq8(CurrentNElectrons,CurrentIterationCounter,CurrentLastValuesForgaIt,ValuesForNextNElectrons)
    #     print('\nCalling EvaluateExportDataComplete:');EvaluateExportDataComplete(ValuesForgaLogComplete, ValuesForNElectronsComplete, StorageOfSecondTermOfNumeratorOfEq8)
    #     print('\nCalling EvaluateImportDataComplete:');CurrentNElectrons, ValuesForgaLogComplete, ValuesForNElectronsComplete, StorageOfSecondTermOfNumeratorOfEq8 = EvaluateImportDataComplete()
    #     print('\nCalling EvaluatePlotForNElectronsComplete:');EvaluatePlotForNElectronsComplete(ValuesForgaLogComplete, ValuesForNElectronsComplete)
    #     print('\nCalling EvaluateSmearingOutTheStepViaFitting:');ValuesForgaLogSmoothed, ValuesForNElectronsSmoothed, CurrentNElectrons = EvaluateSmearingOutTheStepViaFitting(ValuesForgaLogComplete, ValuesForNElectronsComplete, CurrentNElectrons, CurrentIterationCounter, CurrentLastValuesForgaIt)
    #     print('\nCalling EvaluatePlotForNElectronsSmoothed:');EvaluatePlotForNElectronsSmoothed(ValuesForgaLogSmoothed, ValuesForNElectronsSmoothed)
    
    if UsedIterationScheme == 'IterationStepByStepPointsFromLToR':
        ValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, ValuesFornHEPhotonsSpectralInxga = EvaluateEquation5AndnHEPhotons(CurrentIterationCounter, CurrentValuesForgaIt, ValuesForNextNElectrons, None, None)
    elif UsedIterationScheme == 'PointsFromRToLIterationPointwise':
        ValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, ValuesFornHEPhotonsSpectralInxga = EvaluateEquation5AndnHEPhotons(1, CurrentValuesForgaIt, CurrentValuesForNElectrons, None, None)
    
    if UsedIterationScheme == 'IterationStepByStepPointsFromLToR':
        print('\nIC-created HE-photon energy-density:', ICCreatedEnergyDensity(Usedn0,CurrentIterationCounter,CurrentLastValuesForgaIt,ValuesForNextNElectrons,None,None,0.1)*me*c**2, "J/m^3")
        print('\nIC-created HE-photon power-density:', ICCreatedPowerDensity(Usedn0,CurrentIterationCounter,HEPhotonsxga1,CurrentLastValuesForgaIt,ValuesForNextNElectrons,None,None,0.1), "1/(m^3*s)")
    elif UsedIterationScheme == 'PointsFromRToLIterationPointwise':
        CurrentElectronsEnergyDensity = ElectronsEnergyDensity(CurrentLastValuesForgaIt,CurrentValuesForNElectrons)*me*c**2
        print("\nElectrons' energy-density:                        ", CurrentElectronsEnergyDensity, "J/m^3")
        print("    Corresponding electrons' luminosity:          ", CurrentElectronsEnergyDensity*FactorEnergyDensityToLuminosity('InteractionRadius')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('InteractionRadius')[0])
        print("    Corresponding electrons' luminosity:          ", CurrentElectronsEnergyDensity*FactorEnergyDensityToLuminosity('RClouds')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('RClouds')[0])
        CurrentICCreatedEnergyDensity = ICCreatedEnergyDensity(Usedn0,1,CurrentLastValuesForgaIt,CurrentValuesForNElectrons,None,None,0.1)*me*c**2
        print('\nIC-created HE-photon energy-density:              ', CurrentICCreatedEnergyDensity, "J/m^3")
        print('    Corresponding IC-created HE-photon luminosity:', CurrentICCreatedEnergyDensity*FactorEnergyDensityToLuminosity('InteractionRadius')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('InteractionRadius')[0])
        print('    Corresponding IC-created HE-photon luminosity:', CurrentICCreatedEnergyDensity*FactorEnergyDensityToLuminosity('RClouds')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('RClouds')[0])
        CurrentICCreatedPowerDensity = ICCreatedPowerDensity(Usedn0,1,10,CurrentLastValuesForgaIt,CurrentValuesForNElectrons,None,None,0.1)
        print('\nIC-created HE-photon power-density:               ', CurrentICCreatedPowerDensity, "1/(m^3*s)")
        print('    Corresponding IC-created HE-photon power:     ', SphericalVolume('InteractionRadius')[1]*CurrentICCreatedPowerDensity*me*c**2, "W (using R=%s)" % SphericalVolume('InteractionRadius')[0])
        print('    Corresponding IC-created HE-photon power:     ', SphericalVolume('RClouds')[1]*CurrentICCreatedPowerDensity*me*c**2, "W (using R=%s)" % SphericalVolume('RClouds')[0])
        
        if IncludeSynchrotron:
            TrueSyncPsEnergyDensity = SyncPsEnergyDensity(1,5,CurrentValuesForgaIt,CurrentValuesForNElectrons,0.1)*me*c**2 # The whole energy-density of synchrotron-radiation in units of J/m^3.
            LostSyncPsEnergyDensity = SyncPsEnergyDensity(0.1/SynchrotronCriticalxSync(CurrentLowestga),5,CurrentValuesForgaIt,CurrentValuesForNElectrons,0.1)*me*c**2 # The energy-density of synchrotron-radiation, that is lost by resetting the upper cut-off, in units of J/m^3.
            print("\nSyncPs' total, true energy-density:", TrueSyncPsEnergyDensity,"J/m^3")
            print('Fraction of energy-density, that is lost by resetting the upper cut-off:', LostSyncPsEnergyDensity/TrueSyncPsEnergyDensity)
                
            TrueSyncPsNumberDensity = SyncPsNumberDensity(10**(-6),5,CurrentValuesForgaIt,CurrentValuesForNElectrons,0.1) # The whole number-density of synchrotron-radiation in units of 1/m^3.
            LostSyncPsNumberDensity = SyncPsNumberDensity(0.1/SynchrotronCriticalxSync(CurrentLowestga),5,CurrentValuesForgaIt,CurrentValuesForNElectrons,0.1) # The number-density of synchrotron-radiation, that is lost by resetting the upper cut-off, in units of 1/m^3.
            print("SyncPs' total, true number-density:", TrueSyncPsNumberDensity,"/m^3")
            print("Fraction of number-density, that is lost by resetting the upper cut-off:", LostSyncPsNumberDensity/TrueSyncPsNumberDensity)
    
    print('\nCalling EvaluateConversionOfData:');ValuesForEnergyGeV, ValuesForObsFluxGeVPercm2s, InterpolatedObjectForValuesForObsFluxGeVPercm2s, ValuesForTransmittedObsFluxGeVPercm2s = EvaluateConversionOfData(ValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, ValuesFornHEPhotonsSpectralInxga)
    print('\nCalling EvaluateExportEnergyFluxDensity:');EvaluateExportEnergyFluxDensity(ValuesForEnergyGeV,ValuesForObsFluxGeVPercm2s,ValuesForTransmittedObsFluxGeVPercm2s)
    
    print('\nCalling EvaluateMrk501_35Bins:');DataEnergyGeV_35Bins, DataObsFluxGeVPercm2s_35Bins, DataObsFluxErrorGeVPercm2s_35Bins = EvaluateMrk501_35Bins(ValuesForEnergyGeV, ValuesForObsFluxGeVPercm2s, InterpolatedObjectForValuesForObsFluxGeVPercm2s)
    print('\nCalling EvaluateImportSSCModel:');SSCFrequencyHz, SSCFluxergPercm2s = EvaluateImportSSCModel()
    print('\nCalling EvaluateMrk501_35Bins_IncludeSSC:');EvaluateMrk501_35Bins_IncludeSSC(CurrentNElectrons, DataEnergyGeV_35Bins, DataObsFluxGeVPercm2s_35Bins, DataObsFluxErrorGeVPercm2s_35Bins, ValuesForEnergyGeV, ValuesForObsFluxGeVPercm2s, SSCFrequencyHz, SSCFluxergPercm2s, InterpolatedObjectForValuesForObsFluxGeVPercm2s)
    
    EvaluateLowerBorderOfNonVanRangeOfCOfA1AndCOfA21AndElectronLossRate(None,None,CurrentValuesForgaIt,1)
    EvaluatePOfB11AndHEPhotonLossRate(None,None)
    
    EndingTimeOfEntireEvaluation = time.time() # Now, determine the point of time after the end of the evaluation.
    TimeIntervalOfEntireEvaluation = (EndingTimeOfEntireEvaluation-StartingTimeOfEntireEvaluation)/3600.0 # Determine the time in hours, that was spent by the evaluation. 
    print("\nTime spent on the entire evaluation:", TimeIntervalOfEntireEvaluation, "hours\n")

if __name__ == "__main__" and PerformanceMode == 'SuperIteration' and UsedIterationScheme == 'PointsFromRToLIterationPointwise':      # Test block

    print('\n\n-------------------------- Main evaluation with super-iteration --------------------------')
    time.sleep(10)
    StartingTimeOfEntireEvaluation = time.time() # Measure the spent time of the entire evaluation.
    
    print('\nCalling EvaluateExportInputValues:');EvaluateExportInputValues()
    
    # print('\nCalling EvaluateImportDataOfIteration:')
    # CurrentNElectrons, CurrentValuesForgaIt, CurrentLastValuesForgaIt, CurrentValuesForNElectrons, NElectronsSuperIt, CurrentNElectronsga0, CurrentLowestga, SuperItCounter, nSyncPs, ValuesForxSyncSuperIt, ValuesFornSyncPs, nSyncPsSuperIt, CurrentxSync0, CurrentxSync1, CurrentxSync0Used, CurrentxTotal0 = EvaluateImportDataOfIteration('N(gamma) versus gamma from 2020-03-24 10-23 (data).dat','nSyncPs versus xSync from 2020-02-19 19-23 (SuperIt=1, data).dat')
    
    print('\nCalling EvaluateSuperIteration and EvaluatePlotForIteration:')
    CurrentNElectrons, CurrentValuesForgaIt, CurrentLastValuesForgaIt, CurrentValuesForNElectrons, NElectronsSuperIt, CurrentNElectronsga0, CurrentLowestga, SuperItCounter, nSyncPs, ValuesForxSyncSuperIt, ValuesFornSyncPs, nSyncPsSuperIt, CurrentxSync0, CurrentxSync1, CurrentxSync0Used, CurrentxTotal0 = EvaluateSuperIteration(CurrentNElectrons, CurrentValuesForgaIt, CurrentLastValuesForgaIt, CurrentValuesForNElectrons, NElectronsSuperIt, CurrentNElectronsga0, CurrentLowestga, SuperItCounter, nSyncPs, ValuesForxSyncSuperIt, ValuesFornSyncPs, nSyncPsSuperIt, CurrentxSync0, CurrentxSync1, CurrentxSync0Used, CurrentxTotal0) # If an imported SuperIt is to be resumed, one has to plug in True as last argument, here.

    NumberOfTasksOfFinalCalculation = 13 # Some additional results (numbers as well as plots are determined now in a multiprocessing pool (called FinalCalculation). This number of tasks has to be equal to the number of .apply_async items, that are appended to FinalCalculationResultQueue in the following.
    NumberOfProcessesOfFinalCalculation=min(NumberOfTasksOfFinalCalculation,int(np.floor(multiprocessing.cpu_count()/2))) # multiprocessing.cpu_count()/2 gives the number of physical CPUs available on this machine. 
    PoolOfFinalCalculation = multiprocessing.Pool(processes=NumberOfProcessesOfFinalCalculation) # Initialise the pool of workers for multiple processes.
    
    FinalCalculationResultQueue = [] # Create the queue of tasks and put the following tasks into it:
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(EvaluateEquation5AndnHEPhotons, args=[1,CurrentValuesForgaIt,CurrentValuesForNElectrons,ValuesForxSyncSuperIt,ValuesFornSyncPs])) # EvaluateEquation5AndnHEPhotons. This will plot and save the HEP-distribution and append a tuple of three arrays to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(ElectronsEnergyDensity, args=[CurrentLastValuesForgaIt,CurrentValuesForNElectrons])) # ElectronsEnergyDensity. This will append a number to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(ICCreatedEnergyDensity, args=[Usedn0,1,CurrentLastValuesForgaIt,CurrentValuesForNElectrons,ValuesForxSyncSuperIt,ValuesFornSyncPs], kwds={'RelativeError': 0.1})) # ICCreatedEnergyDensity. This will append a number to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(ICCreatedPowerDensity, args=[Usedn0,1,10,CurrentLastValuesForgaIt,CurrentValuesForNElectrons,ValuesForxSyncSuperIt,ValuesFornSyncPs], kwds={'RelativeError': 0.1})) # ICCreatedPowerDensity. This will append a number to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(SyncPsEnergyDensity, args=[1,5,CurrentValuesForgaIt,CurrentValuesForNElectrons], kwds={'RelativeError': 0.1})) # SyncPsEnergyDensity. This will append a number to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(SyncPsEnergyDensity, args=[0.1/SynchrotronCriticalxSync(CurrentLowestga),5,CurrentValuesForgaIt,CurrentValuesForNElectrons], kwds={'RelativeError': 0.1})) # SyncPsEnergyDensity. This will append a number to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(SyncPsEnergyDensity, args=[x0/SynchrotronCriticalxSync(CurrentLowestga),5,CurrentValuesForgaIt,CurrentValuesForNElectrons], kwds={'RelativeError': 0.1})) # SyncPsEnergyDensity. This will append a number to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(SyncPsNumberDensity, args=[10**(-6),5,CurrentValuesForgaIt,CurrentValuesForNElectrons], kwds={'RelativeError': 0.1})) # SyncPsNumberDensity. This will append a number to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(SyncPsNumberDensity, args=[0.1/SynchrotronCriticalxSync(CurrentLowestga),5,CurrentValuesForgaIt,CurrentValuesForNElectrons], kwds={'RelativeError': 0.1})) # SyncPsNumberDensity. This will append a number to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(SyncPsNumberDensity, args=[x0/SynchrotronCriticalxSync(CurrentLowestga),5,CurrentValuesForgaIt,CurrentValuesForNElectrons], kwds={'RelativeError': 0.1})) # SyncPsNumberDensity. This will append a number to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(EvaluateLowerBorderOfNonVanRangeOfCOfA1AndCOfA21AndElectronLossRate, args=[ValuesForxSyncSuperIt,ValuesFornSyncPs,CurrentValuesForgaIt,SuperItCounter])) # EvaluateLowerBorderOfNonVanRangeOfCOfA1AndCOfA21AndElectronLossRate. This will plot and save the electron-loss-rate and append None to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(EvaluatePOfB11AndHEPhotonLossRate, args=[ValuesForxSyncSuperIt,ValuesFornSyncPs])) # EvaluatePOfB11AndHEPhotonLossRate. This will plot and save the HEP-loss-rate and append None to the result-list.
    FinalCalculationResultQueue.append(PoolOfFinalCalculation.apply_async(EvaluateSyncPsLossRate, args=[CurrentValuesForgaIt,CurrentValuesForNElectrons])) # EvaluateSyncPsLossRate. This will plot and save the SyncP-loss-rate and append None to the result-list.
    
    FinalCalculationResult = [r.get() for r in FinalCalculationResultQueue] # Compile the results into a list.
    PoolOfFinalCalculation.close() # Clear the Pool-object to exit the processes and to...
    PoolOfFinalCalculation.join() # ...retrieve the used RAM.

    ValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, ValuesFornHEPhotonsSpectralInxga = FinalCalculationResult[0]
    print("\nElectrons' energy-density:                        ", FinalCalculationResult[1]*me*c**2, "J/m^3")
    print("    Corresponding electrons' luminosity:          ", FinalCalculationResult[1]*me*c**2*FactorEnergyDensityToLuminosity('InteractionRadius')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('InteractionRadius')[0])
    print("    Corresponding electrons' luminosity:          ", FinalCalculationResult[1]*me*c**2*FactorEnergyDensityToLuminosity('RClouds')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('RClouds')[0])
    print('\nIC-created HE-photon energy-density:              ', FinalCalculationResult[2]*me*c**2, "J/m^3")
    print('    Corresponding IC-created HE-photon luminosity:', FinalCalculationResult[2]*me*c**2*FactorEnergyDensityToLuminosity('InteractionRadius')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('InteractionRadius')[0])
    print('    Corresponding IC-created HE-photon luminosity:', FinalCalculationResult[2]*me*c**2*FactorEnergyDensityToLuminosity('RClouds')[1], "W (using R=%s)" % FactorEnergyDensityToLuminosity('RClouds')[0])
    print('\nIC-created HE-photon power-density:               ', FinalCalculationResult[3], "1/(m^3*s)")   
    print('    Corresponding IC-created HE-photon power:     ', SphericalVolume('InteractionRadius')[1]*FinalCalculationResult[3]*me*c**2, "W (using R=%s)" % SphericalVolume('InteractionRadius')[0])
    print('    Corresponding IC-created HE-photon power:     ', SphericalVolume('RClouds')[1]*FinalCalculationResult[3]*me*c**2, "W (using R=%s)" % SphericalVolume('RClouds')[0])  
    TrueSyncPsEnergyDensity = FinalCalculationResult[4]*me*c**2 # The whole energy-density of sync-radiation in units of J/m^3.
    LostSyncPsEnergyDensity1 = FinalCalculationResult[5]*me*c**2 # The energy-density of sync-radiation, that is lost by resetting the upper cut-off to 0.1, in units of J/m^3.
    LostSyncPsEnergyDensity2 = FinalCalculationResult[6]*me*c**2 # The energy-density of sync-radiation, that is lost by resetting the upper cut-off to x0, in units of J/m^3.
    print("\nSyncPs' total, true energy-density:", TrueSyncPsEnergyDensity,"J/m^3")
    print('Fraction of energy-density, that is lost by resetting the upper cut-off to 0.1:', LostSyncPsEnergyDensity1/TrueSyncPsEnergyDensity)        
    print('Fraction of energy-density, that is lost by resetting the upper cut-off to x0:', LostSyncPsEnergyDensity2/TrueSyncPsEnergyDensity)        
    TrueSyncPsNumberDensity = FinalCalculationResult[7] # The whole number-density of sync-radiation in units of 1/m^3.
    LostSyncPsNumberDensity1 = FinalCalculationResult[8] # The number-density of sync-radiation, that is lost by resetting the upper cut-off to 0.1, in units of 1/m^3.
    LostSyncPsNumberDensity2 = FinalCalculationResult[9] # The number-density of sync-radiation, that is lost by resetting the upper cut-off to x0, in units of 1/m^3.
    print("SyncPs' total, true number-density:", TrueSyncPsNumberDensity,"/m^3")
    print("Fraction of number-density, that is lost by resetting the upper cut-off to 0.1:", LostSyncPsNumberDensity1/TrueSyncPsNumberDensity)
    print("Fraction of number-density, that is lost by resetting the upper cut-off to x0:", LostSyncPsNumberDensity2/TrueSyncPsNumberDensity)
    
    #print('\nCalling EvaluateConversionOfData:');ValuesForEnergyGeV, ValuesForObsFluxGeVPercm2s, InterpolatedObjectForValuesForObsFluxGeVPercm2s, ValuesForTransmittedObsFluxGeVPercm2s = EvaluateConversionOfData(ValuesForxgammanHEPhotons, ValuesFornHEPhotonsSumSpectralInxga, ValuesFornHEPhotonsSpectralInxga)
    #print('\nCalling EvaluateExportEnergyFluxDensity:');EvaluateExportEnergyFluxDensity(ValuesForEnergyGeV,ValuesForObsFluxGeVPercm2s,ValuesForTransmittedObsFluxGeVPercm2s)
        
    EndingTimeOfEntireEvaluation = time.time() # Now, determine the point of time after the end of the evaluation.
    TimeIntervalOfEntireEvaluation = (EndingTimeOfEntireEvaluation-StartingTimeOfEntireEvaluation)/3600.0 # Determine the time in hours, that was spent by the evaluation. 
    print("\nTime spent on the entire evaluation:", TimeIntervalOfEntireEvaluation, "hours\n")

if __name__ == "__main__": # This prevents evaluation in child-processes.
    print('\n_________________________________________________________________________________________')
    print('_________________________________________________________________________________________\n')
    print('-----------------------------      End of programme       -------------------------------')
    print('_________________________________________________________________________________________')
    print('_________________________________________________________________________________________\n\n')


# Documentation of version 2: IterationCounter was added to the arguments of SecondTermOfNumeratorOfEq8 and to NElectronsOfEq8. The interpolation within _EvaluateNElectronsOfEq8 was added. Furthermore _ImportDataOfIteration() was added.

# Documentation of version 3: The fitting around the discontinuity was added. ValuesForgaLogComputation was added to the arguments of SecondTermOfNumeratorOfEq8 and to NElectronsOfEq8.

# Documentation of version 4: Multiprocessing was added in SecondTermOfNumeratorOfEq8. The part for the evaluation in the shell was added. Furthermore, for all functions, beginning with underscore, the underscore was deleted and the function-name was prefixed with "Evaluate", if not done yet. Additionally, the functions EvaluateNElectronsOfEq8, EvaluateExportDataComplete, EvaluateImportDataComplete, EvaluatePlotForNElectronsComplete, EvaluateSmearingOutTheStepViaFitting, EvaluatePlotForNElectronsSmoothed and EvaluateEquation5 have been adapted in such a way as to be usable via function-call and to make them able to access and alter the global variables. Minor changes.

# Documentation of version 5: The multiprocessing implementation didn't work in version 11. To be precise, the InterpolatedObject could not be pickled. This means that it could not be imported into the child-processes. Therefore the program had to be changed such that not the InterpolatedObject but only the values, that are interpolated, are imported into the child-processes. Then, the interpolation has to be done within the child-processes.

# Documentation of version 6: IntelligentStorage and ApplyIntelligentStorage have been added. Instead, the complicated and intricate procedure within SecondTermOfNumeratorOfEq8 to use AuxiliaryQuantityOfSecondTermOfNumeratorOfEq8 could be erased. The if-test in EvaluateExportInputValueshas been deleted as it is not needed any more.

# Documentation of version 7: In SecondTermOfNumeratorOfEq8, LowestIntegrationBorder was lowered to exactly 1.0/(4*x0). In the smoothing-subroutine, the default value for the finding of the intersection was adjusted to 1.2/(4*x0). A time-measurement of the entire evaluation was added.

# Documentation of version 8: Up to now, there occurred differences in the results between python2- and python3-evaluations. These were due to python2 treating 1/2 as integer-division. To find remedy / make the code also valid for python2, all numbers that are definitively no integer are rendered with a decimal point. Furthermore, quantities that can be simplified, like e. g. 1/2, are simplified, e. g. to 0.5.
# SearchAFile0rDirectory was generalised for directories.
# The epsrel-keyword-argument has been introduced in some integrals.
# nHEPhotons and EvaluatenHEPhotons were added.
# Minor changes.

# Documentation of v1_0 (valid for all three parts):
# The module Concerning_1997ApJ477585M_v1_0.py is included. To to this, it is imported and superfluous import-statements and constant-definitions are removed.
# Due to import-, searching-a-file- and multiprocessing-problems, the function GoToFile was renamed to SearchAFile0rDirectoryEverywhere and the GoToFile-call, to find the file that is to be imported, had to be removed. Now, all the four files have to be situated in the same directory. Furthermore, Pyzo can do multiprocessing, too, now.
# RelativeError is introduced as an argument (given to epsrel in the end) into those functions that build upon an integration. In the iteration, RelativeError is successively decreased from high to lower values.
# The complete structure with IntelligentStorage and ApplyIntelligentStorage in part3 is removed. Instead, the old way of storing SecondTermOfNumeratorOfEq8 in a global variable is used. Why was this done? If an input value in part1 is changed, then, in any case, the complete session has to be restarted and the complete code has to be reimported. By this, the purpose of ApplyIntelligentStorage to realise whether the input values have been modified can not be used anyway.
# n0AF for a general accretion flow was introduced. COfA1 and POfB1 were adapted.
# The definitions of NElectronsga0 and of UpperCutOffOf4thTermOfEq1 have been imprecise up to now and have been revised now based on "Determining an upper cut-off of NElectrons (revised, more exact).png" and on "Determining an upper cut-off of the fourth term (revised, more exact).png".
# The section "Further input values" was added.
# The cascade equation was altered such that it includes an escape-term for the HE photons. Essentially, HEPhotonDestructionRate and NormalisedSpectralPPProbability have been added and POfB11 and pOfB18 have been substituted by them in the terms of the kinetic equation and in the computation of nHEPhotons. An if-test has been added to POfB11 such that it adopts the value 0.0 for all xga<1/x0, where no pair-production is possible.
# It was realised that HighestgaOfEq8 in the function EvaluateNElectronsOfEq8 is still at 1/(4*x0). It is lowered to 1/(8*x0), so that it matches to Lowestga.
# The computation of the total energy-density-rate of the injected photons and of the IC-created ones was introduced.
# n0MultiDelta was implemented.
# The computation of IntegralOfDenominatorOfEq8 was shifted from section "Considering the denominator of equation 8" of part 3 to part 1, where a new section "Energy-density of soft background photons" was created.
# IntegrandOfA1 and IntegrandOfB1 have been refined, such that the value of the function is =0 everywhere the function would otherwise adopt negative values or where it isn't defined at all.
# SecondHighestgaLog is introduced as a fix of the sparse sampling problem (The converging NElectrons changes if the sampling range is changed.). The sampling range is now constructed such that SecondHighestgaLog is always the second highest point and separated from Highestga by a constant distance.
# The optical depth determinations were added.
# Electron escape was included. Especially, ElectronDisappearanceRate, NormalisedSpectralICScatteringProbability and ElectronSpectralDisappearanceDensityRate were added and NextNElectrons was modified.
# In EvaluateSmearingOutTheStepViaFitting the fitting function was changed from a smoothly broken power-law with linear change of the exponent to logarithmic change of exponent.
# In the iteration, the rounding of the number of sampling points was removed. Furthermore, IntermediateDivisionPointgaLog and ValuesForgaLogVeryUpperRange were introduced.

# Documentation of v1_1:
# The new iteration scheme 'PointsFromRToLIterationPointwise' was implemented for evaluation without multiprocessing and almost the whole code was revised accordingly. Iteration is performed into the Thomson-regime (hence below 1/(8*x0)).

# Documentation of v1_2:
# The new iteration scheme 'PointsFromRToLIterationPointwise' was implemented for evaluation with multiprocessing.

# Documentation of v1_3:
# All evaluating function-calls and prints and the 'Evaluate...'-function-definitions have been shifted into a test-block, so that they are not evaluated by child-processes.
# An additional test was introduced in MoreIterationsNecessary, so that there is always at least one additional iteration step evaluated at a certain ga after the iteration at higher values of ga was finished.
# Minor changes.

# Documentation of v1_4:
# Multiprocessing was introduced in FourthTermOfEq1.
# The sampling of ValuesForgaLogComputation was changed.

# Documentation of v1_5:
# The codes Concerning_1997ApJ477585M.py (ADAF-model), Estimates_5_after_MiniWorkshop.py (about Dorit's gap-jet-SSC-toy-model), Mrk501_estimates_version1.py (Toy-model, linking a gap-model with my IC-pair-cascade-model for Mrk501) and the 3 parts Concerning_1988ApJ335786Z_v1_4.py (IC-pair-cascade-model based on Zdz.) were linked to one complete code. The necessary imports have been set and additionally, all the files have been moved to the directory "Astro - Python Code". Overlapping, superfluous or contradicting definitions have been revise or renamed.

# Documentation of v1_6:
# In the main iteration (EvaluateIteration) the multiprocessing implementation via Pool is substituted by Process. This is necessary because Pool is used in FourthTermOfEq1, too, and a Pool-child-process cannot spawn subordinate Pool-processes but a Process-child-process can spawn subordinate Pool-processes. This modification is implemented only in the 'PointsFromRToLIterationPointwise' iteration-scheme.
# The sampling of ValuesForgaLogComputation was modified.

# Documentation of v2_0:
# Synchrotron-radiation has been incorporated, in respect to as a mechanism of electron energy-losses, however not yet as a source of soft photons.

# Documentation of v2_1:
# In the functions FourthTermOfEq1Alternative1, FourthTermOfEq1Alternative2 and FourthTermOfEq1, the intelligent sampling of integration borders has been introduced via UsedSampleIntBordersIn4thTerm. The various realisations of SampleIntegrationBorders are roughly said designed such that each additional order of magnitude of the integration range is covered by one additional integration subrange. Furthermore, the lowest subrange and the highest subrange are narrower, because the integrand can be very steep here. If multiprocessing is used, the number of used subprocesses is limited by the minimum of the number of available CPUs as well as the number of supranges (subtasks).
# In SynchrotronSpectralProductionRate and in SynchrotronSpectralLossRate, the if-tests have been modified to prevent the comparison of an array with None.

# Documentation of v2_2:
# SynchrotronPhotonEscapeRate, NonRelativisticGyroFrequency, SynchrotronCriticalFrequencyCoefficient, SynchrotronCriticalFrequency, SynchrotronCriticalxSyncCoefficient, SynchrotronCriticalxSync, SynchrotronEmissivityIntegrand, SynchrotronEmissivityCoefficient, SynchrotronEmissivity and EvaluateSynchrotronEmissivity were added.
# part6 was renamed to part7, and a new part6 was created for the kinetic equation of the synchrotron-photons. Therein SyncPsSpectralProductionRateIntegrand, SyncPsSpectralProductionRate, SyncPsSpectralICScatteringRate, SyncPsICScatteringRate, SyncPsSpectralICScatteringRateAlt, SyncPsICScatteringRateAlt, SyncPsTotalICScatteringRateIntegrand, SyncPsTotalICScatteringRate, SyncPsLossRate, SyncPsOpticalDepth, SyncPsSpectralNumberDensity, SyncPsNumberDensity, SyncPsSpectralEnergyDensity, SyncPsEnergyDensity, EvaluateSyncPsSpectralProductionRateAndEnergyDensity, EvaluateSyncPsSpectralICScatteringRate, EvaluateSyncPsICScatteringRate, EvaluateSyncPsSpectralICScatteringRateAlt, EvaluateSyncPsICScatteringRateAlt and EvaluateSyncPsLossRate were added.
# The abbreviation EAstOverE was added in IntegrandOfA1. The abbreviation EAstOverE was added in IntegrandOfB1. From these integrands, an additional division by 4 was absorbed into IntegrandOfA1Coefficient. Furthermore, in these integrands, the explicit function expressions were substituted by previous definitions.
# The super-iteration was began to be introduced.

# Documentation of v2_3:
# The building of the super-iteration is resumed.
# CreateANewNElectronsIterated was introduced. IntegrandOfA1WithSynchrotron and COfA1WithSynchrotron were created. COfA1 was renamed to COfA1WithoutSynchrotron and a new, generally usable COfA1 was introduced. All functions, that build upon COfA1, are revised. IntegrandOfA21WithSynchrotron and COfA21WithSynchrotron were created. IntegrandOfCOfA21 was renamed to IntegrandOfCOfA21WithoutSynchrotron and a new, generally usable IntegrandOfCOfA21 was introduced.COfA21 was renamed to COfA21WithoutSynchrotron and a new, generally usable COfA21 was introduced. Similarly, the pair-production equations were generalised for the incorporation of synchrotron-radiation.
# Minor changes.
# In the functions COfA21WithoutSynchrotron, COfA21WithSynchrotron, COfA21, POfB11WithoutSynchrotron, POfB11WithSynchrotron and POfB11, the integration border sampling algorithm was changed from the preliminary logarithmic sampler to SampleIntegrationBorders3. An error message, that occurred now and then, does not occur now any more.
# COfA21Alternative was added and seems to be slightly faster than COfA21. POfB11Alternative was added.
# Multiprocessing was applied in the final calculation.
# In the IC-scattering section, everything concerning Dotgamma (energy-loss-rate of the electrons) was introduced.

# Documentation of v2_4:
# ValuesForgaLogComputation is substituted by ValuesForgaIt.
# Up to now, xSync0Used was always set equal to x0. Now, xSync0Used is determined by xSync0. This might lead to a different NElectronsga0 and this might again lead to a changed xSync0. Thus, the complete code is revised to implement this. Especially, a changing sampling of ValuesForgaIt is introduced.
# Furthermore xgammaLimitCleaned was introduced.
# It was realised that the usage of the Gaussian in the definition of xgammaLimit can be abolished again, because the function pOfB18 (which is now NormalisedSpectralPPProbability due to the inclusion of escape) has no singularity any more. Thus, in xgammaLimitClean (which is used in the 4th term) xgammaLimitOriginal is used again instead of xgammaLimit.
# Minor changes.
# The structure of the code was rearranged such that the temporary results of an interrupted SuperIt can be used to resume the SuperIt. To import the intermediate belay results, one has to manually insert the temporary results of that common iteration (say the n. common iteration), that has been interrupted, into the latest saved NElectronsSuperIt-dictionary (say SuperItCounter = n-1)

# Documentation of v2_5:
# It was realised, that NormalisedSpectralPPProbability in dependence on xga (integration direction of 4th term) has jumps (additional contributions or troughs) at all the points xga=xgammaLimitClean(i,gamma) and at all xga=1/i with i in [xSync0Used,x0MultiDelta,x0]. These values are added as additional integration borders of the integration range of the 4th term.
# It was realised, that the computation of the first integration range of the 4th term sometimes takes a lot of time, while its result (contribution to the whole integral) seems almost negligible. Thus, from now on two possibilities can be realised: If DefaultTake1stIntRangeOf4thTermIntoAccount=False, the 1st range of the 4th term is only determined for those iteration steps, whose relative change is below a certain threshold or whose number of integration ranges is less than 4. This change was (mainly) realised via modifying CreateANewNElectronsIterated, MoreIterationsNecessary and all functions that build upon the 4th term. If DefaultTake1stIntRangeOf4thTermIntoAccount=True, the 1st range is taken into account in all integrations.
# It was realised that COfA1 has kinks due to the discontinuous contributions of n0(x). The kinks of COfA1 have been inserted as additional integration borders into the ThirdTerm as well as in COfA21 as well as in Dotgamma as well as in IntegralOfBracketsOfEq1.
# IntegrandOfA21WithoutSynchrotron and IntegrandOfA21WithSynchrotron have been added.
# ArtificialPushUpOfInitialisation was added.
# In COfA1WithSynchrotron and in POfB1WithSynchrotron, the NumericalBiggestIntegrationBorder was introduced to make the integration range smaller and thus save computation time.
# ADAFFOf28Correct was corrected.
# Minor changes.

# Things to do:
# Check the difference between COfA21Alternative and COfA21.