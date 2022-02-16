## ------------------------------------------------------------------------------------
## Spectra of advection-dominated accretion flows
##
## according to 1997ApJ...477..585M
## ------------------------------------------------------------------------------------

import os
import sys
import operator
import multiprocessing
import types
import time
from datetime import datetime
import pylab as pl
import matplotlib
matplotlib.use('svg') # This had to be inserted, to prevent conflicting GUIs in the plotting functions in the parallelised FinalCalculation. The placement of this call is decisive. It can create a UserWarning. However, the placement is valid like this.
import numpy as np
from scipy import integrate
from scipy import optimize
from scipy.interpolate import interp1d
from scipy.special import kn as ModifiedBesselSecondKindIntegerOrder
from scipy.special import kv as ModifiedBesselSecondKindRealOrder
from scipy.special import erfc as ComplementaryErrorFunction
from scipy.special import exp1 as ExpIntegralFunction
from queue import Empty


## Auxiliary definitions

class ImportFailedError(Exception):
    def __init__(self):
        print('Unable to import or search the specified file.')

class WrongFunctionCallError(Exception):
    def __init__(self, FunctionName):
        print('Wrong call of %s.' % FunctionName)

class InvalidOrderOfADAFParametersError(Exception):
    def __init__(self):
        print("Attention: The ADAF's borders are in an invalid order.")

class InvalidSetOfParametersError(Exception):
    def __init__(self, Specification):
        print('Attention: The combination of the given parameters or their specific values are invalid. %s' % Specification)

def Heavi(LowerStep,x,UpperStep):
    '''Creation of some kind of Heaviside-function. Below LowerStep and above UpperStep it is equal to zero. Between LowerStep and UpperStep it is equal to 1. At LowerStep and at UpperStep it is equal to 1.0, too. x can be an array of numbers or just one number.'''
    return np.where(operator.and_(LowerStep <= x,x <= UpperStep), 1.0, 0.0)

def SampleIntegrationBorders1(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. It creates a comparably moderate number of borders with comparably (with respect to SampleIntegrationBorders4) broad outer integration ranges.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    if LogDifference <= 1: # If the borders are separated by less than one or by one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),3) # ...sample three borders and two integration ranges.
    else: # If the borders are separated by more than one order of magnitude...
        NumberOfBorders = np.ceil(LogDifference)+1 # ...sample a number of borders, which increases with increasing logarithmic difference...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        AdditionalBorder1 = 10**((2.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/3.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 2 in logarithmic space.
        AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+2.0*np.log10(ListOfIntegrationBorders[-1]))/3.0) # And one, which is situated between the second highest and the highest border in the ratio 2 to 1 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def SampleIntegrationBorders2(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. It creates a comparably small number of borders with comparably (with respect to SampleIntegrationBorders3) broad outer integration ranges.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    if LogDifference <= 1: # If the borders are separated by less than one or by one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),3) # ...sample three borders and two integration ranges.
    else: # If the borders are separated by more than one order of magnitude...
        NumberOfBorders = np.ceil(LogDifference) # ...sample a number of borders, which increases with increasing logarithmic difference...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        AdditionalBorder1 = 10**((2.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/3.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 2 in logarithmic space.
        AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+2.0*np.log10(ListOfIntegrationBorders[-1]))/3.0) # And one, which is situated between the second highest and the highest border in the ratio 2 to 1 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def SampleIntegrationBorders3(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. It creates a comparably small number of borders with comparably (with respect to SampleIntegrationBorders2) narrow outer integration ranges.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    if LogDifference <= 1: # If the borders are separated by less than one or by one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),3) # ...sample three borders and two integration ranges.
    else: # If the borders are separated by more than one order of magnitude...
        NumberOfBorders = np.ceil(LogDifference) # ...sample a number of borders, which increases with increasing logarithmic difference...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        AdditionalBorder1 = 10**((3.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/4.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 3 in logarithmic space.
        AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+3.0*np.log10(ListOfIntegrationBorders[-1]))/4.0) # And one, which is situated between the second highest and the highest border in the ratio 3 to 1 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def SampleIntegrationBorders4(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. It creates a comparably moderate number of borders with comparably (with respect to SampleIntegrationBorders1) narrow outer integration ranges.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    NumberOfBorders = np.ceil(LogDifference)+1 # ...sample a number of borders, which increases with increasing logarithmic difference...
    ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
    AdditionalBorder1 = 10**((3.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/4.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 3 in logarithmic space.
    AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+3.0*np.log10(ListOfIntegrationBorders[-1]))/4.0) # And one, which is situated between the second highest and the highest border in the ratio 3 to 1 in logarithmic space.
    ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def SampleIntegrationBorders5(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. It creates a comparably big number of borders with comparably (with respect to SampleIntegrationBorders6) broad outer integration ranges.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    if LogDifference <= 0.2: # If the borders are separated by far less than one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),3) # ...sample three borders and two integration ranges.
    elif LogDifference <= 1: # If the borders are separated by less than one or by one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),2) # ...sample 2 borders.
        AdditionalBorder1 = 10**((3.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/4.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 3 in logarithmic space.
        AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+3.0*np.log10(ListOfIntegrationBorders[-1]))/4.0) # And one, which is situated between the second highest and the highest border in the ratio 3 to 1 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    else: # If the borders are separated by more than one order of magnitude...
        NumberOfBorders = np.ceil(LogDifference)+2 # ...sample a number of borders, which increases with increasing logarithmic difference...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        AdditionalBorder1 = 10**((2.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/3.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 2 in logarithmic space.
        AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+2.0*np.log10(ListOfIntegrationBorders[-1]))/3.0) # And one, which is situated between the second highest and the highest border in the ratio 2 to 1 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def SampleIntegrationBorders6(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. It creates a comparably big number of borders with comparably (with respect to SampleIntegrationBorders5) narrow outer integration ranges.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    if LogDifference <= 0.3: # If the borders are separated by far less than one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),3) # ...sample three borders and two integration ranges.
    elif LogDifference <= 1: # If the borders are separated by less than one or by one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),2) # ...sample 2 borders.
        AdditionalBorder1 = 10**((3.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/4.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 3 in logarithmic space.
        AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+3.0*np.log10(ListOfIntegrationBorders[-1]))/4.0) # And one, which is situated between the second highest and the highest border in the ratio 3 to 1 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    else: # If the borders are separated by more than one order of magnitude...
        NumberOfBorders = np.ceil(LogDifference)+2 # ...sample a number of borders, which increases with increasing logarithmic difference...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        AdditionalBorder1 = 10**((3.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/4.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 3 in logarithmic space.
        AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+3.0*np.log10(ListOfIntegrationBorders[-1]))/4.0) # And one, which is situated between the second highest and the highest border in the ratio 3 to 1 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.
    
def SampleIntegrationBorders7(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. It creates a comparably big number of borders with comparably (with respect to SampleIntegrationBorders5) narrow outer integration ranges. The difference to SampleIntegrationBorders6 is that it uses more borders at the smallest differences between LowestIntegrationBorder and BiggestIntegrationBorder.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    NumberOfBorders = np.ceil(LogDifference)+2 # Sample a big number of borders, which increases with increasing logarithmic difference...
    ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
    AdditionalBorder1 = 10**((3.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/4.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 3 in logarithmic space.
    AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+3.0*np.log10(ListOfIntegrationBorders[-1]))/4.0) # And one, which is situated between the second highest and the highest border in the ratio 3 to 1 in logarithmic space.
    ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def SampleIntegrationBorders8(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. The difference to SampleIntegrationBorders7 is that it uses additional borders immediately below BiggestIntegrationBorder and that it does not use the additional border above LowestIntegrationBorder. It is tailored for SyncPsICScatteringRate, because the integrand SyncPsSpectralICScatteringRate has a sharp peak below BiggestIntegrationBorder. The peak gets sharper with increasing BiggestIntegrationBorder'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    NumberOfBorders = np.ceil(LogDifference)+2 # Sample a big number of borders, which increases with increasing logarithmic difference...
    ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
    AdditionalBorder1 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+3.0*np.log10(ListOfIntegrationBorders[-1]))/4.0) # Insert an additional border, which is situated between the second highest and the highest border in the ratio 3 to 1 in logarithmic space.
    if np.log10(BiggestIntegrationBorder)>3:
        AdditionalList=[10**((1.0*np.log10(ListOfIntegrationBorders[-2])+(10**(i+1)-1)*np.log10(ListOfIntegrationBorders[-1]))/10**(i+1)) for i in range(0,int(np.ceil(np.log10(BiggestIntegrationBorder)))-3)] # This list samples additional borders. It samples the more borders, the higher BiggestIntegrationBorder is. They are situated between the second highest and the highest border in the ratio 10^(i+1)-1 to 1 in logarithmic space.
        ListOfIntegrationBorders=np.append(np.asarray(AdditionalList),ListOfIntegrationBorders)
    ListOfIntegrationBorders=np.append(np.asarray([AdditionalBorder1]),ListOfIntegrationBorders)
    ListOfIntegrationBorders = np.sort(ListOfIntegrationBorders)
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def SampleIntegrationBorders9(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. It creates a comparably big number of borders. The difference to SampleIntegrationBorders7 is that it uses more borders above LowestIntegrationBorder.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    NumberOfBorders = np.ceil(LogDifference)+2 # Sample a big number of borders, which increases with increasing logarithmic difference...
    ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
    AdditionalBorder1 = 10**((2.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/3.0) # Insert three additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 2 in logarithmic space.
    AdditionalBorder2 = 10**((8.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/9.0) # One, which is situated between the lowest and the second lowest border in the ratio 1 to 8 in logarithmic space.
    AdditionalBorder3 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+3.0*np.log10(ListOfIntegrationBorders[-1]))/4.0) # And one, which is situated between the second highest and the highest border in the ratio 3 to 1 in logarithmic space.
    ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2,AdditionalBorder3]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def SampleIntegrationBorders10(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. Up to v2_4, there was the necessity to use tendentially a big number of small ranges in the 4th term of equation 4, especially immediately above LowestIntegrationBorder. In v2_5 additional, physically motivated borders are used in the 4th term. Thus the division of the remaining ranges does not need to be that dense. This is realised now. This function creates a comparably small number of borders. The difference to SampleIntegrationBorders3 is that it uses no border for a small difference between BiggestIntegrationBorder and LowestIntegrationBorder.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    if LogDifference <= 0.4: # If the borders are separated by far less than one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),2) # ...sample two borders and one integration range.
    elif LogDifference <= 1: # If the borders are separated by less than one or by one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),3) # ...sample three borders and two integration ranges.
    else: # If the borders are separated by more than one order of magnitude...
        NumberOfBorders = np.ceil(LogDifference) # ...sample a number of borders, which increases with increasing logarithmic difference...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        AdditionalBorder1 = 10**((3.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/4.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 3 in logarithmic space.
        AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+3.0*np.log10(ListOfIntegrationBorders[-1]))/4.0) # And one, which is situated between the second highest and the highest border in the ratio 3 to 1 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def SampleIntegrationBorders11(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. This is based on SampleIntegrationBorders10. However, there are two modifications: First, the difference between BiggestIntegrationBorder and LowestIntegrationBorder, where no additional border is introduced, is made even smaller. Second, the lowest integration range is made narrower. The motivation for these changes was the observation with SampleIntegrationBorders10, that for big gamma, the list of integration borders is e.g. like this [5790055.135, 5793196.974, 5799650.972, 5800643.581, 5802598.772, 11244641.201] and the highest range takes a lot of time. It was then found that division of the highest range in the ratio 1:5 or 1:6 is helpful.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    if LogDifference <= 0.1: # If the borders are separated by far less than one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),2) # ...just sample two borders and one integration range.
    elif LogDifference <= 1: # If the borders are separated by less than one or by one order of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),2) # ...sample two borders and...
        AdditionalBorder1 = 10**((6.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[-1]))/7.0) # ...add one border, which is situated between the lowest border and the highest border in the ratio 1 to 6 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1]),ListOfIntegrationBorders))
    else: # If the borders are separated by more than one order of magnitude...
        NumberOfBorders = np.ceil(LogDifference) # ...sample a number of borders, which increases with increasing logarithmic difference...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        AdditionalBorder1 = 10**((4.0*np.log10(ListOfIntegrationBorders[0])+1.0*np.log10(ListOfIntegrationBorders[1]))/5.0) # Insert two additional borders. One, which is situated between the lowest and the second lowest border in the ratio 1 to 4 in logarithmic space.
        AdditionalBorder2 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+3.0*np.log10(ListOfIntegrationBorders[-1]))/4.0) # And one, which is situated between the second highest and the highest border in the ratio 3 to 1 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1,AdditionalBorder2]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def SampleIntegrationBorders12(LowestIntegrationBorder,BiggestIntegrationBorder):
    '''This divides the integration range reaching from LowestIntegrationBorder to BiggestIntegrationBorder into a number of narrower integration ranges. It creates a very small number of borders and has one additional border at the upper end. It is tailored for the sampling of UsedSampleIntBordersOfxSync.'''
    LogDifference = np.log10(BiggestIntegrationBorder)-np.log10(LowestIntegrationBorder)
    if LogDifference <= 3.0: # If the borders are separated by less than three or by three orders of magnitude...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),3) # ...sample three borders and two integration ranges.
    else: # If the borders are separated by more than three orders of magnitude...
        NumberOfBorders = np.floor(LogDifference) # ...sample a number of borders, which increases with increasing logarithmic difference...
        ListOfIntegrationBorders = np.logspace(np.log10(LowestIntegrationBorder),np.log10(BiggestIntegrationBorder),NumberOfBorders)
        AdditionalBorder1 = 10**((1.0*np.log10(ListOfIntegrationBorders[-2])+1.0*np.log10(ListOfIntegrationBorders[-1]))/2.0) # Add one border, which is situated between the second highest and the highest border in the ratio 1 to 1 in logarithmic space.
        ListOfIntegrationBorders = np.sort(np.append(np.asarray([AdditionalBorder1]),ListOfIntegrationBorders))
    return ListOfIntegrationBorders, len(ListOfIntegrationBorders)-1 # Return the integration borders as well as the number of ranges, which is also the number of subtasks that are used for the total integration.

def GetCurrentDate():
    '''This function returns a string which represents the current date and time. It is of the format'YYYY-MM-DD HH-MM'.'''
    CurrentYear = datetime.ctime(datetime.today())[-4:]
    MonthToNumber = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04', 'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08', 'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}
    CurrentMonth = MonthToNumber[datetime.ctime(datetime.today())[4:7]]
    CurrentDay = datetime.ctime(datetime.today())[8:10].replace(' ', '0')
    CurrentTime = datetime.ctime(datetime.today())[11:16].replace(':', '-')
    CurrentDate = '-'.join([CurrentYear,CurrentMonth,CurrentDay])
    CurrentDateTime = ' '.join([CurrentDate,CurrentTime])
    return CurrentDateTime

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

def ImportFileToString(NameOfFileToImport):
    '''This function reads everything out of file NameOfFileToImport straight into a string.'''
    PathOfFileToImport=SearchAFile0rDirectoryEverywhere(NameOfFileToImport,False)
    if PathOfFileToImport==None:
        raise ImportFailedError
    else:
        ImportedFile = open('%s' % os.path.join(PathOfFileToImport,NameOfFileToImport), 'r')
        OpenedLines = ImportedFile.readlines()
        ImportedFile.close()
        if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
            print('    Imported this file.')
        return OpenedLines

def DeleteFile(PartOfFilenameToDelete,WhereToSearch='.'):
    '''Search in the directory WhereToSearch and in all its subfolders for files that contain the string PartOfFilenameToDelete in their filename. An example for Filename would be "temporary file" and an example for WhereToSearch would be "C:\\Messy Dir\\". "." denotes the current directory. Then, the files that have been found are deleted.'''
    print('\nSearching and deleting files including the string "%s":' % PartOfFilenameToDelete)
    Filecounter=0
    for dirpath, subdirnames, filenames in os.walk(WhereToSearch):
        for Filename in filenames:
            if PartOfFilenameToDelete in Filename:
                Filecounter+=1
                PathContainingAFileWithPartOfFilenameToDelete=os.path.join(dirpath,Filename)
                os.remove(PathContainingAFileWithPartOfFilenameToDelete)
                print(PathContainingAFileWithPartOfFilenameToDelete)
    print("    Finished with deleting", Filecounter ,"files.")
    
def ColourRedToBlack(i,N):
    '''This returns the rgb tuple where the value of red is decreasing linearly as i goes from 1 to N.'''
    return (np.abs((i-N)/(1-N)),0,0)


## Constants, input-parameters and initialisation

np.set_printoptions(threshold=np.inf, precision=15)

c = 299792458.0   # The velocity of light in units m/s
cGaussian = 29979245800.0   # The velocity of light in units cm/s
h = 6.626070*10**(-34)  # The Planck's Constant in units J*s
hGaussian = 6.626070*10**(-27)  # The Planck's Constant in units erg*s
me = 9.109384*10**(-31)  # The electron's mass in units kg
meGaussian = 9.109384*10**(-28)  # The electron's mass in units g
mProton = 1.672621*10**(-27)  # The proton's mass in units kg
kB = 1.380649*10**(-23)  # The Boltzmann's Constant in units J/K
kBGaussian = 1.380649*10**(-16)  # The Boltzmann's Constant in units erg/K
e = 1.602177*10**(-19)  # The elementary charge in units of C
eGaussian = 4.803204*10**(-10)  # The elementary charge in units of Fr = g^0.5*cm^1.5*s^(-1)
G = 6.674*10**(-11) # The gravitational constant in units of m^3/(kg*s^2)
GGaussian = 6.674*10**(-8) # The gravitational constant in units of cm^3/(g*s^2)
aRadiation = (8*np.pi**5*kB**4)/(15*c**3*h**3) # The radiation-constant in units of J/(m^3*K^4)
mu0 = 4.0*np.pi*10**(-7) # The permeability of vacuum (magnetic constant) in units of N/A^2 = kg*m/(A^2*s^2)
epsilon0 = 1.0/(mu0*c**2.0) # The permittivity of vacuum (electric constant) in units of A^2*s^4/(kg*m^3)
rClassical = e**2.0/(4.0*np.pi*epsilon0*me*c**2.0) # The classical electron radius in units of m. Cf. e.g. Malcolm Longair's "High Energy Astrophysics", chapter 8.1.
sigmaT = 8.0*np.pi*rClassical**2.0/3.0  # The Thomson-cross-section in units m^2
sigmaTGaussian = 8*np.pi*eGaussian**4/(3*meGaussian**2*cGaussian**4)  # The Thomson-cross-section in units cm^2
alphaFine = e**2.0/(2.0*c*epsilon0*h) # The fine-structure constant. It is dimensionless.
OmegaR = 0.0 # The dimensionless energy-density (density-parameter) of radiation.
OmegaM = 0.3 # The dimensionless energy-density (density-parameter) of matter.
OmegaL = 0.7 # The dimensionless energy-density (density-parameter) of the vacuum.
Hubble0 = 70.0 # The contemporary Hubble-constant in units of km/(s*Mpc).

# In what follows, the following parameters are used:
alphaViscosity=0.3 # alphaViscosity is the viscosity parameter. It is without units. INPUT VALUE!
betaPressure=0.5 # betaPressure is the ratio of gas pressure to total pressure as defined in equation 1. It is without units. INPUT VALUE!
ADAFrInner=3.0 # ADAFrInner is the inner border of the accretion disk in units of Schwarzschild-radii. INPUT VALUE!
ADAFrOuter=1000.0 # ADAFrOuter is the outer border of the accretion disk in units of Schwarzschild-radii. INPUT VALUE!
ADAFTemperatureElectron=6.4*10**9 # ADAFTemperatureElectron is the electrons' temperature In the ADAF in units of Kelvin. INPUT VALUE!
M9=1 # M9 is the mass of the central compact object in units of 10^9 solar masses. INPUT VALUE!
Dotm=0.0008 # Dotm is the accretion rate in units of the Eddington-accretion rate. INPUT VALUE!
UseApproxxM=False # UseApproxxM specifies whether the approximated formula for xM is used or whether the more exact numerical determination is used. INPUT VALUE!

etaff = 0.1 # A dimensionless efficiency parameter for the liberation of mass energy in an accretion process. INPUT VALUE!
Toy1Distance = 700 # The distance to the accretion disk in units of Schwarzschild-radii, that is used in the gap-jet-SSC-toy-model (toy-model 1). INPUT VALUE!
B4 = 1.0*10**(-4) # This is the magnetic flux-density in units of 10**4 Gauss, which is one T. INPUT VALUE!
MagneticEnergyDensity = B4**2/(2*mu0) # The energy-density, that is stored in the magnetic field. The unit is J/m^3.
hGap = 0.5 # The vacuum-gap-height in units of Schwarzschild-radii. INPUT VALUE!
rCurvature = 1 # The curvature-radius of the magnetic field-lines in units of Schwarzschild-radii. INPUT VALUE!

# Initialisation of global objects:
ADAFIntelligentStorage = {} # This object stores values for certain quantities. Especially, it stores the result-values of ADAFxMNumerical, ADAFnuMin, ADAFnuPeak and ADAFSpectralLuminosityCyclosynchrotronPeak as well as the input-values which have been used for determining these result-values. This storage is intelligent because when it is applied via ApplyADAFIntelligentStorage, it realises by itself whether it can just recall the saved result-values (with the advantage of saving time) or whether it has to determine the result-values anew.

def ApplyADAFIntelligentStorage(ReturnedQuantity):
    '''This function treats the ADAFIntelligentStorage. Its purpose is to save time, because the evaluation of ADAFxMNumerical, ADAFnucOf21 or of ADAFSpectralLuminosityCyclosynchrotron is very time-intensive. Hence, it is reasonable to store the values as soon as they have been determined and to recall the stored values when they are needed instead of recalling the time-intensive functions again and again.'''
    # This case is evaluated either if the ADAFIntelligentStorage is empty or if one or more of the input-values have been altered. 
    if ADAFIntelligentStorage == {} or ADAFIntelligentStorage['alphaViscosity'] != alphaViscosity or ADAFIntelligentStorage['betaPressure'] != betaPressure or ADAFIntelligentStorage['ADAFrInner'] != ADAFrInner or ADAFIntelligentStorage['ADAFrOuter'] != ADAFrOuter or ADAFIntelligentStorage['ADAFTemperatureElectron'] != ADAFTemperatureElectron or ADAFIntelligentStorage['M9'] != M9 or ADAFIntelligentStorage['Dotm'] != Dotm or ADAFIntelligentStorage['UseApproxxM'] != UseApproxxM:
        # In this case all the input-values are stored in the ADAFIntelligentStorage and...
        ADAFIntelligentStorage['alphaViscosity'] = alphaViscosity
        ADAFIntelligentStorage['betaPressure'] = betaPressure
        ADAFIntelligentStorage['ADAFrInner'] = ADAFrInner
        ADAFIntelligentStorage['ADAFrOuter'] = ADAFrOuter
        ADAFIntelligentStorage['ADAFTemperatureElectron'] = ADAFTemperatureElectron
        ADAFIntelligentStorage['M9'] = M9
        ADAFIntelligentStorage['Dotm'] = Dotm
        ADAFIntelligentStorage['UseApproxxM'] = UseApproxxM
        # ...the result-values are determined and stored in the ADAFIntelligentStorage.
        ADAFIntelligentStorage['ADAFxMNumerical'] = ADAFxMNumerical(alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm)
        ADAFIntelligentStorage['ADAFnuMin'] = ADAFnucOf21(alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,UseApproxxM)
        ADAFIntelligentStorage['ADAFnuPeak'] = ADAFnucOf21(alphaViscosity,betaPressure,ADAFrInner,ADAFTemperatureElectron,M9,Dotm,UseApproxxM)
        ADAFIntelligentStorage['ADAFSpectralLuminosityCyclosynchrotronPeak'] = ADAFSpectralLuminosityCyclosynchrotron(ADAFIntelligentStorage['ADAFnuPeak'],alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,UseApproxxM)
        #print('1st if-test done.')
    # In this part, a stored result-value is returned, depending on the ReturnedQuantity-argument.
    if ReturnedQuantity=='ADAFxMNumerical':
        #print('2nd if-test done, returned', ADAFIntelligentStorage)
        return ADAFIntelligentStorage['ADAFxMNumerical']
    elif ReturnedQuantity=='ADAFnuMin':
        return ADAFIntelligentStorage['ADAFnuMin']
    elif ReturnedQuantity=='ADAFnuPeak':
        return ADAFIntelligentStorage['ADAFnuPeak']
    elif ReturnedQuantity=='ADAFSpectralLuminosityCyclosynchrotronPeak':
        return ADAFIntelligentStorage['ADAFSpectralLuminosityCyclosynchrotronPeak']

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    print('\n')
    print('\n_________________________________________________________________________________________\n')
    print('-----------------------------------    ADAF-model    ------------------------------------')
    print('_________________________________________________________________________________________\n\n')


## Conversions and parameters

def MpcTom(lengthInMpc):
    '''This converts a distance/length from Megaparsec to metres.'''
    return lengthInMpc*10.0**6*3.085678*10.0**16

HubbleSI = Hubble0*1000/MpcTom(1.0) # The contemporary Hubble-constant in units of 1/s.

def mToMpc(lengthInm):
    '''This converts a distance/length from metres to Megaparsec.'''
    return lengthInm/(10.0**6*3.085678*10.0**16)

def OpeningAngleToSolidAngle(OpeningAngleInDegree):
    '''This computes the solid angle (in units of steradians) of a cone with given opening angle (in units of °).'''
    OpeningAngleInRad = 2*np.pi*OpeningAngleInDegree/360
    return 4*np.pi*(np.sin(OpeningAngleInRad/4))**2

ADAFc1 = 0.5 # This is the coefficient c1 from 1995ApJ...452..710N which has this value according to 1997ApJ...477..585M
ADAFc3 = 0.3 # This is the coefficient c3 from 1995ApJ...452..710N which has this value according to 1997ApJ...477..585M

def ADAFs1Of5(alphaViscosity,betaPressure):
    '''The parameter s_1 according to equation 5. It is dimensionless.'''
    return 1.42*10**9*np.sqrt((1-betaPressure)*ADAFc3/(alphaViscosity*ADAFc1))

def ADAFs2Of22(xM):
    '''The parameter s_2 according to equation 22. It is dimensionless. xM is something like a normalised frequency.'''
    return 1.19*10**(-13)*xM

ADAFs3Of23 = 1.05*10**(-24) # This coefficient s_3 is dimensionless according to equation 23.

def Theta(Temperature):
    '''Theta is a dimensionless temperature, namely the energy-equivalent of the temperature Temperature in units of the electrons' rest-energy. In the following, Theta with capital T is a dimensionless temperature and has to be distinguished between theta as an angle.'''
    return kB*Temperature/(me*c**2)

def ToThe9SolarMassesTokg(M9):
    '''Convert the mass M9 to SI-units, i. e. to kg.'''
    return M9*(10**9)*(1.99*10**30)

def rSchwarzschild(M9):
    '''Determine the Schwarzschild-radius of an object of mass M9. It is in units of m.'''
    return 2*G*ToThe9SolarMassesTokg(M9)/c**2

def rSchwarzschildTom(d,M9):
    '''Conversion of a distance d in units of Schwarzschild-radii to that distance in units of m.'''
    return d*rSchwarzschild(M9)

def FrequencyToEnergy(nu):
    '''Conversion of the frequency nu to energy. nu has to be given in units of Hz, while the returned energy is in units of J.'''
    return nu*h

def DotMEddington(M9,etaff):
    '''The Eddington accretion rate in units of kg/s.'''
    return (4*np.pi*G*ToThe9SolarMassesTokg(M9)*mProton)/(sigmaT*etaff*c)

def DotMEddingtonTokgPers(M9,etaff,Dotm):
    '''This converts an accretion rate in units of DotMEddington to kg/s.'''
    return Dotm*DotMEddington(M9,etaff)

def ergPercm3ToJPerm3(EnergyDensityInergPercm3):
    '''This converts an energy-density from unit erg/cm^3 to unit J/m^3'''
    return EnergyDensityInergPercm3/10


## Bremsstrahlung

def ADAFFOf28(ThetaElectron):
    '''This is the auxiliary function of equation 28. ThetaElectron is the electrons' dimensionless temperature. It is dimensionless. Attention: For ADAFTemperatureElectron>1.4*10**10 K, this gets negative and consequently ADAFSpectralLuminosityBremsstrahlungOf30 gets negative, too.'''
    if ThetaElectron<=1:
        return (9*ThetaElectron/(2*np.pi))*(np.log(1.123*ThetaElectron+0.48)+1.5)+2.3*ThetaElectron*(np.log(1.123*ThetaElectron)+1.28)
    elif ThetaElectron>1:
        return 4*(2*ThetaElectron/np.pi**3)**0.5*(1+1.781*ThetaElectron**1.34)+1.73*ThetaElectron**(3.0/2)*(1+1.1*ThetaElectron+ThetaElectron**2-1.25*ThetaElectron**(5.0/2))

def ADAFFOf28Correct(ThetaElectron):
    '''This is the auxiliary function of equation 28. ThetaElectron is the electrons' dimensionless temperature. It is dimensionless. This is the correct definition.'''
    if ThetaElectron<=1:
        return 4*(2*ThetaElectron/np.pi**3)**0.5*(1+1.781*ThetaElectron**1.34)+1.73*ThetaElectron**(3.0/2)*(1+1.1*ThetaElectron+ThetaElectron**2-1.25*ThetaElectron**(5.0/2))
    elif ThetaElectron>1:
        return (9*ThetaElectron/(2*np.pi))*(np.log(1.123*ThetaElectron+0.48)+1.5)+2.3*ThetaElectron*(np.log(1.123*ThetaElectron)+1.28)

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.    
    def EvaluateADAFFOf28Correct():
        X=np.linspace(0.1,2.9,10000)
        Y=[ADAFFOf28(i) for i in X]
        pl.plot(X,Y)
        Y2=[ADAFFOf28Correct(i) for i in X]
        pl.plot(X,Y2)

def ADAFSpectralLuminosityBremsstrahlungOf30Coefficient(alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm):
    '''This is the coefficient of the bremsstrahlung-luminosity-function, in other words it just not includes the exp-term. It is in units of W/Hz. Whether ADAFFOf28Correct or ADAFFOf28 is used makes no difference in the computations for the A&A Mrk501 bump paper.'''
    return 2.29*10**24*np.log(ADAFrOuter/ADAFrInner)*ADAFFOf28Correct(Theta(ADAFTemperatureElectron))*(10**9*M9)*Dotm**2/(alphaViscosity**2*ADAFc1**2*ADAFTemperatureElectron)*10**(-7)

def ADAFSpectralLuminosityBremsstrahlungOf30(nu,alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm):
    '''This is the spectral bremsstrahlung-luminosity of the disk in units of W/Hz. nu is the frequency in Hz. Attention, here, M9 is used in contrast to the original equation.'''
    return ADAFSpectralLuminosityBremsstrahlungOf30Coefficient(alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)*np.exp(-h*nu/(kB*ADAFTemperatureElectron))


## Cyclosynchrotron

def ADAFxMApproxOfB3(Dotm):
    '''This is equation B3 solved for xM. It is dimensionless and is similar to a normalised frequency.'''
    xM = 10**((3.6+0.25*np.log10(Dotm)))
    #print('Approximate xM =', xM)
    return xM

def ADAFB2(y,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm):
    '''This is equation B2, where the left-hand side was shifted to the right, such that one has the form 0=... It is dimensionless.'''
    return 10.36+0.26*(np.log((10**9*M9)*Dotm)-np.log((Theta(ADAFTemperatureElectron))**3*ModifiedBesselSecondKindIntegerOrder(2.0,1.0/Theta(ADAFTemperatureElectron)))-np.log(alphaViscosity*ADAFc1*ADAFc3*(1-betaPressure)/0.0225))-y-1.852*np.log(y)

def ADAFxMNumerical(alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm):
    '''This solves the equation 0=ADAFB2 numerically for y and then substitutes according to y=xM^(1/3). Hence, it yields xM, which is a dimensionless quantity. For reasonable parameters, this is about twice the approximate value. However, the change in the spectrum is quite small.'''
    Root = optimize.fsolve(ADAFB2, 1, args=(alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm))
    Root = Root[0]
    xM = Root**3
    #print('Numerically found xM =', xM)
    return xM

def ADAFnucOf21(alphaViscosity,betaPressure,r,ADAFTemperatureElectron,M9,Dotm,UseApproxxM=False):
    '''The critical frequency or cut-off frequency as defined in eq. 21. It is in units of Hz. Attention, this is not what is commonly defined as the critical frequency of synchrotron-radiation. With the keyword argument UseApproxxM, one can choose whether to use the numerical value for xM or the approximate value. Using UseApproxxM==False is however faster due to the intelligent storage.'''
    if UseApproxxM==False:
        xM = ApplyADAFIntelligentStorage('ADAFxMNumerical')
    elif UseApproxxM==True:
        xM = ADAFxMApproxOfB3(Dotm)
    return ADAFs1Of5(alphaViscosity,betaPressure)*ADAFs2Of22(xM)*np.sqrt(Dotm/(10**9*M9))*ADAFTemperatureElectron**2*r**(-1.25)

def ADAFSpectralLuminosityCyclosynchrotronOf25(nu,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm,UseApproxxM=False):
    '''This is the spectral cyclosynchrotron-luminosity of the disk in units of W/Hz, however this expression is valid only above the lower cut-off nu_min and below nu_peak. nu is the frequency in Hz. Attention, here, M9 is used in contrast to the original equation. With the keyword argument UseApproxxM, one can choose whether to use the numerical value for xM or the approximate value.'''
    if UseApproxxM==False:
        xM = ApplyADAFIntelligentStorage('ADAFxMNumerical')
    elif UseApproxxM==True:
        xM = ADAFxMApproxOfB3(Dotm)
    return ADAFs3Of23*(ADAFs1Of5(alphaViscosity,betaPressure)*ADAFs2Of22(xM))**1.6*(10**9*M9)**1.2*Dotm**0.8*ADAFTemperatureElectron**4.2*nu**0.4*10**(-7)

def ADAFSpectralLuminosityCyclosynchrotronBelownuMin(nu,alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,UseApproxxM=False):
    '''This is the spectral cyclosynchrotron-luminosity of the disk in units of W/Hz, however this expression is valid only below the lower cut-off nu_min. nu is the frequency in Hz. It is constructed such that is has slope 22/13 and that it coincides with the other part exactly at nu_min. Attention, here, M9 is used in contrast to the original equation. With the keyword argument UseApproxxM, one can choose whether to use the numerical value for xM or the approximate value.'''
    if UseApproxxM==False:
        ADAFnuMin = ApplyADAFIntelligentStorage('ADAFnuMin')
        return ADAFSpectralLuminosityCyclosynchrotronOf25(ADAFnuMin,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm)*(nu/ADAFnuMin)**(22.0/13)
    elif UseApproxxM==True:
        ADAFnuMin = ADAFnucOf21(alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,True)
        return ADAFSpectralLuminosityCyclosynchrotronOf25(ADAFnuMin,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm,True)*(nu/ADAFnuMin)**(22.0/13)

def ADAFSpectralLuminosityCyclosynchrotron(nu,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,UseApproxxM=False):
    '''This is the spectral cyclosynchrotron-luminosity of the disk in units of W/Hz. nu is the frequency in Hz. Attention, here, M9 is used in contrast to the original equation. With the keyword argument UseApproxxM, one can choose whether to use the numerical value for xM or the approximate value.'''
    if UseApproxxM==False:
        ADAFnuMin = ApplyADAFIntelligentStorage('ADAFnuMin')
        ADAFnuPeak = ApplyADAFIntelligentStorage('ADAFnuPeak')
        if nu<ADAFnuMin:
            return ADAFSpectralLuminosityCyclosynchrotronBelownuMin(nu,alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)
        elif nu<=ADAFnuPeak:
            return ADAFSpectralLuminosityCyclosynchrotronOf25(nu,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm)
        else:
            return 0
    if UseApproxxM==True:
        ADAFnuMin = ADAFnucOf21(alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,True)
        ADAFnuPeak = ADAFnucOf21(alphaViscosity,betaPressure,ADAFrInner,ADAFTemperatureElectron,M9,Dotm,True)
        if nu<ADAFnuMin:
            return ADAFSpectralLuminosityCyclosynchrotronBelownuMin(nu,alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,True)
        elif nu<=ADAFnuPeak:
            return ADAFSpectralLuminosityCyclosynchrotronOf25(nu,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm,True)
        else:
            return 0


## Comptonised spectrum

def ADAFtauesOf31(alphaViscosity,ADAFrInner,Dotm):
    '''The mean optical depth to inverse-Compton-scattering according to equation 31. It is dimensionless.'''
    return 6.2*Dotm/(alphaViscosity*ADAFc1*np.sqrt(ADAFrInner))

def ADAFAOf32(ADAFTemperatureElectron):
    '''The mean amplification factor of inverse-Compton-scattering according to equation 32, which is dimensionless.'''
    return 1+4.0*Theta(ADAFTemperatureElectron)+16*(Theta(ADAFTemperatureElectron))**2

def ADAFalphacOf34(alphaViscosity,ADAFrInner,ADAFTemperatureElectron,Dotm):
    '''The slope of the comptonised spectrum according to equation 34. It is without units.'''
    return -np.log(ADAFtauesOf31(alphaViscosity,ADAFrInner,Dotm))/np.log(ADAFAOf32(ADAFTemperatureElectron))

def ADAFSpectralLuminosityComptonised(nu,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,UseApproxxM=False):
    '''This is the spectral luminosity due to comptonisation (equation 38) of the cyclosynchrotron-photons in units of W/Hz. nu is the frequency in Hz. Attention, here, M9 is used in contrast to the original equation. With the keyword argument UseApproxxM, one can choose whether to use the numerical value for xM or the approximate value.'''
    if UseApproxxM==False:
        ADAFnuPeak = ApplyADAFIntelligentStorage('ADAFnuPeak')
        if nu<=ADAFnuPeak:
            return 0
        elif nu<=3*kB*ADAFTemperatureElectron/h:
            ADAFSpectralLuminosityCyclosynchrotronPeak = ApplyADAFIntelligentStorage('ADAFSpectralLuminosityCyclosynchrotronPeak')
            return ADAFSpectralLuminosityCyclosynchrotronPeak*(ADAFnuPeak/nu)**ADAFalphacOf34(alphaViscosity,ADAFrInner,ADAFTemperatureElectron,Dotm)
        else:
            return 0
    if UseApproxxM==True:
        ADAFnuPeak = ADAFnucOf21(alphaViscosity,betaPressure,ADAFrInner,ADAFTemperatureElectron,M9,Dotm,True)
        if nu<=ADAFnuPeak:
            return 0
        elif nu<=3*kB*ADAFTemperatureElectron/h:
            ADAFSpectralLuminosityCyclosynchrotronPeak = ADAFSpectralLuminosityCyclosynchrotron(ADAFnuPeak,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,True)
            return ADAFSpectralLuminosityCyclosynchrotronPeak*(ADAFnuPeak/nu)**ADAFalphacOf34(alphaViscosity,ADAFrInner,ADAFTemperatureElectron,Dotm)
        else:
            return 0


## Total accretion flow

def ADAFSpectralLuminosity(nu,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,UseApproxxM=False):
    '''This spectral luminosity is the sum of the bremsstrahlung-luminosity, the comptonised luminosity and the cyclosynchrotron-luminosity. It is in units of W/Hz, too.'''
    return ADAFSpectralLuminosityBremsstrahlungOf30(nu,alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)+ADAFSpectralLuminosityCyclosynchrotron(nu,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,UseApproxxM)+ADAFSpectralLuminosityComptonised(nu,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,UseApproxxM)

def ADAFSurface(ADAFrOuter,ADAFrInner,M9):
    '''This gives the surface of the accretion disk according to my geometry in units of m^2. Remark however, that 1997ApJ...477..585M assume a spherical accretion flow.'''
    return 2*np.pi*rSchwarzschildTom(ADAFrOuter,M9)**2 + 2*np.pi*rSchwarzschildTom(ADAFrInner,M9)**2 + 2*np.pi*(rSchwarzschildTom(ADAFrOuter,M9)*np.sqrt(1.25*rSchwarzschildTom(ADAFrOuter,M9)**2) - rSchwarzschildTom(ADAFrInner,M9)*np.sqrt(1.25*rSchwarzschildTom(ADAFrInner,M9)**2))

def ADAFEnergyDensitySpectralInFrequency(nu,SpectralLuminosity,*PositionalArgumentsOfSpectralLuminosityExceptnu):
    '''Assume a disk that has the spectral luminosity SpectralLuminosity (in W/Hz) and the outer and inner borders ADAFrOuter and ADAFrInner with the above geometry. Then this gives the spectral energy-density of photons (at the disk's borders/location) in units of J/(m^3*Hz). Hence, it is still spectral in respect to the frequency. nu is the photon's frequency in Hz.
    *PositionalArgumentsOfSpectralLuminosityExceptnu should be in the order alphaViscosity, betaPressure, ADAFrInner, ADAFrOuter, ADAFTemperatureElectron, M9, Dotm, depending on the arguments that are taken by the function SpectralLuminosity.'''
    return SpectralLuminosity(nu,*PositionalArgumentsOfSpectralLuminosityExceptnu)/(c*ADAFSurface(ADAFrOuter,ADAFrInner,M9))

def ADAFEnergyDensitySpectralInEnergy(epsilon,SpectralLuminosity,*PositionalArgumentsOfSpectralLuminosityExceptnu):
    '''The same as ADAFEnergyDensitySpectralInFrequency, however spectral in respect to and dependent on the energy epsilon of the emitted radiation. Hence, this is in 1/(m^3). epsilon is measured in J.'''
    return ADAFEnergyDensitySpectralInFrequency(epsilon/h,SpectralLuminosity,*PositionalArgumentsOfSpectralLuminosityExceptnu)/h

def ADAFNumberDensitySpectralInEnergy(epsilon,SpectralLuminosity,*PositionalArgumentsOfSpectralLuminosityExceptnu):
    '''This gives the spectral number-density of the ADAF-radiation from the disk at the disk's borders/location. It is in 1/(m^3*J).'''
    return ADAFEnergyDensitySpectralInEnergy(epsilon,SpectralLuminosity,*PositionalArgumentsOfSpectralLuminosityExceptnu)/epsilon

def ADAFNumberDensitySpectralInx(x,SpectralLuminosity,*PositionalArgumentsOfSpectralLuminosityExceptnu):
    '''This gives the spectral number-density of the ADAF-radiation from the disk at the disk's borders/location. It is however spectral in respect to and dependent on the dimensionless energy x. It is in units of 1/(m^3).'''
    return ADAFNumberDensitySpectralInEnergy(x*me*c**2,SpectralLuminosity,*PositionalArgumentsOfSpectralLuminosityExceptnu)*me*c**2

def ADAFNumberDensitySpectralInFrequency(nu,SpectralLuminosity,*PositionalArgumentsOfSpectralLuminosityExceptnu):
    '''This gives the spectral number-density of the ADAF-radiation from the disk at the disk's borders/location. It is however spectral in respect to and dependent on the frequency nu. It is in units of 1/(m^3*Hz).'''
    return ADAFEnergyDensitySpectralInFrequency(nu,SpectralLuminosity,*PositionalArgumentsOfSpectralLuminosityExceptnu)/(h*nu)

def ADAFTotalEnergyDensity(alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,Method):
    '''The total energy-density of the ADAF-radiation from the disk at the disk's borders/location. It is in J/m^3. It is yielded by integration of ADAFEnergyDensitySpectralInFrequency over the frequency nu. Method specifies whether the integration is to be performed 'Analytically' or 'Numerically'.'''
    # At performing the integration, it is reasonable to integrate the bremsstrahlung-, the cyclosynchrotron- and the comptonised contribution separately and to sum up the contributions afterwards.
    ADAFnuMin = ApplyADAFIntelligentStorage('ADAFnuMin') # Here, there is a kink. Left of it ADAFSpectralLuminosityCyclosynchrotronBelownuMin is valid, while ADAFSpectralLuminosityCyclosynchrotronOf25 is valid right of it.
    ADAFnuPeak = ApplyADAFIntelligentStorage('ADAFnuPeak') # Here, there is again a kink. Left of it ADAFSpectralLuminosityCyclosynchrotronOf25 is valid, while ADAFSpectralLuminosityComptonised is valid right of it.
    if Method == 'Analytically':
        BremsstrahlungContribution = ADAFEnergyDensitySpectralInFrequency(0,ADAFSpectralLuminosityBremsstrahlungOf30,alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)*kB*ADAFTemperatureElectron/h # This is the analytical integration of ADAFSpectralLuminosityBremsstrahlungOf30 from 0 to infinity.
        CyclosynchrotronContributionBelownuMin = ADAFEnergyDensitySpectralInFrequency(ADAFnuMin,ADAFSpectralLuminosityCyclosynchrotronBelownuMin,alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)*ADAFnuMin*13.0/35.0 # This is the analytical integration of ADAFSpectralLuminosityCyclosynchrotron from 0 to nu_min, where there is a kink.
        CyclosynchrotronContributionAbovenuMin = ADAFEnergyDensitySpectralInFrequency(ADAFnuPeak,ADAFSpectralLuminosityCyclosynchrotronOf25,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm)*ADAFnuPeak*5.0/7.0-ADAFEnergyDensitySpectralInFrequency(ADAFnuMin,ADAFSpectralLuminosityCyclosynchrotronOf25,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm)*ADAFnuMin*5.0/7.0 # This is the analytical integration of ADAFSpectralLuminosityCyclosynchrotron from nu_min, to nu_peak.
        ComptonisedContribution = ADAFEnergyDensitySpectralInFrequency(3*kB*ADAFTemperatureElectron/h,ADAFSpectralLuminosityComptonised,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)*(3*kB*ADAFTemperatureElectron/h)/(1-ADAFalphacOf34(alphaViscosity,ADAFrInner,ADAFTemperatureElectron,Dotm))-ADAFEnergyDensitySpectralInFrequency('ADAFSpectralLuminosityCyclosynchrotronPeak',ApplyADAFIntelligentStorage)*(ADAFnuPeak)/(1-ADAFalphacOf34(alphaViscosity,ADAFrInner,ADAFTemperatureElectron,Dotm)) # Analytical integration of ADAFSpectralLuminosityComptonised from nu_peak (where it adopts the value ApplyADAFIntelligentStorage('ADAFSpectralLuminosityCyclosynchrotronPeak')) up to 3*kB*ADAFTemperatureElectron/h.
    elif Method == 'Numerically':
        BremsstrahlungContribution = integrate.quad(ADAFEnergyDensitySpectralInFrequency, 0, 300*kB*ADAFTemperatureElectron/h, args=(ADAFSpectralLuminosityBremsstrahlungOf30,alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0] # This is the numerical integration of ADAFSpectralLuminosityBremsstrahlungOf30 from 0 to infinity. The upper integration border should actually be infinity. However, when it is too high a number, then the integration fails. Therefore the upper border was reduced to a value, that was found such that the result is very similar to the analytical one.
        CyclosynchrotronContributionBelownuMin = integrate.quad(ADAFEnergyDensitySpectralInFrequency, 0, ADAFnuMin, args=(ADAFSpectralLuminosityCyclosynchrotronBelownuMin,alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0] # Numerical integration of ADAFSpectralLuminosityCyclosynchrotron from 0 to nu_min, where there is a kink.
        CyclosynchrotronContributionAbovenuMin = integrate.quad(ADAFEnergyDensitySpectralInFrequency, ADAFnuMin, ADAFnuPeak, args=(ADAFSpectralLuminosityCyclosynchrotronOf25,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm))[0] # Numerical integration of ADAFSpectralLuminosityCyclosynchrotron from nu_min to nu_peak.
        ComptonisedLowestBorder = ADAFnuPeak # For the numerical integration of ADAFSpectralLuminosityComptonised one has to integrate from nu_peak up to 3*kB*ADAFTemperatureElectron/h. However, this range has to be subdivided into smaller ranges, for the integration not to fail.
        ComptonisedHighestBorder = 3*kB*ADAFTemperatureElectron/h
        ComptonisedNumberOfBorders = max(4,round(np.log10(ComptonisedHighestBorder)-np.log10(ComptonisedLowestBorder))) # ...Make at least 4 borders, that is 3 intervals.  ...
        ComptonisedListOfBorders = np.logspace(np.log10(ComptonisedLowestBorder),np.log10(ComptonisedHighestBorder),ComptonisedNumberOfBorders)
        ComptonisedContribution = 0
        for LeftBorder, RightBorder in zip(ComptonisedListOfBorders[:-1], ComptonisedListOfBorders[1:]):
            ComptonisedContribution += integrate.quad(ADAFEnergyDensitySpectralInFrequency, LeftBorder, RightBorder, args=(ADAFSpectralLuminosityComptonised,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0] # Numerical integration of ADAFSpectralLuminosityComptonised for each part-range.
    SumOfContributions = CyclosynchrotronContributionBelownuMin+CyclosynchrotronContributionAbovenuMin+ComptonisedContribution+BremsstrahlungContribution
    if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
        print(Method, 'determined total energy-density of ADAF-photons:          ', SumOfContributions, 'J/m^3')
        print('   ', Method, 'determined bremsstrahlung-contribution:               ', BremsstrahlungContribution, 'J/m^3')
        print('   ', Method, 'determined cyclosynchrotron-contribution below nu_min:', CyclosynchrotronContributionBelownuMin, 'J/m^3')
        print('   ', Method, 'determined cyclosynchrotron-contribution above nu_min:', CyclosynchrotronContributionAbovenuMin, 'J/m^3')
        print('   ', Method, 'determined comptonised contribution:                  ', ComptonisedContribution, 'J/m^3')
    return SumOfContributions

# Determine the total energy-density of the ADAF in units of J/m^3.
ADAFTotalEnergyDensityNum = ADAFTotalEnergyDensity(alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,'Numerically')
ADAFTotalEnergyDensityAna = ADAFTotalEnergyDensity(alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,'Analytically')

def ADAFTotalNumberDensity(alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,Method):
    '''The total number-density of the ADAF-radiation from the disk at the disk's borders/location. It is in 1/m^3. It is yielded by integration of ADAFNumberDensitySpectralInFrequency over the frequency nu. Method specifies whether the integration is to be performed 'Analytically' or 'Numerically'.'''
    # At performing the integration, it is reasonable to integrate the bremsstrahlung-, the cyclosynchrotron- and the comptonised contribution separately and to sum up the contributions afterwards.
    ADAFnuMin = ApplyADAFIntelligentStorage('ADAFnuMin') # Here, there is a kink. Left of it ADAFSpectralLuminosityCyclosynchrotronBelownuMin is valid, while ADAFSpectralLuminosityCyclosynchrotronOf25 is valid right of it.
    ADAFnuPeak = ApplyADAFIntelligentStorage('ADAFnuPeak') # Here, there is again a kink. Left of it ADAFSpectralLuminosityCyclosynchrotronOf25 is valid, while ADAFSpectralLuminosityComptonised is valid right of it.
    if Method == 'Analytically':
        BremsstrahlungWithoutExpTerm = ADAFSpectralLuminosityBremsstrahlungOf30Coefficient(alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)/(c*ADAFSurface(ADAFrOuter,ADAFrInner,M9)*h) # This is the number-density (spectral in frequency) without the term np.exp(-h*nu/(kB*ADAFTemperatureElectron))/nu. This term has to be integrated along nu from 0 to infinity. Fortunately, the error is not big, if it is integrated only from 1 to infinity.
        BremsstrahlungContribution = BremsstrahlungWithoutExpTerm*ExpIntegralFunction(h/(kB*ADAFTemperatureElectron)) # Integration of np.exp(-h*nu/(kB*ADAFTemperatureElectron))/nu from 0 to infinity is realised via the exponential integral.
        CyclosynchrotronContributionBelownuMin = ADAFEnergyDensitySpectralInFrequency(ADAFnuMin,ADAFSpectralLuminosityCyclosynchrotronBelownuMin,alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)*13.0/22.0/h # This is the analytical integration of ADAFSpectralLuminosityCyclosynchrotron from 0 to nu_min, where there is a kink.
        CyclosynchrotronContributionAbovenuMin = ADAFEnergyDensitySpectralInFrequency(ADAFnuPeak,ADAFSpectralLuminosityCyclosynchrotronOf25,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm)*5.0/(2.0*h)-ADAFEnergyDensitySpectralInFrequency(ADAFnuMin,ADAFSpectralLuminosityCyclosynchrotronOf25,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm)*5.0/(2.0*h) # This is the analytical integration of ADAFSpectralLuminosityCyclosynchrotron from nu_min, to nu_peak.
        ComptonisedContribution = ADAFEnergyDensitySpectralInFrequency(3*kB*ADAFTemperatureElectron/h,ADAFSpectralLuminosityComptonised,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)/(-h*ADAFalphacOf34(alphaViscosity,ADAFrInner,ADAFTemperatureElectron,Dotm))-ADAFEnergyDensitySpectralInFrequency('ADAFSpectralLuminosityCyclosynchrotronPeak',ApplyADAFIntelligentStorage)/(-h*ADAFalphacOf34(alphaViscosity,ADAFrInner,ADAFTemperatureElectron,Dotm)) # Analytical integration of ADAFSpectralLuminosityComptonised from nu_peak (where it adopts the value ApplyADAFIntelligentStorage('ADAFSpectralLuminosityCyclosynchrotronPeak')) up to 3*kB*ADAFTemperatureElectron/h.
    elif Method == 'Numerically':
        BremsstrahlungLowestBorder = 1#me**2*c**4/(300*kB*ADAFTemperatureElectron*h) # For the numerical integration of ADAFSpectralLuminosityBremsstrahlung one has to integrate from 0 up to infinity. However, this range has to be subdivided into smaller ranges, for the integration not to fail. As 1 is used as lowest border in the analytical computation, it is used here, too.
        BremsstrahlungHighestBorder = 300*kB*ADAFTemperatureElectron/h # The upper integration border should actually be infinity. However, when it is too high a number, then the integration fails.
        BremsstrahlungNumberOfBorders = max(5,round(np.log10(BremsstrahlungHighestBorder)-np.log10(BremsstrahlungLowestBorder))) # ...Make at least 5 borders, that is 4 intervals.  ...
        BremsstrahlungListOfBorders = np.logspace(np.log10(BremsstrahlungLowestBorder),np.log10(BremsstrahlungHighestBorder),BremsstrahlungNumberOfBorders)
        BremsstrahlungContribution = 0
        for LeftBorder, RightBorder in zip(BremsstrahlungListOfBorders[:-1], BremsstrahlungListOfBorders[1:]):
            BremsstrahlungContribution += integrate.quad(ADAFNumberDensitySpectralInFrequency, LeftBorder, RightBorder, args=(ADAFSpectralLuminosityBremsstrahlungOf30,alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0] # This is the numerical integration of ADAFSpectralLuminosityBremsstrahlungOf30. Therefore the upper border was reduced to a value, that was found such that the result is very similar to the analytical one.
        CyclosynchrotronContributionBelownuMin = integrate.quad(ADAFNumberDensitySpectralInFrequency, 0, ADAFnuMin, args=(ADAFSpectralLuminosityCyclosynchrotronBelownuMin,alphaViscosity,betaPressure,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0] # Numerical integration of ADAFSpectralLuminosityCyclosynchrotron from 0 to nu_min, where there is a kink.
        CyclosynchrotronContributionAbovenuMin = integrate.quad(ADAFNumberDensitySpectralInFrequency, ADAFnuMin, ADAFnuPeak, args=(ADAFSpectralLuminosityCyclosynchrotronOf25,alphaViscosity,betaPressure,ADAFTemperatureElectron,M9,Dotm))[0] # Numerical integration of ADAFSpectralLuminosityCyclosynchrotron from nu_min to nu_peak.
        ComptonisedLowestBorder = ADAFnuPeak#me**2*c**4/(300*kB*ADAFTemperatureElectron*h) # For the numerical integration of ADAFSpectralLuminosityComptonised one has to integrate from nu_peak up to 3*kB*ADAFTemperatureElectron/h. However, this range has to be subdivided into smaller ranges, for the integration not to fail.
        ComptonisedHighestBorder = 3*kB*ADAFTemperatureElectron/h
        ComptonisedNumberOfBorders = max(4,round(np.log10(ComptonisedHighestBorder)-np.log10(ComptonisedLowestBorder))) # ...Make at least 4 borders, that is 3 intervals.  ...
        ComptonisedListOfBorders = np.logspace(np.log10(ComptonisedLowestBorder),np.log10(ComptonisedHighestBorder),ComptonisedNumberOfBorders)
        ComptonisedContribution = 0
        for LeftBorder, RightBorder in zip(ComptonisedListOfBorders[:-1], ComptonisedListOfBorders[1:]):
            ComptonisedContribution += integrate.quad(ADAFNumberDensitySpectralInFrequency, LeftBorder, RightBorder, args=(ADAFSpectralLuminosityComptonised,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0] # Numerical integration of ADAFSpectralLuminosityComptonised for each part-range.
        # ComptonisedContribution = integrate.quad(ADAFNumberDensitySpectralInFrequency, ComptonisedLowestBorder, ComptonisedHighestBorder, args=(ADAFSpectralLuminosityComptonised,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm))[0]
    SumOfContributions = ComptonisedContribution+BremsstrahlungContribution#+CyclosynchrotronContributionBelownuMin+CyclosynchrotronContributionAbovenuMin
    if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
        print(Method, 'determined total number-density of ADAF-photons:          ', SumOfContributions, '1/m^3')
        print('   ', Method, 'determined bremsstrahlung-contribution:               ', BremsstrahlungContribution, '1/m^3')
        print('   ', Method, 'determined cyclosynchrotron-contribution below nu_min:', CyclosynchrotronContributionBelownuMin, '1/m^3')
        print('   ', Method, 'determined cyclosynchrotron-contribution above nu_min:', CyclosynchrotronContributionAbovenuMin, '1/m^3')
        print('   ', Method, 'determined comptonised contribution:                  ', ComptonisedContribution, '1/m^3')
    return SumOfContributions
    # Sometimes, it is assumed that ADAF-photons pair-produce electrons with themselves. Then, it is useful to determine the total number-density of ADAF-photons above the PP-threshold. The maximum energy of the ADAF-spectrum is 300*kB*ADAFTemperatureElectron. Hence, the threshold is me**2*c**4/(300*kB*ADAFTemperatureElectron). Thus, one has to assign me**2*c**4/(300*kB*ADAFTemperatureElectron*h) to ComptonisedLowestBorder and to BremsstrahlungLowestBorder. The contributions CyclosynchrotronContributionBelownuMin and CyclosynchrotronContributionAbovenuMin are then =0 and have to be removed from SumOfContributions.

# Determine the total number-density of the ADAF in units of 1/m^3.
ADAFTotalNumberDensityNum = ADAFTotalNumberDensity(alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,'Numerically')
ADAFTotalNumberDensityAna = ADAFTotalNumberDensity(alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,'Analytically')


## Evaluation

if "MainProcess" in multiprocessing.current_process().name: # This prevents evaluation in child-processes.
    
    def EvaluateADAFSpectralLuminosityBremsstrahlungOf30():
        ValuesFornu = np.logspace(8,23,100)
        ValuesForADAFSpectralLuminosityBremsstrahlungOf30 = ADAFSpectralLuminosityBremsstrahlungOf30(ValuesFornu,alphaViscosity,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm)
        pl.figure(figsize=(12, 9), num="Spectral bremsstrahlung-luminosity of an ADAF versus frequency")
        pl.title('Spectral bremsstrahlung-luminosity of an ADAF versus frequency', fontsize=24)
        pl.loglog(ValuesFornu,ValuesFornu*ValuesForADAFSpectralLuminosityBremsstrahlungOf30, label='$\dot m = %g$' % Dotm)
        pl.legend(loc="best", fontsize=20)
        pl.xlabel('frequency $f / {\mathrm{Hz}}$', fontsize=24)
        pl.ylabel('frequency times spectral luminosity $f \cdot L(f) / \mathrm{W}$', fontsize=24)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=24)
        ax.yaxis.set_tick_params(labelsize=24)
        pl.ylim(1.0e28,1.0e36)
    
    def EvaluateADAFSpectralLuminosityCyclosynchrotron():
        ValuesFornu = np.logspace(6,15,1000)
        ValuesForADAFSpectralLuminosityCyclosynchrotron = np.asarray([ADAFSpectralLuminosityCyclosynchrotron(i,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,UseApproxxM) for i in ValuesFornu])
        pl.figure(figsize=(12, 9), num="Spectral cyclosynchrotron-luminosity of an ADAF versus frequency")
        pl.title('Spectral cyclosynchrotron-luminosity of an ADAF versus frequency', fontsize=24)
        pl.loglog(ValuesFornu,ValuesFornu*ValuesForADAFSpectralLuminosityCyclosynchrotron, label='$\dot m = %g$' % Dotm)
        pl.legend(loc="best", fontsize=20)
        pl.xlabel('frequency $f / {\mathrm{Hz}}$', fontsize=24)
        pl.ylabel('frequency times spectral luminosity $f \cdot L(f) / \mathrm{W}$', fontsize=24)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=24)
        ax.yaxis.set_tick_params(labelsize=24)
        pl.ylim(1.0e28,1.0e36)
    
    def EvaluateADAFSpectralLuminosityComptonised():
        ValuesFornu = np.logspace(9,23,1000)
        ValuesForADAFSpectralLuminosityComptonised = np.asarray([ADAFSpectralLuminosityComptonised(i,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm,UseApproxxM) for i in ValuesFornu])
        pl.figure(figsize=(12, 9), num="Spectral comptonised luminosity of an ADAF versus frequency")
        pl.title('Spectral comptonised-luminosity of an ADAF versus frequency', fontsize=24)
        pl.loglog(ValuesFornu,ValuesForADAFSpectralLuminosityComptonised, label='$\dot m = %g$' % Dotm)
        pl.legend(loc="best", fontsize=20)
        pl.xlabel('frequency $f / {\mathrm{Hz}}$', fontsize=24)
        pl.ylabel('frequency times spectral luminosity $f \cdot L(f) / \mathrm{W}$', fontsize=24)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=24)
        ax.yaxis.set_tick_params(labelsize=24)
        #pl.ylim(1.0e28,1.0e36)
    
    def EvaluateADAFSpectralLuminosityForComparison():
        '''I could approximately reproduce figure 1 of Mahadevan using an electron temperature of 5 cdot 10^9 K. Furthermore, I could very approximately reproduce figure 2 of Esin using an electron temperature of 3.5 cdot 10^9 K and setting ADAFrOuter equal to Esin's r_tr.'''
        ValuesFornu = np.logspace(6,22,2000)
        ValuesDotm=[0.0001,0.0002,0.0004,0.0008]
        ValuesT=[[1.7*10**10,1.3*10**10],[1.4*10**10,1.0*10**10],[1.1*10**10,8.1*10**9],[8.4*10**9,6.4*10**9]]
        for i in [0,1,2,3]:
            for j in [0,1]:
                Dotm=ValuesDotm[i]
                ADAFTemperatureElectron=ValuesT[i][j]
                if j==0:
                    Linestyle='-'
                else:
                    Linestyle='--'
                if i==0:
                    Color='blue'
                elif i==1:
                    Color='red'
                elif i==2:
                    Color='green'
                elif i==3:
                    Color='orange'
                ValuesForADAFSpectralLuminosity = np.asarray([ADAFSpectralLuminosity(i,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm) for i in ValuesFornu])
                pl.figure(figsize=(12, 9), num="Spectral luminosity of an ADAF versus frequency")
                #pl.title('Spectral luminosity of an ADAF versus frequency', fontsize=24)
                pl.loglog(ValuesFornu,ValuesFornu*ValuesForADAFSpectralLuminosity, label='$\dot m = %g$' % Dotm, color=Color, linestyle=Linestyle)
        ValuesForADAFSpectralLuminosity = np.asarray([ADAFSpectralLuminosity(i,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm) for i in ValuesFornu])
        pl.figure(figsize=(12, 9), num="Spectral luminosity of an ADAF versus frequency")
        #pl.title('Spectral luminosity of an ADAF versus frequency', fontsize=24)
        pl.loglog(ValuesFornu,ValuesFornu*ValuesForADAFSpectralLuminosity, color='black')
        #pl.legend(loc="best", fontsize=20)
        pl.xlabel('frequency $\\nu / {\mathrm{Hz}}$', fontsize=24)
        pl.ylabel('freq. times spec. lum. $\\nu \cdot L_{\mathrm{ADAF}}(\\nu) / \mathrm{W}$', fontsize=24)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=24, pad=8)
        ax.yaxis.set_tick_params(labelsize=24)
        pl.xlim(1.0e8,1.0e21)
        pl.ylim(1.0e26,1.0e34)
        pl.subplots_adjust(bottom=0.13)
    
    def EvaluateADAFNumberDensitySpectralInx():
        ValuesForx = np.logspace(-19,3,20000)
        ValuesForADAFNumberDensitySpectralInx = np.asarray([ADAFNumberDensitySpectralInx(i,ADAFSpectralLuminosity,alphaViscosity,betaPressure,ADAFrInner,ADAFrOuter,ADAFTemperatureElectron,M9,Dotm) for i in ValuesForx])
        pl.figure(figsize=(16, 12), num="Spectral number-density of an ADAF versus dimensionless energy")
        pl.title('Spectral number-density of an ADAF versus dimensionless energy', fontsize=24)
        pl.loglog(ValuesForx,ValuesForADAFNumberDensitySpectralInx, label='Total')
        pl.legend(loc="best", fontsize=20)
        pl.xlabel('Dimensionless energy $x$', fontsize=24)
        pl.ylabel('Spectral number-density in $m^{-3}$', fontsize=24)
        ax = pl.gca()
        ax.xaxis.set_tick_params(labelsize=24)
        ax.yaxis.set_tick_params(labelsize=24)
        pl.ylim(1.0e-20,1.0e25)
