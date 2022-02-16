# PhD-project
This is the python code developed and used during my PhD time.

The code consists out of seven python files, as outlined in the dissertation.
Subject of the code is:
1. Initialisation of the code, basic physical definitions, basic auxiliary functions, treatment of ADAF spectra
2. A toy model of AGN jets (not discussed in the dissertation)
3. A toy model of Mrk 501 (not discussed in the dissertation)
4. Numerical solution of steady-state IC pair cascades: Specification of input parameters, definition of interaction rates and initialisation of the iterative scheme
5. Numerical solution of steady-state IC pair cascades: Definitions of the terms of the electron kinetic equation, iterative scheme, plot of N, export of input parameters
6. Numerical solution of steady-state IC pair cascades: Experimental extension to synchrotron IC pair cascades, specification of the synchrotron photons kinetic equation, super iteration
7. Numerical solution of steady-state IC pair cascades: Obsolete part to determine N in the Thomson regime, computation of the HE photon spectrum, conversion to observed flux density, comparison with observational data, main file executions

The code was developed with python 3.4.1 with the Interactive Editor for Python (IEP) 3.6 of Pyzo. It is optimised for usage on 64 bit Windows systems, but it might also run on other systems. Multi-core applications might not be working within IEP. To execute the multi-core-taylored algorithms, Windows shell, power shell or shell on unix systems has to be used. 
