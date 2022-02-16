# PhD-project
This is the python code developed and used during my PhD time.

The code consists out of seven python files and several auxiliary files, as outlined in the dissertation.

Subject of the code is:
- part1.py: Initialisation of the code, basic physical definitions, basic auxiliary functions, treatment of ADAF spectra
- part2.py: A toy model of AGN jets (not discussed in the dissertation)
- part3.py: A toy model of Mrk 501 (not discussed in the dissertation)
- part4.py: Numerical solution of steady-state IC pair cascades: Specification of input parameters, definition of interaction rates and initialisation of the iterative scheme
- part5.py: Numerical solution of steady-state IC pair cascades: Definitions of the terms of the electron kinetic equation, iterative scheme, plot of N, export of input parameters
- part6.py: Numerical solution of steady-state IC pair cascades: Experimental extension to synchrotron IC pair cascades, specification of the synchrotron photons kinetic equation, super iteration
- part7.py: Numerical solution of steady-state IC pair cascades: Obsolete part to determine N in the Thomson regime, computation of the HE photon spectrum, conversion to observed flux density, comparison with observational data, main file executions

The following auxiliary data files are provided:
- Four .csv files that contain the used emission lines and their relative strengths, three files for the applications to 3C 279 and one file for the Mrk 501 case.
- The file "Mrk 501 observational data.txt" containing the used MAGIC data points, as described in chapter 5 and as shown in figure 5.5.
- The file "SSC model for Mrk 501.dat" that contains the added SSC-modelled flux density used in chapter 5.
- The file "3C 279 observational data.txt" containing the used Fermi LAT data points, as described in chapter 6 and as shown in figure 6.4.
- Seven .dat files that include the theoretically determined SEDs in the case of 3C 279, as described in chapter 6 and as shown in figure 6.4.

The code was developed with python 3.4.1 with the Interactive Editor for Python (IEP) 3.6 of Pyzo. It is optimised for usage on 64 bit Windows systems, but it might also run on other systems. Multi-core applications might not be working within IEP. To execute the multi-core-taylored algorithms, Windows shell, power shell or shell on unix systems has to be used. 
