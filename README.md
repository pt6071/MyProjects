# MyProjects

A potpourri of various projects of profesional and personal interest, mostly acoustics, beamforming, and DSP-related.  Note that all Matlab code has been developed to run on version R2017a and that all C code has been developed for Microsoft Visual Studio 2012.

### Acoustics
* **scatterPlot.m** - Matlab script that computes and plots the magnitude of the pressure field resulting from either a perfectly rigid sphere or cylinder scattering a plane wave propagating up the y-axis.
* **cylindricalSolver.m** - Matlab function that performs the number crunching for **cylindricalScatter.m**.
* **sphericalSolver.m** - Matlab function that performs the number crunching for **sphericalScatter.m**.

### Beamforming
* **cardioidLinearArray.m** - Matlab script for calculating the beampattern of a uniform linear array (ULA) of directional cardioid-type elements.  Summing cardioids into a ULA is a great combination since cardioids have good low-frequency directivity but suffer SNR loss while linear summing arrays have good high-frequency directivity, improve SNR, and can be electronically steered to listen in different directions.

### DSP
* **myShelf.m** - Matlab function for calculating the taps for a digital shelf filter.
* **blms.m** - Matlab implementation of a Block LMS adaptive filter.
* **blmsdemo.m** - Matlab demonstration script of a BLMD filter identifying a system.
