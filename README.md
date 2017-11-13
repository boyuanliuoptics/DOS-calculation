# Generalized-Gilat-Raubenheimer-method

1. Introduction-----------------------------------------------------------------------------

A method for density-of-states calculation in band theory, especially in photonic crystals.

The DOS calculation is in 'DOS_GR(n).m' file for Matlab. The necessary input data is illustrated in the '*.m' file. We can also use MPB to calculate the band strucutre. 

The script file 'dompb.sh' will help you to run the mpb and extract the data in right format for Matlab file.

'DG.ctl' is an example of a double gyroid photonic crsytal sturcture. In the ctl file it gives required parameters to describe the structure. We can also change the sampling set in ctl files.

Important notice: parameters of '*.ctl' file and Matlab file should agree with each other. For example, the reciprocal lattice vectors, the number of k points in each dimension in mesh, the inter k-points number between two high symmetry points in line, etc.


2. Learn about mpb------------------------------------------------------------------------

You can download mpb from the website: https://mpb.readthedocs.io/en/latest/Download/. Some introduction documents for mpb are also in the website.

You need to change some presets to use some convenient usage like multithread or running on a sever.


3. Details in the programs of GGR method and tetrahedron (Tr) method-----------------------

In general, GGR method is an extrapolation method and Tr method is an interpolation method. The details could refer to our article. As for the coding aspect, there are some notable details. 

