# Generalized-Gilat-Raubenheimer-method
A method for density-of-states calculation in band theory, especially in photonic crystals.

The DOS calculation is in 'DOS_GR(n).m' file for Matlab. The necessary input data is illustrated in the m file. We can also use MPB to calculate the band strucutre. 

The script file 'dompb.sh' will help you to run the mpb and extract the data in right format for Matlab file.

'DG.ctl' is an example of a double gyroid photonic crsytal sturcture. In the ctl file it gives required parameters to describe the structure. We can also change the sampling set in ctl files.

Important notice: parameters of ctl file and Matlab file should agree with each other. For example, the reciprocal lattice vectors, the number of k points in each dimension in mesh, the inter k-points number between two high symmetry points in line, etc.
