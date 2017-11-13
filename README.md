# Generalized-Gilat-Raubenheimer-method

1. -------------------------------------Introduction----------------------------------------

A method for density-of-states calculation in band theory, especially in photonic crystals.

The DOS calculation is in 'DOS_GR(n).m' file for Matlab. The necessary input data is illustrated in the '*.m' file. We can also use MPB to calculate the band strucutre. 

The script file 'dompb.sh' will help you to run the mpb and extract the data in right format for Matlab file.

'DG.ctl' is an example of a double gyroid photonic crsytal sturcture. In the ctl file it gives required parameters to describe the structure. We can also change the sampling set in ctl files.

Important notice: parameters of '*.ctl' file and Matlab file should agree with each other. For example, the reciprocal lattice vectors, the number of k points in each dimension in mesh, the inter k-points number between two high symmetry points in line, etc.


2. -------------------------------------Learn about mpb-----------------------------------

You can download mpb from the website: https://mpb.readthedocs.io/en/latest/Download/. Some introduction documents for mpb are also in the website.

You need to change some presets to use some convenient usage like multithread (mpb-split) or running on a sever:

(a) Multithread computation. Use 'which mpb-split' to find the path and edit the file 'mpb-split'. A part of the code is as below:

-----------------coding part starts------------------
i=1
while test `expr $i \< $1` = 1; do
    /usr/local/bin/mpb interactive?=false k-split-index=$i k-split-num=$* > $tmpname.$i &
    subprocesses="$subprocesses $!"
    i=`expr $i + 1`
done
------------------coding part ends-------------------

and change 'i=1' to 'i=0' in the first line.

(b) Running on a sever. Add an extra expression before running mpb: 'OMP_NUM_THREADS=1'. For example, in 'dompb.sh' line 41,

-----------------coding part starts------------------
OMP_NUM_THREADS=1 mpb-split 8 Zone?=false $file_band_ctl > $data_bandline
------------------coding part ends-------------------

and with this, one thread will be processed in one cpu core.

3. ----------------Details in the programs of GGR method and tetrahedron (Tr) method-------

In general, GGR method is an extrapolation method and Tr method is an interpolation method. The details could refer to our article. As for the coding aspect, there are some notable details. Please check the example 'DG.ctl' for detailed infomation.

(a) The k-points mesh of the two methods are different. Both of two should cover the Brillouin zone (or half of BZ) completely with the subcells, parallelepipeds or tetrahedra. 

(b) GGR method would perform better with a k-points mesh avoiding Gamma point, whose group velocities are zero in many cases. Therefore we suggest to use an even number of sampling k points along one dimension, when the range of this dimension is [-bi/2, bi/2].
