# Density-of-states calculation in photonic crystals


## Authorship

Written by Boyuan Liu (zkdlby@mail.ustc.edu.cn) in [L01 group](http://l01.iphy.ac.cn/L01web-English/html/index-english.html) in Institute of Physics CAS, under instruction of [Prof. Ling Lu](http://l01.iphy.ac.cn/linglu/). 

## Brief Introduction


We present two open-source programs to calculate density of states (DOS) in photonic crystals using Generalized-Gilat-Raubenheimer (GGR) method and tetrahedron (Tr) method respectively. GGR method program needs both the frequency band data and group velocity data to calculate DOS and tetrahedron method program only needs band data while has a relative less accuracy compared with GGR. 

We suggest to use MIT Photonic-Bands (MPB) to calculate the frequency band and group velocity. Input these data files into the DOS calculation program 'DOS_GR.m' or 'DOS_Tr.m' to obtain the DOS of the structure. You can also directly input the data files from other band-computing softwares to the DOS program, on the condition that the data is arranged in the right format in the files.

## Citation


We propose the GGR method in our article ['Generalized Gilat-Raubenheimer method for density-of-states calculation in photonic crystals'](https://doi.org/10.1088/2040-8986/aaae52). If you use our code for your research, please cite it properly, like  "The density of states is calculated by the program in  \cite{Liu2017GGR}. "

Reference

    @article{10.1088/2040-8986/aaae52,
      author={Boyuan Liu and J D Joannopoulos and Steven G Johnson and Ling Lu},
        title={Generalized Gilat-Raubenheimer method for density-of-states calculation in photonic crystals},
      journal={Journal of Optics},

      url={https://doi.org/10.1088/2040-8986/aaae52},
      year={2018},
      abstract={Abstract Efficient numeric algorithm is the key for accurate evaluation of density of states (DOS) in band theory. Gilat-Raubenheimer (GR) method proposed in 1966 is an efficient linear extrapolation method which was limited in specific lattices. Here, using an affine transformation, we provide a new generalization of the original GR method to any Bravais lattices and show that it is superior to the tetrahedron method and the adaptive Gaussian broadening method. Finally, we apply our generalized GR (GGR) method to compute DOS of various gyroid photonic crystals of topological degeneracies.}
    }

## Usage


There are two ways to use our programs to calculate DOS.  The first one is to use MPB compute the band data and input them into the DOS calculation programs and the other one is to input the band data directly to the DOS calculation programs. Both ways require user to set the parameters of the photonic crystals correctly in the DOS calculation programs and to adjust the input band data files in the right format. We provide two DOS calculation programs, *DOS_GGR.m* and *DOS_GGR.m*, with different algorithms, and their input files have different requirement as well. The following is the guide of the two ways.

### Using MPB to obtain band data

Firstly, run script file *dompb.sh* in Linux system,

    ./dompb.sh
    
and it will run MPB with the structure-data file *3Dexample.ctl* in the same folder and process the band data and group velocity data into correct format.

Then put the data files and DOS calculation file *DOS_GGR.m* in the same directory and use Matlab to run the file *DOS_GGR.m*. It will output the DOS data and the corresponding plot as below (the structure figure is added afterwards).

![example](https://github.com/boyuanliuoptics/DOS-calculation/blob/master/example.png)

You can also change the settings in the the ctl file at the beginning so that the script file *dompb.sh* will output the data for *DOS_Tr.m*. Then run the file *DOS_Tr.m* with the band data files in the same directory. It will output the DOS data and the corresponding plot.

### Using other band-computing softwares to obtain band data

Firstly, make the input data files into the same format as the output files of *dompb.sh*. Then put the data files and DOS calculation file *DOS_GGR.m* in the same directory and use Matlab to run the file *DOS_GGR.m* or *DOS_Tr.m* according to the input files and the algorithm.

The input files are *band.txt*, *frequency_GGR.txt* and *velocity_GGR.txt* for *DOS_GGR.m*, or *band.txt* and *frequency_Tr.txt* for *DOS_Tr.m*. The right format of each file is as below.

#### *band.txt* (optional)

Storing the frequency band for drawing band structure in line.

        0.5 -0.5 0.5 0.418621 0.418708 0.418778 ... 0.799886 0.800063
        0.483871 -0.483871 0.483871 0.417886 0.417935 0.418915 ... 0.799816 0.800123
        0.467742 -0.467742 0.467742 0.415525 0.415591 0.419394 ... 0.799351 0.79984
        ...
        0 0 0 0 0 0.479697 ... 0.788944 0.798401

In one row, the first three numbers are the coordinate components of one point in k-space, whose units are the reciprocal vectors. The following numbers are band frequencies at this point in ascending order.

Position and band information of different k-points is in the same line. The values (including coordinates and band frequencies) in the same line should be divided by one space.

#### *frequency_GGR.txt*

Storing the frequency band for GGR DOS calculation.

        0.342301
        0.34426
        0.39272
        ...
        0.801733
        
The number in the one row is one band frequency of some k-point. Each frequency is divided by a '\n'. You only need to ensure that the frequency sequence is in accordance with that of *velocity_GGR.txt* (the data of the same line describes the information of the same k-point and of the same band).

#### *velocity_GGR.txt*

Storing the group velocity for GGR DOS calculation.

        0.03878239605480597 0.050501326087597484 0.05049987902906693
        0.05418725891376245 0.05776810445673396 0.05776682796316916
        0.005300927843952232 0.021269420720856042 0.021267892167480044
        ...
        0.02864868228729669 -0.030525708954528805 -0.02997708531049537

The numbers in one row are three orthogonal components of group velocity of some k-point in some band. Each group velocity (three numbers) is divided by a '\n'. You only need to ensure that the group velocity sequence is in accordance with that of *frequency_GGR.txt*.

#### *frequency_Tr.txt*

Storing the frequency band and position of sampling k points for Tr DOS calculation.

        0 -0.5 -0.5 0.336365 0.336982 0.390348 ... 0.809612 0.812481
        0.1 -0.5 -0.5 0.343023 0.343637 0.394022 ... 0.809883 0.81058
        0.2 -0.5 -0.5 0.36122 0.361826 0.402864 ... 0.807251 0.808129
        ...
        0.5 0.5 0.5 0.418717 0.418717 0.41889 ... 0.799788 0.799807

The data format is the same as that of *band.txt*.

## 2D cases

We expand GGR and Tr methods into in 2D structures and find that the DOS of 2D structures has bad smoothness (discontinuity in its first derivative) which do not appear in 3D. The reason for this zigzag DOS plot is that the formula of DOS contribution of each subcell in 2D cases is not continous in the first derivative.
