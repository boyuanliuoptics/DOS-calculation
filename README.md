Density-of-states calculation in photonic crystals
=====================================================

Authorship
----------------------------------------
Written by Boyuan Liu (zkdlby@mail.ustc.edu.cn) in [L01 group](http://l01.iphy.ac.cn/L01web-English/html/index-english.html) in Institute of Physics CAS leadered by [Prof. Ling Lu](http://l01.iphy.ac.cn/linglu/). 

Brief Introduction
----------------------------------------

We present two open-source programs to calculate density of states (DOS) in photonic crystals using Generalized-Gilat-Raubenheimer (GGR) method and tetrahedron (Tr) method respectively. GGR method program needs both the frequency band data and group velocity data to calculate DOS and tetrahedron method program only needs band data while has a relative less accuracy compared with GGR. 

We suggest to use MIT Phtonic-Bands (MPB) to calculate the frequency band and group velocity. Input these data files into the DOS calculation program 'DOS_GR.m' or 'DOS_Tr.m' to obtain the DOS of the structure. You can also directly input the data files from other band-compuating softwares to the DOS program, on the condition that the data is arranged in the right format in the files.

Citation
----------------------------------------

We propose the GGR method in our article 'Generalized Gilat-Raubenheimer method for density-of-states calculation in photonic crystals'. If you use our code for your research, please cite it properly, like  "The density of states is calculated by the program in  \cite{Liu2017GGR}. "

Reference

    @article{Liu2017GGR,
    archivePrefix = {arXiv},
    arxivId = {1711.07993},
    author = {Liu, Boyuan and Johnson, Steven G. and Joannopoulos, John D. and Lu, Ling},
    eprint = {1711.07993},
    title = {{Generalized Gilat-Raubenheimer method for density-of-states calculation in photonic crystals}},
    url = {http://arxiv.org/abs/1711.07993},
    year = {2017}
    }

Usage
----------------------------------------
Firstly, run script file *dompb.sh* in Linux system,

    ./dompb.sh
    
and it will run MPB with the structure-data file *3Dexample.ctl* in the same folder and process the band data and group velocity data into correct format.

Then use Matlab to run the file *DOS_GGR.m*. It will output the DOS data and the corresponding plot as below.



You can change the settings in the the ctl file to use *DOS_Tr.m* to calculate the DOS. It is ok to use band data from other band-computing softwares, as long as the band data is in the right form.

Notice and Future Plans
----------------------------------------
The two methods are simple to expand to 2D structures. However, the DOS in 2D structures is not as smooth as that in 3D due to the inherent characteristic of the two methods. Thus we recommend Gaussian method in 2D to obtain a better DOS (will be available soon).
