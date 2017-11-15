Generalized-Gilat-Raubenheimer-method
=====================================================

Authorship
----------------------------------------
Written by Boyuan Liu in Matlab, scheme and Linux shell (zkdlby@mail.ustc.edu.cn). 

Brief Introduction
----------------------------------------

We present an open-source program to calculate density of states (DOS) using Generalized-Gilat-Raubenheimer (GGR) method for band theory, especially in photonic crystals. You can use the code to calculate the DOS from the frequency band data and group velocity data. We also provide another DOS calculation program using tetrahedron method which has a relative less accuracy compared with GGR while it only needs band data without demand for group velocity data. 

We suggest to use MIT Phtonic-Bands (MPB) to calculate the frequency band and group velocity. Input these data into the DOS calculation program 'DOS_GR.m' or 'DOS_Tr.m' to obtain DOS of the structure. You can also directly input the data files from other band compuation softwares to the DOS program, as long as the data is arranged in the right format in the files.

Citation
----------------------------------------
We propose the GGR method in our article 'Generalized Gilat-Raubenheimer method for density-of-states calculation in photonic crystals'. If you use our code for your research, please cite it properly. It will be post on arxiv soon. To be continue

Usage
----------------------------------------

