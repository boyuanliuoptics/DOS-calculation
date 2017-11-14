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

References
----------------------------------------
[1] Andrew J Morris, Rebecca J Nicholls, Chris J Pickard, and Jonathan R Yates. Optados: A tool for obtaining density of states, core-level and optical spectra from electronic structure codes. Computer Physics Communications, 185(5):1477–1485, 2014. 

[2] G. Gilat and L. J. Raubenheimer. Accurate numerical method for calculating frequency-distribution functions in solids. Physical Review, 144(2):390–395, 1966. 

[3] L. J Raubenheimer and G Gilat. Accurate numerical method of calculating frequency distribution functions in solids. ii. extension to hcp crystals. Physical Review, 157(3):586–599, 1967.

[4] Z Kam and G Gilat. Accurate numerical method for calculating frequency distribution functions in solids. iii. extension to tetragonal crystals. Physical Review, 175(3):1156, 1968. 

[5] E Finkman, Z Kam, E Cohen, and G Gilat. Accurate numerical method for calculating spectra in solidsiv. extension to trigonal crystals. Journal of Physics and Chemistry of Solids, 32(10):2423–2427, 1971. 

[6] H Bross. On the eﬃciency of diﬀerent schemes for the evaluation of the density of states and related properties in solids. physica status solidi (b), 179(2):429–439, 1993. 

[7] Jonathan R Yates, Xinjie Wang, David Vanderbilt, and Ivo Souza. Spectral and fermi surface properties from wannier interpolation. Physical Review B, 75(19):195121, 2007. 

[8] CJ Pickard and MC Payne. Extrapolative approaches to brillouin-zone integration. Physical Review B, 59(7):4685, 1999. 

[9] CJ Pickard and MC Payne. Second-order k p perturbation theory with vanderbilt pseudopotentials and plane waves. Physical Review B, 62(7):4383, 2000. 

[10] G Lehmann and M Taut. On the numerical calculation of the density of states and related properties. physica status solidi (b), 54(2):469–477, 1972. 

[11] O Jepson and OK Anderson. The electronic structure of hcp ytterbium. Solid State Communications, 9(20):1763–1767, 1971. 

[12] Peter E Bl¨ochl, Ove Jepsen, and Ole Krogh Andersen. Improved tetrahedron method for brillouin-zone integrations. Physical Review B, 49(23):16223, 1994. 

[13] MS Methfessel, MH Boon, and FM Mueller. Analyticquadratic method of calculating the density of states. Journal of Physics C: Solid State Physics, 16(27):L949, 1983. 

[14] MH Boon, MS Methfessel, and FM Mueller. Singular integrals over the brillouin zone: the analytic-quadratic method for the density of states. Journal of Physics C: Solid State Physics, 19(27):5337, 1986. 

[15] MS Methfessel, MH Boon, and FM Mueller. Singular integrals over the brillouin zone: inclusion of k-dependent matrix elements. Journal of Physics C: Solid State Physics, 20(8):1069, 1987. 

[16] Jorge E Mu¨ller and John W Wilkins. Band-structure approach to the x-ray spectra of metals. Physical Review B, 29(8):4331, 1984.

[17] Kurt Busch and Sajeev John. Photonic band gap formation in certain self-organizing systems. Physical Review E, 58(3):3896, 1998. 

[18] Patrick M Johnson, A Femius Koenderink, and Willem L Vos. Ultrafast switching of photonic density of states in photonic crystals. Physical Review B, 66(8):081102, 2002. 

[19] Ivan S Nikolaev, Willem L Vos, and A Femius Koenderink. Accurate calculation of the local density of optical states in inverse-opal photonic crystals. JOSA B, 26(5):987– 997, 2009. 

[20] P Kano, D Barker, and M Brio. Analysis of the analytic dispersion relation and density of states of a selected photonic crystal. Journal of Physics D: Applied Physics, 41(18):185106, 2008. 

[21] Victor Liu and Shanhui Fan. Eﬃcient computation of equifrequency surfaces and density of states in photonic crystals using dirichlet-to-neumann maps. JOSA B, 28(8):1837–1843, 2011. 

[22] G Gilat. Analysis of methods for calculating spectral properties in solids. Journal of Computational Physics, 10(3):432–465, 1972.

[23] G Wiesenekker and EJ Baerends. Quadratic integration over the three-dimensional brillouin zone. Journal of Physics: Condensed Matter, 3(35):6721, 1991. 

[24] Ling Lu, Liang Fu, John D Joannopoulos, and Marin Soljaˇci´c. Weyl points and line nodes in gyroid photonic crystals. Nature photonics, 7(4):294–299, 2013. 

[25] Luyang Wang, Shao Kai Jian, and Hong Yao. Topological photonic crystal with equifrequency weyl points. Physics, (6), 2015. 

[26] Steven G. Johnson and J. D. Joannopoulos. Block-iterative frequency-domain methods for maxwell’s equations in a planewave basis. Opt. Express, 8(3):173–190, 2001. 

[27] R. I. Saye. High-order quadrature methods for implicitly deﬁned surfaces and volumes in hyperrectangles. Siam Journal on Scientiﬁc Computing, 37(2):A993–A1019, 2015.
