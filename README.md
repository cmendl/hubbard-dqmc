Determinant quantum Monte Carlo implementation
==============================================

C implementation and test files of the determinant quantum Monte Carlo (DQMC) method applied to Hubbard-type models.

A *Makefile* is available in the *bin* subfolder, assuming that the Intel C compiler and MKL are available.

To run the code, call `hubbard_dqmc <paramfile>`; some example parameter files are provided in the *bin* subfolder.

The Mathematica unit test notebooks can be opened by [Mathematica](https://www.wolfram.com/mathematica) or the free [CDF player](https://www.wolfram.com/cdf-player).


References
----------
1. R. Blankenbecler, D. J. Scalapino, R. L. Sugar  
   Monte Carlo calculations of coupled boson-fermion systems. I  
   Phys. Rev. D 24, 2278 (1981) [DOI](https://doi.org/10.1103/PhysRevD.24.2278)
2. S. R. White, D. J. Scalapino, R. L. Sugar, E. Y. Loh, J. E. Gubernatis, and R. T. Scalettar  
   Numerical study of the two-dimensional Hubbard model  
   Phys. Rev. B 40, 506-516 (1989) [DOI](https://doi.org/10.1103/PhysRevB.40.506)
3. Z. Bai, C.-R. Lee, R.-C. Li, S. Xu  
   Stable solutions of linear systems involving long chain of matrix multiplications  
   Linear Algebra Appl. 435, 659-673 (2011) [DOI](https://doi.org/10.1016/j.laa.2010.06.023)
4. A. Tomas, C.-C. Chang, R. Scalettar, Z. Bai  
   Advancing large scale many-body QMC simulations on GPU accelerated multicore systems  
   IEEE 26th International Parallel & Distributed Processing Symposium (IPDPS) 308-319 (2012) [DOI](https://doi.org/10.1109/IPDPS.2012.37)
5. S. Gogolenko, Z. Bai, R. Scalettar  
   Structured orthogonal inversion of block p-cyclic matrices on multicore with GPU accelerators  
   Euro-Par 2014 Parallel Processing, LNCS 8632, pages 524-535 (2014) [DOI](https://doi.org/10.1007/978-3-319-09873-9_44)
6. C. Jiang, Z. Bai, R. Scalettar  
   A fast selected inversion algorithm for Green's function calculation in many-body quantum Monte Carlo simulations  
   IEEE International Parallel and Distributed Processing Symposium, 2016 [DOI](https://doi.org/10.1109/IPDPS.2016.69)
