# HSS-matfun
Matlab code for the computation of f(A), where A is a Hermitian HSS matrix. 

The main algorithm is "hss2_funm_symm_telescop" which takes as an input a matrix in hss2 format. To convert a full matrix into hss2 one can use the routine "full_to_hss2".

This repository contains the code to reproduce all the numerical experiments from the preprint ``A. A. Casulli, D. Kressner, and L. Robol, Computing Functions of Symmetric Hierarchically Semiseparable Matrices, 2024, https://arxiv.org/abs/2402.17369''.

-- Correspondence between scripts and tables in the paper --

demo_Comparison_inv -> Table 1

demo_Comparison_inv_banded -> Table 2

demo_expm.m -> Table 3

demo_expm_laplace.m -> Table 4

testSampling2.m -> Table 5

demo_Comparison_inv_sqrt.m -> Table 6

demo_sign.m -> Table 7

-- Dependencies -- The following toolboxes are needed to run (some of) the scripts:

hm-toolbox https://github.com/numpi/hm-toolbox

rktoolbox http://guettel.com/rktoolbox/

chebfun https://www.chebfun.org/

