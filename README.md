## Introduction
This paper considers a class of multiblock nonsmooth nonconvex optimization problems. To this end, we propose an accelerated block proximal linear method with adaptive momentum for solving such problems. Compared to the previous methods for solving this class of problems, the advantages of our method are (1) it extends applicability to multivariable optimization problems and allows to randomly picking the block of variables to update,  (2) it allows the extrapolation parameter to be independent of other parameters, (3) it has an adaptive extrapolation parameter strategy.  While maintaining the above advantages, we prove that our method can ensure the convergence of this class of problems, and also prove that it is globally convergent to a first-order critical point as well as establish its convergence rates.  
The numerical experiments, including additional experiments, demonstrate that our  method outperforms state-of-the-art algorithms.  This package contains code for the sparse multiple non-negative matrix factorization with 0 norm ($\ell_0$-smNMF) problem and the sparse non-negative CP tensor decomposition with 0 norm ($\ell_0$-SNCP) problem described in the paper[<sup>1</sup>](#refer-id). 

## Matlab code
A toy example explains how to use the $ABPL^+$ function. For "L0smNMF", just run the function 'main_Run_me' . For "L0SNCP", before running the script, please first add the toolbox 'tensortoolbox'[<sup>2</sup>](#refer-id) to the running path of matlab, and then run the function 'main_Run_me'. More descriptions about these functions can be found in their annotation part.

## Data
This code has built-in the data mentioned in paper[<sup>1</sup>](#refer-id). 

## Reference
<div id="refer-id"></div>
[1] Yang W, Min W. An Accelerated Block Proximal Framework with Adaptive Momentum for Nonconvex and Nonsmooth Optimization[J]. arXiv preprint arXiv:2308.12126, 2023.

[2] Brett W. Bader and Tamara G. Kolda. 2006. Algorithm 862: MATLAB tensor classes for fast algorithm prototyping. ACM Trans. Math. Softw. 32, 4 (December 2006), 635–653. https://doi.org/10.1145/1186785.1186794
