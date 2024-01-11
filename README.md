## Introduction
This paper considers a class of multiblock nonconvex nonsmooth optimization problems. 
To this end, we propose an accelerated block proximal linear method with adaptive momentum ($ABPL^+$) for solving such problems. Compared to the previous methods, the advantages of our method are: (1) it extends applicability to multivariable optimization problems and allows to randomly picking the block of variables to update,  (2) it allows the extrapolation parameter to be independent of other parameters, (3) it has an adaptive extrapolation parameter strategy.  While maintaining the above advantages, we prove that our method can ensure the convergence of this class of problems, and we also prove that the sequence generated by our method is globally convergent to a first-order critical point as well as establish the convergence rate of our method.  Sparse multilayer nonnegative matrix factorization with $\ell_0$-constraints and sparse nonnegative CP decomposition with $\ell_0$-constraints are known to be NP-hard, nonconvex and nonsmooth problems.  To address these challenges, we apply our method to solve these problems. The numerical experimental results on solving these problems demonstrate that our method outperforms the state-of-the-art algorithms. 
This package contains code for the sparse multiple non-negative matrix factorization with $\ell_0$-constraints ($\ell_0$-smNMF) problem and the sparse non-negative CP tensor decomposition with $\ell_0$-constraints ($\ell_0$-SNCP) problem. 

## Matlab code
A toy example explains how to use the $ABPL^+$ function. For "L0smNMF", just run the function 'main_Run_me' . For "L0SNCP", before running the script, please first add the toolbox 'tensortoolbox'[<sup>2</sup>](#refer-id) to the running path of matlab, and then run the function 'main_Run_me'. More descriptions about these functions can be found in their annotation part.

## Data
This code has built-in the data mentioned in paper. 

## Reference
<div id="refer-id"></div>
<!-- [1] Yang W, Min W. An Accelerated Block Proximal Framework with Adaptive Momentum for Nonconvex and Nonsmooth Optimization -->

[1] Brett W. Bader and Tamara G. Kolda. 2006. Algorithm 862: MATLAB tensor classes for fast algorithm prototyping. ACM Trans. Math. Softw. 32, 4 (December 2006), 635–653. https://doi.org/10.1145/1186785.1186794
