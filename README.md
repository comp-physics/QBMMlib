# QBMMlib

## Authors

Spencer H. Bryngelson, Rodney O. Fox, Tim Colonius 

## Abstract

QBMMlib is an open source Mathematica package of quadrature-based moment methods and their algorithms.
Such methods are commonly used to solve fully-coupled disperse flow and combustion problems, though formulating and closing the corresponding governing equations can be complex.
QBMMlib aims to make analyzing these techniques simple and more accessible.
Its routines use symbolic manipulation to formulate the moment transport equations for a population balance equation and a prescribed dynamical system.
However, the resulting moment transport equations are unclosed.
QBMMlib trades the moments for a set of quadrature points and weights via an inversion algorithm, of which several are available.
Quadratures then closes the moment transport equations.

## Files

* `QBMMlib.wl` Package file
* `Examples.nb` Example notebook 

## Acknowledgement
Great thanks is to Professor Alberto Passalacqua (Iowa State University) for his part in developing these algorithms and teaching me the same.
Funding was provided via the U.S. Office of Naval Research under grant numbers N0014-17-1-2676 and N0014-18-1-2625.

## References

* Patel, R. G., Desjardins, O., & Fox, R. O. (2019). Three-dimensional conditional hyperbolic quadrature method of moments. Journal of Computational Physics: X, 1, 100006. https://doi.org/10.1016/j.jcpx.2019.100006
* Yuan, C., & Fox, R. O. (2011). Conditional quadrature method of moments for kinetic equations. Journal of Computational Physics, 230(22), 8216–8246. https://doi.org/10.1016/j.jcp.2011.07.020
* Marchisio, D. L., Fox, R. O., & Fox, R. O. (2010). Computational models for polydisperse particulate and multiphase systems. Cambridge University Press (Vol. 9780521858). https://doi.org/10.1017/CBO9781139016599
* Bryngelson, S. H., Charalampopoulos, A., Sapsis, T. P., & Colonius, T. (2020). A Gaussian moment method and its augmentation via LSTM recurrent neural networks for the statistics of cavitating bubble populations. International Journal of Multiphase Flow, 127, 103262. https://doi.org/10.1016/j.ijmultiphaseflow.2020.103262
