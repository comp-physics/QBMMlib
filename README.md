# QBMMlib

## Authors

* Spencer H. Bryngelson (Caltech) 
    * spencer@caltech.edu
* Tim Colonius (Caltech)
    * colonius@caltech.edu
* Rodney O. Fox (Iowa State)
    * rofox@iastate.edu

## Preprint

Check out the QBMMlib preprint on the [arXiv](https://arxiv.org/abs/2008.05063).

### Cite me!

```
@article{2008.05063v1,
Author        = {Spencer H. Bryngelson and Tim Colonius and Rodney O. Fox},
Title         = {{QBMMlib: A} library of quadrature-based moment methods},
Eprint        = {2008.05063},
ArchivePrefix = {arXiv},
Year          = {2020},
Month         = {Aug},
Url           = {http://arxiv.org/abs/2008.05063v1},
}
```

## Python version!

If Mathematica just isn't for you, we also developed a Python version of QBMMlib called [PyQBMMlib](https://github.com/sbryngelson/PyQBMMlib). Read more about it at the link.

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

## Installation

To install QBMMlib locally, open Wolfram Mathematica and choose the following options: 
* `File` -> `Install...`
    * `Type` -> `Package`
    * `Source` -> `QBMMlib.wl`
    * `InstallName` -> `QBMMlib`
    * `OK`

This places `QBMMlib.wl` where Mathematica can find it. You can then issue:
```
Get["QBMMlib"];
```
Alternatively, it can be loaded into any notebook without installation by copying `QBMMlib.wl` to the notebook directory and issuing
```
SetDirectory[NotebookDirectory[]]
Get[NotebookDirectory[] <> "QBMMlib.wl"];
```

## Acknowledgement
Great thanks is owed to Professor Alberto Passalacqua (Iowa State University) for his part in developing these algorithms and teaching me the same.
Funding was provided via the U.S. Office of Naval Research under grant numbers N0014-17-1-2676 and N0014-18-1-2625.
