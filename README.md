# IsolatedPoints 
### (A Maple library for computing real isolated points of an algebraic hypersurface) 
### (in development)

## About

IsolatedPoints.mla is a Maple library for computing real isolated points of an algebraic hypersurface.

Given a multivariate polynomial `f`, the call of function `isolatedPoints(f)` computes a data `[param,boxes]` where:
* `param` is a zero-dimensional parametrization that encodes a finite set `C` of points containing all the real isolated points of the set `f = 0`,
* `boxes` is a list of boxes that identify the desired real isolated points among the elements of `C`.

For example, the call

```
with(IsolatedPoints):
isolatedPoints(x^2-y^3+y^2);
```
returns
```
[param = [0,0,u], boxes = [[[0,0],[0,0]]]]
```
that identifies the isolated point `(0,0)`.

## Installation

It requires the following software to function correctly
* [Maple 2018+](https://www.maplesoft.com/),
* [FGb](https://www-polsys.lip6.fr/~jcf/FGb/) (v. 1.70) for Grobner basis computations,
* [RAGlib](https://www-polsys.lip6.fr/~safey/RAGLib/) (v. 3.54).

The library is provided as a file `IsolatedPoints.mla` which should be loaded in Maple by:

1. Add the following line in your mapleinit file (usually found at $HOME/.mapleinit)

```
libname := "PATH/TO/IsolatedPoints.mla", libname:
```

2. In each Maple session, you can load the library by

```
with(IsolatedPoints):
```
