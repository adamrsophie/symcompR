# Random matrix ensembles

The Gaussian and circular ensembles of `sympy.stats`.
`joint_eigen_distribution()` and `level_spacing_distribution()` give the
joint distribution of the eigenvalues and the level spacing density;
`JointEigenDistribution()` is the joint distribution of the eigenvalues
of a matrix with random entries.

## Usage

``` r
GaussianEnsemble(sym, dim)

GaussianUnitaryEnsemble(sym, dim)

GaussianOrthogonalEnsemble(sym, dim)

GaussianSymplecticEnsemble(sym, dim)

CircularEnsemble(sym, dim)

CircularUnitaryEnsemble(sym, dim)

CircularOrthogonalEnsemble(sym, dim)

CircularSymplecticEnsemble(sym, dim)

joint_eigen_distribution(mat)

level_spacing_distribution(mat)

JointEigenDistribution(mat)

RandomMatrixSymbol(sym, n, m, pspace = NULL)
```

## Arguments

- sym:

  The name of the random matrix.

- dim:

  Its dimension.

- mat:

  A random matrix symbol, or a matrix of random variables.

- n, m:

  Dimensions of a random matrix symbol.

- pspace:

  A probability space, or `NULL`.

## Value

A random matrix symbol, a `Lambda`, or a distribution.

## Examples

``` r
G <- GaussianUnitaryEnsemble("U", 2)
level_spacing_distribution(G)
#> Lambda(_s, 32*_s^2*exp(-4*_s^2/pi)/pi^2)
```
