# Second quantization

Creation and annihilation operators of bosons and fermions, occupation
number (Fock) states and the algebra of `sympy.physics.secondquant`.
`B(k)` and `Bd(k)` annihilate and create a boson in state `k`;
`sym_F(k)` and `Fd(k)` do the same for a fermion. A symbolic index may
carry the assumption `above_fermi = TRUE` (a particle state) or
`below_fermi = TRUE` (a hole state). `BKet()`, `BBra()`, `FKet()` and
`FBra()` are occupation number states. `sq_Dagger()` is the Hermitian
conjugate, `sq_Commutator()` the commutator, `sq_InnerProduct()` the
inner product of a bra and a ket and `apply_operators()` applies the
operators of an expression to the states on their right.

## Usage

``` r
sq_Dagger(arg)

AntiSymmetricTensor(symbol, upper, lower)

SqOperator(k)

BosonicOperator(k)

AnnihilateBoson(k)

CreateBoson(k)

B(k)

Bd(k)

FermionicOperator(k)

AnnihilateFermion(k)

CreateFermion(k)

sym_F(k)

Fd(k)

FockState(occupations)

BosonState(occupations)

FermionState(occupations, fermi_level = 0)

FockStateKet(occupations)

FockStateBra(occupations)

FockStateBosonKet(occupations)

FockStateBosonBra(occupations)

FockStateFermionKet(occupations, fermi_level = 0)

FockStateFermionBra(occupations, fermi_level = 0)

BKet(occupations)

BBra(occupations)

FKet(occupations, fermi_level = 0)

FBra(occupations, fermi_level = 0)

apply_operators(e)

sq_InnerProduct(bra, ket)

matrix_rep(op, basis)

VarBosonicBasis(n_max)

FixedBosonicBasis(n_particles, n_levels)

sq_Commutator(a, b)
```

## Arguments

- arg, a, b, e:

  Expressions.

- symbol, upper, lower:

  The name, upper and lower indices of an antisymmetric tensor.

- k:

  The state (index) acted upon: an integer or a symbol.

- occupations:

  The occupation numbers (bosons) or the occupied states (fermions), as
  a list.

- fermi_level:

  The number of states below the Fermi level.

- bra, ket:

  A bra and a ket.

- op:

  An operator.

- basis:

  A basis from `VarBosonicBasis()` or `FixedBosonicBasis()`.

- n_max:

  The number of states of a variable particle number basis.

- n_particles, n_levels:

  The number of particles and of levels of a fixed particle number
  basis.

## Value

An expression, a list or, for the bases, an object of class
`VarBosonicBasis` or `FixedBosonicBasis`.

## Examples

``` r
p <- Symbol("p"); q <- Symbol("q")
Fd(p) * sym_F(q)
#> CreateFermion(p)*AnnihilateFermion(q)
sq_Dagger(B(0))
#> b+(0)
apply_operators(Bd(0) * BKet(list(1, 2)))
#> sqrt(2)*FockStateBosonKet((2, 2))
sq_Commutator(B(0), Bd(0))
#> 1
```
