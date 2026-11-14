# Quantum mechanics

This vignette covers `sympy.physics.quantum` (states, operators, Hilbert
spaces, spin, qubits, gates and quantum algorithms),
`sympy.physics.secondquant` (second quantization),
`sympy.physics.paulialgebra`, `sympy.physics.hep.gamma_matrices`, and
the solvable models and angular momentum coefficients of `sympy.physics`
(hydrogen atom, harmonic oscillators, particle on a ring, Wigner
symbols) together with `sympy.physics.matrices`.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y")
f <- Function("f"); g <- Function("g")
```

Several functions return R lists of expressions. The small helper
`tup()` shows such a list compactly as a (nested)
[`Tuple()`](https://adamrsophie.github.io/symcompR/reference/Tuple.md).

``` r

tup <- function(v) {
  do.call(Tuple, lapply(v, function(e) if (is.list(e) && !is_sym(e)) tup(e) else e))
}
```

## States, operators and Hilbert spaces

[`Ket()`](https://adamrsophie.github.io/symcompR/reference/quantum-states.md)
and
[`Bra()`](https://adamrsophie.github.io/symcompR/reference/quantum-states.md)
are abstract states,
[`Operator()`](https://adamrsophie.github.io/symcompR/reference/quantum-operators.md)
is an abstract operator, and
[`Dagger()`](https://adamrsophie.github.io/symcompR/reference/Dagger.md)
is the Hermitian conjugate. Products of a bra and a ket form inner and
outer products.

``` r

A <- Operator("A"); Bop <- Operator("B")
Dagger(A * Bop)
#> Dagger(B)*Dagger(A)
Dagger(Ket("psi"))
#> <psi|
Dagger(sym_I * A)
#> -I*Dagger(A)
Bra("phi") * Ket("psi")
#> <phi|psi>
Ket("psi") * Bra("phi")
#> |psi><phi|
OuterProduct(Ket("a"), Bra("b"))
#> |a><b|
Dagger(OuterProduct(Ket("a"), Bra("b")))
#> |b><a|
A * Bop - Bop * A
#> A*B - B*A
```

Quantum objects are Python classes in SymPy.
[`qattr()`](https://adamrsophie.github.io/symcompR/reference/quantum-attributes.md)
reads an attribute,
[`qcall()`](https://adamrsophie.github.io/symcompR/reference/quantum-attributes.md)
calls a method,
[`qclass()`](https://adamrsophie.github.io/symcompR/reference/quantum-attributes.md)
gives the class name, and
[`qisinstance()`](https://adamrsophie.github.io/symcompR/reference/quantum-attributes.md)
tests the class hierarchy (`QExpr`, `State`, `StateBase`, `KetBase`,
`BraBase` and so on).

``` r

k <- Ket("psi")
qattr(k, "label")
#> [[1]]
#> psi
qattr(k, "dual")
#> <psi|
qattr(k, "hilbert_space")
#> H
qcall(k, "dual_class")
#> [1] "Bra"
qclass(k)
#> [1] "Ket"
qisinstance(k, "QExpr")
#> [1] TRUE
qisinstance(k, "State")
#> [1] TRUE
qisinstance(k, "StateBase")
#> [1] TRUE
qisinstance(k, "KetBase")
#> [1] TRUE
qisinstance(Bra("psi"), "BraBase")
#> [1] TRUE
```

Time dependent states carry a time label, and orthogonal states have
inner products that evaluate to a Kronecker delta.

``` r

TimeDepKet("psi", "t")
#> |psi;t>
Dagger(TimeDepKet("psi", "t"))
#> <psi;t|
qattr(TimeDepBra("psi", "t"), "time")
#> t
qisinstance(TimeDepKet("psi", "t"), "TimeDepState")
#> [1] TRUE
OrthogonalBra(1) * OrthogonalKet(2)
#> <1|2>
doit(InnerProduct(OrthogonalBra(1), OrthogonalKet(2)))
#> 0
doit(InnerProduct(OrthogonalBra(1), OrthogonalKet(1)))
#> 1
qisinstance(OrthogonalKet(1), "OrthogonalState")
#> [1] TRUE
```

Special operators:
[`HermitianOperator()`](https://adamrsophie.github.io/symcompR/reference/quantum-operators.md),
[`UnitaryOperator()`](https://adamrsophie.github.io/symcompR/reference/quantum-operators.md),
[`IdentityOperator()`](https://adamrsophie.github.io/symcompR/reference/quantum-operators.md),
[`AntiCommutator()`](https://adamrsophie.github.io/symcompR/reference/AntiCommutator.md),
and the reduced Planck constant
[`HBar()`](https://adamrsophie.github.io/symcompR/reference/HBar.md)
(also available as `hbar`).

``` r

qcall(HermitianOperator("H"), "inverse")
#> H^(-1)
Dagger(HermitianOperator("H"))
#> H
Dagger(UnitaryOperator("U")) * UnitaryOperator("U")
#> 1
IdentityOperator() * A
#> A
represent(IdentityOperator(2))
#> Matrix(2, 2)
#> [1, 0]
#> [0, 1]
AntiCommutator(Bop, A)
#> {A,B}
doit(AntiCommutator(A, Bop))
#> A*B + B*A
HBar()
#> hbar
2 * hbar
#> 2*hbar
```

Hilbert spaces combine by direct sums (`+`), tensor products (`*`) and
tensor powers (`^`);
[`hs_contains()`](https://adamrsophie.github.io/symcompR/reference/hilbert-spaces.md)
tests whether a state lives in a space.

``` r

ComplexSpace(2) * FockSpace()
#> C(2)*F
ComplexSpace(2)^3
#> C(2)^3
qattr(ComplexSpace(2)^3, "dimension")
#> 8
qattr(ComplexSpace(2) + FockSpace(), "dimension")
#> oo
HilbertSpace()
#> H
L2(Interval(0, 1))
#> L2(Interval(0, 1))
DirectSumHilbertSpace(ComplexSpace(2), FockSpace())
#> C(2)+F
TensorProductHilbertSpace(ComplexSpace(2), ComplexSpace(3))
#> C(2)*C(3)
TensorPowerHilbertSpace(ComplexSpace(2), 3)
#> C(2)^3
hs_contains(ComplexSpace(2), Qubit("0"))
#> [1] TRUE
hs_contains(ComplexSpace(2), Ket("a"))
#> [1] FALSE
```

Tensor products of operators multiply factor by factor;
[`tensor_product_simp()`](https://adamrsophie.github.io/symcompR/reference/combined_tensor_printing.md)
performs that simplification and
[`combined_tensor_printing()`](https://adamrsophie.github.io/symcompR/reference/combined_tensor_printing.md)
switches the combined printing of tensor products of kets.

``` r

TensorProduct(A, Bop) * TensorProduct(Operator("C"), Operator("D"))
#> (A*C)x(B*D)
tensor_product_simp(TensorProduct(A, Bop) * TensorProduct(Operator("C"), Operator("D")))
#> (A*C)x(B*D)
combined_tensor_printing(FALSE)
TensorProduct(Ket("a"), Ket("b"))
#> |a>|b>
combined_tensor_printing(TRUE)
```

## Wavefunctions and representations

[`Wavefunction()`](https://adamrsophie.github.io/symcompR/reference/Wavefunction.md)
is a function of position together with its domain;
[`DifferentialOperator()`](https://adamrsophie.github.io/symcompR/reference/quantum-operators.md)
acts on wavefunctions.

``` r

L <- Symbol("L", positive = TRUE)
psi <- Wavefunction(sqrt(2 / L) * sin(sym_pi * x / L), list(x, 0, L))
qattr(psi, "norm")
#> 1
qattr(Wavefunction(x^2, x), "expr")
#> x^2
qcall(Wavefunction(x^2, list(x, 0, 1)), "normalize")
#> Wavefunction(sqrt(5)*x^2, (x, 0, 1))
qcall(Wavefunction(x^2, x), "__call__", 3)
#> 9
qcall(Wavefunction(x^2, list(x, 0, 1)), "__call__", 2)
#> 0
D <- DifferentialOperator(Derivative(f(x), x), f(x))
D
#> DifferentialOperator(Derivative(f(x), x),f(x))
qapply(D * Wavefunction(x^3, x))
#> Wavefunction(3*x^2, x)
```

[`qapply()`](https://adamrsophie.github.io/symcompR/reference/qapply.md)
applies operators to states, and
[`represent()`](https://adamrsophie.github.io/symcompR/reference/represent.md)
gives an operator or state in a basis. The position and momentum states
and operators of `sympy.physics.quantum.cartesian` are the usual
continuous bases.

``` r

qapply(OuterProduct(Ket("a"), Bra("b")) * Ket("c"))
#> <b|c>*|a>
qapply(XOp() * XKet(x))
#> x*|x>
XBra(x) * PxKet(y)
#> <x|y>
PxBra(y)
#> <y|
YOp()
#> Y
ZOp()
#> Z
PositionKet3D("x", "y", "z")
#> |xyz>
Dagger(PositionKet3D("x", "y", "z"))
#> <xyz|
PositionBra3D("x", "y", "z")
#> <xyz|
qisinstance(PositionKet3D("x", "y", "z"), "PositionState3D")
#> [1] TRUE
operators_to_state(XOp())
#> |x>
operators_to_state(PxOp())
#> |px>
state_to_operators(XKet())
#> X
represent(XOp(), basis = XKet())
#> x_1*DiracDelta(x_1 - x_2)
represent(XKet(x) * XBra(y), basis = XKet())
#> DiracDelta(x - x_2)*DiracDelta(x_1 - y)
```

The helpers behind
[`represent()`](https://adamrsophie.github.io/symcompR/reference/represent.md)
are exported too: the basis of an operator or state, enumerated states,
the inner product and expectation value forms, and the integration over
unities.

``` r

get_basis(XKet())
#> |x>
get_basis(PxOp())
#> |px>
tup(enumerate_states(XKet(x), 1, 2))
#> (|x_1>, |x_2>)
rep_innerproduct(XKet(x))
#> DiracDelta(x - x_1)
rep_expectation(XOp())
#> x_1*DiracDelta(x_1 - x_2)
x_1 <- Symbol("x_1"); x_2 <- Symbol("x_2")
integrate_result(XOp() * XKet(), x * DiracDelta(x - x_1) * DiracDelta(x_1 - x_2))
#> x*DiracDelta(x - x_1)*DiracDelta(x_1 - x_2)
integrate_result(XOp() * XKet(), x * DiracDelta(x - x_1) * DiracDelta(x_1 - x_2),
                 unities = list(1))
#> x*DiracDelta(x - x_2)
```

## Bosons, fermions and simple systems

[`BosonOp()`](https://adamrsophie.github.io/symcompR/reference/boson.md)
and
[`FermionOp()`](https://adamrsophie.github.io/symcompR/reference/fermion.md)
are annihilation operators; their daggers create particles in Fock
states.
[`normal_ordered_form()`](https://adamrsophie.github.io/symcompR/reference/operator-ordering.md)
moves creation operators to the left using the commutation relations,
while
[`normal_order()`](https://adamrsophie.github.io/symcompR/reference/operator-ordering.md)
only reorders.

``` r

a <- BosonOp("a")
Dagger(a)
#> Dagger(a)
qapply(a * BosonFockKet(2))
#> sqrt(2)*|1>
qapply(Dagger(a) * BosonFockKet(2))
#> sqrt(3)*|3>
doit(Commutator(a, Dagger(a)))
#> 1
BosonFockBra(1) * BosonFockKet(1)
#> <1|1>
qapply(a * BosonCoherentKet(x))
#> x*|x>
Dagger(BosonCoherentKet(x))
#> <x|
doit(BosonCoherentBra(x) * BosonCoherentKet(y))
#> exp(y*conjugate(x) - Abs(x)^2/2 - Abs(y)^2/2)
cf <- FermionOp("c")
doit(AntiCommutator(cf, Dagger(cf)))
#> 1
qapply(Dagger(cf) * FermionFockKet(0))
#> |1>
qapply(Dagger(cf) * FermionFockKet(1))
#> 0
FermionFockBra(0)
#> <0|
normal_ordered_form(a * Dagger(a))
#> 1 + Dagger(a)*a
normal_order(a * Dagger(a))
#> Dagger(a)*a
normal_ordered_form(cf * Dagger(cf))
#> 1 - Dagger(c)*c
```

The particle in a box and the one dimensional harmonic oscillator:

``` r

n <- Symbol("n"); m <- Symbol("m")
qapply(PIABHamiltonian("H") * PIABKet(n))
#> hbar^2*pi^2*n^2*|n>/(2*L^2*m)
PIABBra(n) * PIABKet(m)
#> <n|m>
qapply(RaisingOp("a") * SHOKet(1))
#> sqrt(2)*|2>
qapply(LoweringOp("a") * SHOKet(1))
#> |0>
qapply(NumberOp("N") * SHOKet(3))
#> 3*|3>
qapply(Hamiltonian("H") * SHOKet(Symbol("k")))
#> hbar*k*omega*|k> + hbar*omega*|k>/2
represent(NumberOp("N"), ndim = 3)
#> Matrix(3, 3)
#> [0, 0, 0]
#> [0, 1, 0]
#> [0, 0, 2]
represent(RaisingOp("a"), ndim = 3)
#> Matrix(3, 3)
#> [0, 0      , 0]
#> [1, 0      , 0]
#> [0, sqrt(2), 0]
doit(Commutator(LoweringOp("a"), RaisingOp("a")))
#> 1
SHOBra(1) * SHOKet(1)
#> <1|1>
doit(SHOBra(1) * SHOKet(1))
#> 1
```

## Pauli matrices

[`SigmaX()`](https://adamrsophie.github.io/symcompR/reference/pauli.md),
[`SigmaY()`](https://adamrsophie.github.io/symcompR/reference/pauli.md),
[`SigmaZ()`](https://adamrsophie.github.io/symcompR/reference/pauli.md),
[`SigmaPlus()`](https://adamrsophie.github.io/symcompR/reference/pauli.md)
and
[`SigmaMinus()`](https://adamrsophie.github.io/symcompR/reference/pauli.md)
are the Pauli operators, which act on
[`SigmaZKet()`](https://adamrsophie.github.io/symcompR/reference/pauli.md)
states;
[`qsimplify_pauli()`](https://adamrsophie.github.io/symcompR/reference/pauli.md)
simplifies their products.

``` r

sx <- SigmaX(1); sy <- SigmaY(1)
doit(Commutator(sx, sy))
#> 2*I*SigmaZ()
qsimplify_pauli(sx * sy)
#> I*SigmaZ()
qapply(sx * SigmaZKet(0))
#> |1>
qapply(SigmaZ() * SigmaZKet(1))
#> -|1>
qapply(SigmaMinus() * SigmaZKet(0))
#> |1>
qapply(SigmaPlus() * SigmaZKet(0))
#> 0
represent(SigmaY())
#> Matrix(2, 2)
#> [0, -I]
#> [I, 0 ]
represent(SigmaPlus())
#> Matrix(2, 2)
#> [0, 1]
#> [0, 0]
SigmaZBra(0) * SigmaZKet(0)
#> <0|0>
```

`sympy.physics.paulialgebra` has its own
[`Pauli()`](https://adamrsophie.github.io/symcompR/reference/pauli-algebra.md)
symbols with the product rule \sigma_i \sigma_j = \delta\_{ij} + i
\epsilon\_{ijk} \sigma_k.

``` r

Pauli(1) * Pauli(2)
#> I*sigma3
Pauli(3)^2
#> 1
evaluate_pauli_product(sym_I * Pauli(2) * Pauli(3))
#> -sigma1
```

## Spin and angular momentum

The spin operators `Jx`, `Jy`, `Jz`, `J2`, `Jplus` and `Jminus` act on
the eigenstates
[`JxKet()`](https://adamrsophie.github.io/symcompR/reference/spin-states.md),
[`JyKet()`](https://adamrsophie.github.io/symcompR/reference/spin-states.md)
and
[`JzKet()`](https://adamrsophie.github.io/symcompR/reference/spin-states.md).
[`JxOp()`](https://adamrsophie.github.io/symcompR/reference/spin-operators.md)
and its siblings create operators with other labels.

``` r

JzKet(1, 0)
#> |1,0>
qapply(Jz * JzKet(1, 1))
#> hbar*|1,1>
qapply(Jplus * JzKet(1, 0))
#> sqrt(2)*hbar*|1,1>
qapply(Jminus * JzKet(1, 0))
#> sqrt(2)*hbar*|1,-1>
qapply(J2 * JzKet(1, 0))
#> 2*hbar^2*|1,0>
doit(Commutator(Jx, Jy))
#> hbar*I*Jz
represent(Jx, j = S(1) / 2)
#> Matrix(2, 2)
#> [0     , hbar/2]
#> [hbar/2, 0     ]
represent(Jy, j = S(1) / 2)
#> Matrix(2, 2)
#> [0       , -hbar*I/2]
#> [hbar*I/2, 0        ]
represent(JzOp("J"), j = 1)
#> Matrix(3, 3)
#> [hbar, 0, 0    ]
#> [0   , 0, 0    ]
#> [0   , 0, -hbar]
JxOp("J"); JyOp("J"); J2Op("J"); JplusOp("J"); JminusOp("J")
#> Jx
#> Jy
#> J2
#> J+
#> J-
represent(JxKet(S(1) / 2, S(1) / 2), basis = Jz)
#> Matrix(2, 1)
#> [sqrt(2)/2]
#> [sqrt(2)/2]
rewrite(JxKet(1, 0), "Jz")
#> sqrt(2)*|1,-1>/2 - sqrt(2)*|1,1>/2
rewrite(JyKet(1, 1), "Jz")
#> -|1,-1>/2 + sqrt(2)*I*|1,0>/2 + |1,1>/2
JxBra(1, 0)
#> <1,0|
JyBra(1, 0)
#> <1,0|
JzBra(1, 0) * JzKet(1, 0)
#> <1,0|1,0>
doit(JzBra(1, 0) * JzKet(1, 0))
#> 1
qisinstance(JzKet(1, 0), "SpinState")
#> [1] TRUE
tup(m_values(S(3) / 2))
#> (4, (3/2, 1/2, -1/2, -3/2))
```

Coupled states:
[`couple()`](https://adamrsophie.github.io/symcompR/reference/spin-states.md)
and
[`uncouple()`](https://adamrsophie.github.io/symcompR/reference/spin-states.md)
convert between product states and coupled states.

``` r

half <- S(1) / 2
couple(TensorProduct(JzKet(half, half), JzKet(half, half)))
#> |1,1,j1=1/2,j2=1/2>
uncouple(JzKetCoupled(1, 0, list(half, half)))
#> sqrt(2)*|1/2,-1/2>|1/2,1/2>/2 + sqrt(2)*|1/2,1/2>|1/2,-1/2>/2
JxKetCoupled(1, 0, list(half, half))
#> |1,0,j1=1/2,j2=1/2>
JyKetCoupled(1, 0, list(half, half))
#> |1,0,j1=1/2,j2=1/2>
JzBraCoupled(1, 0, list(half, half))
#> <1,0,j1=1/2,j2=1/2|
JxBraCoupled(1, 0, list(half, half))
#> <1,0,j1=1/2,j2=1/2|
JyBraCoupled(1, 0, list(half, half))
#> <1,0,j1=1/2,j2=1/2|
qisinstance(JzKetCoupled(1, 0, list(half, half)), "CoupledSpinState")
#> [1] TRUE
```

Rotations and the Wigner D matrix:
[`Rotation_D()`](https://adamrsophie.github.io/symcompR/reference/spin-operators.md)
and
[`Rotation_d()`](https://adamrsophie.github.io/symcompR/reference/spin-operators.md)
are SymPy’s `Rotation.D` and `Rotation.d`.

``` r

doit(WignerD(1, 1, 0, 0, sym_pi / 2, 0))
#> -sqrt(2)/2
Rotation_D(1, 1, 0, sym_pi, sym_pi / 2, -sym_pi)
#> WignerD(1, 1, 0, pi, pi/2, -pi)
doit(Rotation_D(1, 1, 0, sym_pi, sym_pi / 2, -sym_pi))
#> sqrt(2)/2
doit(Rotation_d(1, 1, 0, sym_pi / 2))
#> -sqrt(2)/2
Rotation(sym_pi, 0, sym_pi / 2)
#> R(pi,0,pi/2)
qcall(Rotation(sym_pi, 0, sym_pi / 2), "inverse")
#> R(-pi/2,0,-pi)
```

Clebsch Gordan coefficients and Wigner symbols as unevaluated objects,
with
[`cg_simp()`](https://adamrsophie.github.io/symcompR/reference/cg-coefficients.md)
to simplify sums of them:

``` r

doit(CG(1, 0, 1, 0, 2, 0))
#> sqrt(6)/3
doit(CG(half, half, half, -half, 1, 0))
#> sqrt(2)/2
doit(Wigner3j(2, 0, 6, 0, 4, 0))
#> sqrt(715)/143
doit(Wigner6j(1, 2, 3, 2, 1, 2))
#> sqrt(21)/105
doit(Wigner9j(1, 1, 1, 1, 1, 1, 1, 1, 0))
#> 1/18
cg_simp(CG(1, 1, 0, 0, 1, 1) + CG(1, 0, 0, 0, 1, 0) + CG(1, -1, 0, 0, 1, -1))
#> 3
```

## Density matrices

[`Tr()`](https://adamrsophie.github.io/symcompR/reference/Tr.md) is the
trace, also of noncommutative symbols.
[`fidelity()`](https://adamrsophie.github.io/symcompR/reference/density-operators.md)
compares two density matrices.

``` r

Tr(Symbol("a") * Symbol("A", commutative = FALSE))
#> a*Tr(A)
Tr(Matrix(list(list(1, 1), list(2, 2))))
#> 3
doit(Density(list(JzKet(1, 0), 0.5), list(JzKet(1, 1), 0.5)))
#> 0.5*|1,0><1,0| + 0.5*|1,1><1,1|
up <- JzKet(half, half); down <- JzKet(half, -half)
updown <- (up + down) / sqrt(S(2))
fidelity(represent(up * Dagger(up)), represent(updown * Dagger(updown)))
#> sqrt(2)/2
fidelity(represent(up * Dagger(up)), represent(up * Dagger(up)))
#> 1
```

## Qubits and gates

[`Qubit()`](https://adamrsophie.github.io/symcompR/reference/quantum-qubits.md)
is a state of several qubits, written with the most significant qubit
first;
[`IntQubit()`](https://adamrsophie.github.io/symcompR/reference/quantum-qubits.md)
labels the same states by an integer.

``` r

Qubit("01")
#> |01>
QubitBra("01")
#> <01|
qattr(Qubit("01"), "dimension")
#> [1] 2
qubit_to_matrix(Qubit("01"))
#> Matrix(4, 1)
#> [0]
#> [1]
#> [0]
#> [0]
matrix_to_qubit(Matrix(list(0, 0, 1, 0)))
#> |10>
IntQubit(5)
#> |5>
qcall(IntQubit(5), "as_int")
#> [1] 5
IntQubit(2, nqubits = 3)
#> |2>
IntQubitBra(3)
#> <3|
Qubit(IntQubit(5))
#> |101>
qisinstance(Qubit("0"), "QubitState")
#> [1] TRUE
```

Gates act on qubits through
[`qapply()`](https://adamrsophie.github.io/symcompR/reference/qapply.md),
and
[`represent()`](https://adamrsophie.github.io/symcompR/reference/represent.md)
gives their matrices.

``` r

XGate(0)
#> X(0)
qapply(HadamardGate(0) * Qubit("1"))
#> sqrt(2)*|0>/2 - sqrt(2)*|1>/2
qapply(XGate(0) * Qubit("00"))
#> |01>
qapply(YGate(0) * Qubit("0"))
#> I*|1>
qapply(ZGate(0) * Qubit("1"))
#> -|1>
qapply(PhaseGate(0) * Qubit("1"))
#> I*|1>
qapply(TGate(0) * Qubit("1"))
#> exp(I*pi/4)*|1>
qapply(CNOT(1, 0) * Qubit("10"))
#> |11>
qapply(SWAP(0, 1) * Qubit("01"))
#> |10>
represent(CNOT(1, 0), nqubits = 2)
#> Matrix(4, 4)
#> [1, 0, 0, 0]
#> [0, 1, 0, 0]
#> [0, 0, 0, 1]
#> [0, 0, 1, 0]
represent(SWAP(0, 1), nqubits = 2)
#> Matrix(4, 4)
#> [1, 0, 0, 0]
#> [0, 0, 1, 0]
#> [0, 1, 0, 0]
#> [0, 0, 0, 1]
represent(CPHASE(1, 0), nqubits = 2)
#> Matrix(4, 4)
#> [1, 0, 0, 0 ]
#> [0, 1, 0, 0 ]
#> [0, 0, 1, 0 ]
#> [0, 0, 0, -1]
represent(TGate(0), nqubits = 1)
#> Matrix(2, 2)
#> [1, 0          ]
#> [0, exp(I*pi/4)]
represent(IdentityGate(0), nqubits = 1)
#> Matrix(2, 2)
#> [1, 0]
#> [0, 1]
CNotGate(1, 0)
#> CNOT(1,0)
SwapGate(0, 1)
#> SWAP(0,1)
qisinstance(XGate(0), "Gate")
#> [1] TRUE
qisinstance(XGate(0), "OneQubitGate")
#> [1] TRUE
qisinstance(SwapGate(0, 1), "TwoQubitGate")
#> [1] TRUE
```

`normalized(FALSE)` drops the factor 1/\sqrt{2} from the Hadamard gate;
`normalized(TRUE)` restores it.

``` r

normalized(FALSE)
represent(HadamardGate(0), nqubits = 1)
#> Matrix(2, 2)
#> [1, 1 ]
#> [1, -1]
normalized(TRUE)
```

Controlled gates and gates from a matrix:

``` r

CGate(list(0), YGate(1))
#> C((0),Y(1))
represent(CGate(list(0), YGate(1)), nqubits = 2)
#> Matrix(4, 4)
#> [1, 0, 0, 0 ]
#> [0, 0, 0, -I]
#> [0, 0, 1, 0 ]
#> [0, I, 0, 0 ]
CGateS(list(0), YGate(1))
#> C((0),Y(1))
UGate(list(0), Matrix(list(list(0, 1), list(1, 0))))
#> U((0,),Matrix([
#> [0, 1],
#> [1, 0]]))
represent(UGate(list(0), Matrix(list(list(0, 1), list(1, 0)))), nqubits = 1)
#> Matrix(2, 2)
#> [0, 1]
#> [1, 0]
represent_zbasis(list(1), list(0), Matrix(list(list(0, 1), list(1, 0))), 2)
#> Matrix(4, 4)
#> [1, 0, 0, 0]
#> [0, 1, 0, 0]
#> [0, 0, 0, 1]
#> [0, 0, 1, 0]
```

[`gate_simp()`](https://adamrsophie.github.io/symcompR/reference/quantum-gates.md)
and
[`gate_sort()`](https://adamrsophie.github.io/symcompR/reference/quantum-gates.md)
simplify and sort circuits:

``` r

gate_simp(HadamardGate(0) * HadamardGate(0))
#> 1
gate_simp(XGate(0) * XGate(0) * YGate(1))
#> Y(1)
gate_sort(YGate(1) * XGate(0) * ZGate(2))
#> X(0)*Y(1)*Z(2)
gate_sort(XGate(1) * HadamardGate(0)^2 * CNOT(0, 1) * XGate(1) * XGate(0))
#> CNOT(0,1)*X(0)
```

The matrix helpers of `sympy.physics.quantum.matrixutils`:

``` r

matrix_tensor_product(Matrix(list(list(1, 0), list(0, 1))),
                      Matrix(list(list(0, 1), list(1, 0))))
#> Matrix(4, 4)
#> [0, 1, 0, 0]
#> [1, 0, 0, 0]
#> [0, 0, 0, 1]
#> [0, 0, 1, 0]
matrix_dagger(Matrix(list(list(1, sym_I), list(0, 2))))
#> Matrix(2, 2)
#> [1 , 0]
#> [-I, 2]
matrix_to_zero(Matrix(list(list(0, 0), list(0, 0))))
#> 0
matrix_eye(2)
#> Matrix(2, 2)
#> [1, 0]
#> [0, 1]
matrix_zeros(2, 3)
#> Matrix(2, 3)
#> [0, 0, 0]
#> [0, 0, 0]
```

## Measurement

[`measure_all()`](https://adamrsophie.github.io/symcompR/reference/quantum-qubits.md)
and
[`measure_partial()`](https://adamrsophie.github.io/symcompR/reference/quantum-qubits.md)
list the possible outcomes with their probabilities;
[`matrix_to_density()`](https://adamrsophie.github.io/symcompR/reference/quantum-qubits.md)
turns a density matrix into a `Density` object. The `_oneshot` variants
simulate a single measurement, so their results are random.

``` r

tup(measure_all(qapply(HadamardGate(0) * Qubit("0"))))
#> ((|0>, 1/2), (|1>, 1/2))
bell <- qapply(HadamardGate(0) * HadamardGate(1) * Qubit("00"))
tup(measure_partial(bell, list(0)))
#> ((sqrt(2)*|00>/2 + sqrt(2)*|10>/2, 1/2), (sqrt(2)*|01>/2 + sqrt(2)*|11>/2, 1/2))
matrix_to_density(Matrix(list(list(1, 0), list(0, 0))))
#> (|0>, 1)
matrix_to_density(Matrix(list(list(half, 0), list(0, half))))
#> Density((|0>, 1/2),(|1>, 1/2))
set.seed(1)
measure_all_oneshot(bell)
#> |01>
measure_partial_oneshot(bell, list(0))
#> sqrt(2)*|00>/2 + sqrt(2)*|10>/2
```

## Quantum algorithms

The quantum Fourier transform and its inverse decompose into Hadamard,
controlled phase and swap gates.

``` r

QFT(0, 3)
#> QFT(0,3)
qcall(QFT(0, 3), "decompose")
#> SWAP(0,2)*H(0)*C((0),S(1))*H(1)*C((0),T(2))*C((1),S(2))*H(2)
qcall(IQFT(0, 2), "decompose")
#> H(1)*C((0),Rk(1,-2))*H(0)*SWAP(0,1)
represent(QFT(0, 2), nqubits = 2)
#> Matrix(4, 4)
#> [1/2, 1/2 , 1/2 , 1/2 ]
#> [1/2, I/2 , -1/2, -I/2]
#> [1/2, -1/2, 1/2 , -1/2]
#> [1/2, -I/2, -1/2, I/2 ]
RkGate(1, 2)
#> S(1)
Rk(1, 2)
#> S(1)
represent(RkGate(0, 3), nqubits = 1)
#> Matrix(2, 2)
#> [1, 0          ]
#> [0, exp(I*pi/4)]
```

Grover’s search:
[`OracleGate()`](https://adamrsophie.github.io/symcompR/reference/grover.md)
flips the sign of the marked states,
[`WGate()`](https://adamrsophie.github.io/symcompR/reference/grover.md)
is the diffusion operator, and
[`apply_grover()`](https://adamrsophie.github.io/symcompR/reference/grover.md)
runs the whole algorithm.

``` r

oracle_fun <- function(q) isTRUE(q == IntQubit(2, nqubits = 2))
superposition_basis(2)
#> |0>/2 + |1>/2 + |2>/2 + |3>/2
qapply(OracleGate(2, oracle_fun) * IntQubit(2, nqubits = 2))
#> -|2>
qapply(grover_iteration(superposition_basis(2), OracleGate(2, oracle_fun)))
#> |2>
apply_grover(oracle_fun, 2)
#> |2>
WGate(2)
#> W(2)
qapply(WGate(1) * Qubit("0"))
#> |0> + |1> - |0>
```

Shor’s algorithm uses the continued fraction helpers
[`getr()`](https://adamrsophie.github.io/symcompR/reference/shor.md) and
[`ratioize()`](https://adamrsophie.github.io/symcompR/reference/shor.md).
The controlled modular multiplication gate
[`CMod()`](https://adamrsophie.github.io/symcompR/reference/shor.md) is
not finished in SymPy, so
[`CMod()`](https://adamrsophie.github.io/symcompR/reference/shor.md),
[`period_find()`](https://adamrsophie.github.io/symcompR/reference/shor.md)
and [`shor()`](https://adamrsophie.github.io/symcompR/reference/shor.md)
stop with SymPy’s `NotImplementedError`.

``` r

getr(513, 1024, 10)
#> 2
ratioize(list(1, 2, 3), 10)
#> 6
tryCatch(CMod(4, 2, 15), error = conditionMessage)
#> [1] "The CMod gate has not been completed."
tryCatch(period_find(2, 15), error = conditionMessage)
#> [1] "The CMod gate has not been completed."
set.seed(7)
tryCatch(shor(15), error = conditionMessage)
#> [1] "The CMod gate has not been completed."
```

## Circuits, QASM and identities

[`Qasm()`](https://adamrsophie.github.io/symcompR/reference/qasm.md)
reads the QASM circuit description language, from separate lines, from a
string
([`read_qasm()`](https://adamrsophie.github.io/symcompR/reference/qasm.md)),
or from a file
([`read_qasm_file()`](https://adamrsophie.github.io/symcompR/reference/qasm.md)).
The gate constructors of the circuit plotting module create named gates
and measurement gates.

``` r

Qasm("qubit q0", "qubit q1", "h q0", "cnot q0,q1")$get_circuit()
#> CNOT(1,0)*H(1)
read_qasm("qubit q0\nqubit q1\nh q0\ncnot q0,q1")$get_circuit()
#> CNOT(1,0)*H(1)
qasm_file <- tempfile(fileext = ".qasm")
writeLines(c("qubit q0", "qubit q1", "h q0", "cnot q0,q1"), qasm_file)
read_qasm_file(qasm_file)$get_circuit()
#> CNOT(1,0)*H(1)
labeller(3)
#> [1] "q_2" "q_1" "q_0"
Mz(0)
#> Mz(0)
Mx(0)
#> Mx(0)
CreateOneQubitGate("Q")(1)
#> Q(1)
CreateCGate("V")(list(0), 1)
#> C((0),V(1))
```

Circuit utilities find and replace subcircuits, and convert between
numbered and symbolic qubit indices.

``` r

tup(kmp_table(list("A", "B", "C", "D", "A", "B", "D")))
#> (-1, 0, 0, 0, 0, 1, 2)
find_subcircuit(XGate(0) * YGate(0) * XGate(0), YGate(0) * XGate(0))
#> [1] 1
tup(replace_subcircuit(XGate(0) * YGate(0) * XGate(0), list(YGate(0))))
#> (X(0), X(0))
tup(replace_subcircuit(XGate(0) * YGate(0) * XGate(0), list(YGate(0)),
                       list(ZGate(0), ZGate(0))))
#> (X(0), Z(0), Z(0), X(0))
tup(convert_to_symbolic_indices(list(XGate(0), YGate(1), CNOT(1, 0)))[[1]])
#> (X(i0), Y(i1), CNOT(i1,i0))
symbolic <- convert_to_symbolic_indices(list(XGate(5), YGate(7)))
tup(convert_to_real_indices(symbolic[[1]], symbolic[[2]]))
#> (X(5), Y(7))
```

A gate identity is a circuit equal to a scalar multiple of the identity.
[`GateIdentity()`](https://adamrsophie.github.io/symcompR/reference/identitysearch.md)
stores one together with its equivalent rearrangements, and the identity
search functions look for identities among a set of gates. The rules and
equivalent identities are sets, so only their sizes are shown here.

``` r

gi <- GateIdentity(XGate(0), XGate(0))
gi
#> 1
length(generate_gate_rules(list(XGate(0), YGate(0), ZGate(0))))
#> [1] 12
length(generate_equivalent_ids(list(XGate(0), YGate(0), ZGate(0))))
#> [1] 6
is_scalar_matrix(list(XGate(0), XGate(0)), 1, FALSE)
#> [1] TRUE
is_scalar_sparse_matrix(list(XGate(0), XGate(0)), 1, FALSE)
#> [1] TRUE
is_scalar_nonsparse_matrix(list(XGate(0), YGate(0)), 1, FALSE)
#> [1] FALSE
is_reducible(list(XGate(0), YGate(0), ZGate(0)), 1, 0, 3)
#> [1] TRUE
is_degenerate(list(GateIdentity(XGate(0), YGate(0), ZGate(0))),
              GateIdentity(YGate(0), ZGate(0), XGate(0)))
#> [1] FALSE
length(bfs_identity_search(list(XGate(0), YGate(0), ZGate(0)), 1, max_depth = 2))
#> [1] 3
```

The random circuit tools draw gates at random:

``` r

set.seed(7)
random_circuit(4, 2)
#> Y(0)*CNOT(1,0)*Y(0)*Y(1)
tup(random_insert(list(XGate(0), YGate(0)), list(list(XGate(0), XGate(0))), seed = 1))
#> (X(0), X(0), X(0), Y(0))
tup(random_reduce(XGate(0) * YGate(0) * YGate(0),
                  list(GateIdentity(YGate(0), YGate(0))), seed = 1))
#> (X(0),)
tup(random_identity_search(list(XGate(0), YGate(0), ZGate(0)), 2, 1))
#> ()
```

## Second quantization

`sympy.physics.secondquant` has bosonic
([`B()`](https://adamrsophie.github.io/symcompR/reference/secondquant.md),
[`Bd()`](https://adamrsophie.github.io/symcompR/reference/secondquant.md))
and fermionic
([`sym_F()`](https://adamrsophie.github.io/symcompR/reference/secondquant.md),
[`Fd()`](https://adamrsophie.github.io/symcompR/reference/secondquant.md))
annihilation and creation operators acting on occupation number states.
Its `Dagger`, `Commutator` and `InnerProduct` are exported as
[`sq_Dagger()`](https://adamrsophie.github.io/symcompR/reference/secondquant.md),
[`sq_Commutator()`](https://adamrsophie.github.io/symcompR/reference/secondquant.md)
and
[`sq_InnerProduct()`](https://adamrsophie.github.io/symcompR/reference/secondquant.md).

``` r

p <- Symbol("p"); q <- Symbol("q")
B(0)
#> b(0)
Bd(0)
#> b+(0)
sq_Dagger(B(0))
#> b+(0)
BKet(list(1, 2))
#> |(1, 2)>
BBra(list(1, 2))
#> <(1, 2)|
apply_operators(Bd(0) * BKet(list(1, 2)))
#> sqrt(2)*FockStateBosonKet((2, 2))
apply_operators(B(1) * BKet(list(1, 2)))
#> sqrt(2)*FockStateBosonKet((1, 1))
sq_InnerProduct(BBra(list(1, 2)), BKet(list(1, 2)))
#> 1
sq_Commutator(B(0), Bd(0))
#> 1
sq_Commutator(sym_F(p), Fd(q))
#> KroneckerDelta(_a, q)*KroneckerDelta(p, q) - (KroneckerDelta(_i, p)*KroneckerDelta(p, q) + NO(CreateFermion(q)*AnnihilateFermion(p))) - NO(CreateFermion(q)*AnnihilateFermion(p))
Fd(p) * sym_F(q)
#> CreateFermion(p)*AnnihilateFermion(q)
FKet(list(1, 2))
#> |(1, 2)>
FBra(list(1, 2))
#> <(1, 2)|
apply_operators(Fd(3) * FKet(list(1, 2)))
#> |(1, 2, 3)>
apply_operators(sym_F(1) * FKet(list(1, 2)))
#> |(2,)>
```

The underlying classes can also be called directly:

``` r

AnnihilateBoson(0)
#> b(0)
CreateBoson(0)
#> b+(0)
AnnihilateFermion(p)
#> f(p)
CreateFermion(p)
#> f+(p)
FockStateBosonKet(list(1, 2))
#> |(1, 2)>
FockStateBosonBra(list(1, 2))
#> <(1, 2)|
FockStateFermionKet(list(1, 2))
#> |(1, 2)>
FockStateFermionBra(list(1, 2))
#> <(1, 2)|
SqOperator(0)
#> sq(0)
BosonicOperator(0)
#> sq(0)
FermionicOperator(p)
#> sq(p)
FockState(list(1, 2))
#> (1, 2)
FockStateKet(list(1, 2))
#> |(1, 2)>
FockStateBra(list(1, 2))
#> <(1, 2)|
BosonState(list(1, 2))
#> (1, 2)
FermionState(list(1, 2))
#> (1, 2)
```

Bosonic bases and matrix representations:

``` r

VarBosonicBasis(3)
#> [FockState((0,)), FockState((1,)), FockState((2,))]
FixedBosonicBasis(2, 2)
#> [FockState((2, 0)), FockState((1, 1)), FockState((0, 2))]
matrix_rep(Bd(0) * B(0), VarBosonicBasis(3))
#> Matrix(3, 3)
#> [0, 0, 0]
#> [0, 1, 0]
#> [0, 0, 2]
matrix_rep(B(0), VarBosonicBasis(3))
#> Matrix(3, 3)
#> [0, 1, 0      ]
#> [0, 0, sqrt(2)]
#> [0, 0, 0      ]
```

Indices above and below the Fermi level are symbols with the
`above_fermi` and `below_fermi` assumptions.
[`NO()`](https://adamrsophie.github.io/symcompR/reference/secondquant-wicks.md)
is the normal ordered product,
[`contraction()`](https://adamrsophie.github.io/symcompR/reference/secondquant-wicks.md)
the contraction of two operators, and
[`wicks()`](https://adamrsophie.github.io/symcompR/reference/secondquant-wicks.md)
applies Wick’s theorem.

``` r

i <- Symbol("i", below_fermi = TRUE); j <- Symbol("j", below_fermi = TRUE)
ab <- Symbol("a", above_fermi = TRUE); bb <- Symbol("b", above_fermi = TRUE)
FKet(list(i, ab), fermi_level = 1)
#> -FockStateFermionKet((a, i))
NO(sym_F(i) * Fd(ab))
#> -NO(CreateFermion(a)*AnnihilateFermion(i))
NO(Fd(ab) * sym_F(i))
#> :CreateFermion(a)*AnnihilateFermion(i):
contraction(Fd(i), sym_F(j))
#> KroneckerDelta(i, j)
contraction(sym_F(ab), Fd(bb))
#> KroneckerDelta(a, b)
contraction(sym_F(p), Fd(q))
#> KroneckerDelta(_a, q)*KroneckerDelta(p, q)
wicks(sym_F(p) * Fd(q))
#> KroneckerDelta(_a, q)*KroneckerDelta(p, q) - NO(CreateFermion(q)*AnnihilateFermion(p))
wicks(Fd(p) * sym_F(q))
#> KroneckerDelta(_i, q)*KroneckerDelta(p, q) + NO(CreateFermion(p)*AnnihilateFermion(q))
wicks(sym_F(i) * Fd(j) * sym_F(ab) * Fd(bb))
#> KroneckerDelta(a, b)*NO(AnnihilateFermion(i)*CreateFermion(j)) + NO(CreateFermion(b)*AnnihilateFermion(i)*AnnihilateFermion(a)*CreateFermion(j))
evaluate_deltas(KroneckerDelta(i, j) * sym_F(j))
#> f(i)
evaluate_deltas(KroneckerDelta(p, q) * sym_F(q))
#> f(p)
```

Antisymmetric tensors, index permutations and dummy substitution:

``` r

AntiSymmetricTensor("v", list(ab, i), list(bb, j))
#> v((a, i),(b, j))
AntiSymmetricTensor("v", list(i, ab), list(bb, j))
#> -AntiSymmetricTensor(v, (a, i), (b, j))
PermutationOperator(p, q)
#> PermutationOperator(p, q)
qcall(PermutationOperator(p, q), "get_permuted", f(p, q))
#> -f(q, p)
simplify_index_permutations(f(p) * g(q) - f(q) * g(p), list(PermutationOperator(p, q)))
#> f(p)*g(q)*PermutationOperator(p, q)
substitute_dummies(
  AntiSymmetricTensor("t", list(Dummy("a", above_fermi = TRUE)), list(Dummy("i", below_fermi = TRUE))) -
  AntiSymmetricTensor("t", list(Dummy("b", above_fermi = TRUE)), list(Dummy("j", below_fermi = TRUE))))
#> 0
```

## Dirac gamma matrices

[`GammaMatrix()`](https://adamrsophie.github.io/symcompR/reference/gamma-matrices.md)
is the Dirac gamma matrix with a Lorentz index (`LorentzIndex` is the
index type).
[`kahane_simplify()`](https://adamrsophie.github.io/symcompR/reference/gamma-matrices.md)
removes contracted gamma matrices,
[`gamma_trace()`](https://adamrsophie.github.io/symcompR/reference/gamma-matrices.md)
takes traces, and
[`simplify_gpgp()`](https://adamrsophie.github.io/symcompR/reference/gamma-matrices.md)
simplifies products of slashed momenta.

``` r

i0 <- TensorIndex("i0", LorentzIndex); i1 <- TensorIndex("i1", LorentzIndex)
i2 <- TensorIndex("i2", LorentzIndex); i3 <- TensorIndex("i3", LorentzIndex)
kahane_simplify(GammaMatrix(i0) * GammaMatrix(i1) * GammaMatrix(-i0))
#> -2*GammaMatrix(i1)
kahane_simplify(GammaMatrix(i0) * GammaMatrix(-i0))
#> Matrix(4, 4)
#> [4, 0, 0, 0]
#> [0, 4, 0, 0]
#> [0, 0, 4, 0]
#> [0, 0, 0, 4]
kahane_simplify(GammaMatrix(i0) * GammaMatrix(i1) * GammaMatrix(i2) * GammaMatrix(-i0))
#> 2*(GammaMatrix(i1)*GammaMatrix(i2) + GammaMatrix(i2)*GammaMatrix(i1))
gamma_trace(GammaMatrix(i0) * GammaMatrix(i1))
#> 4*metric(i0, i1)
gamma_trace(GammaMatrix(i0) * GammaMatrix(i1) * GammaMatrix(i2) * GammaMatrix(i3))
#> -4*metric(i0, i2)*metric(i1, i3) + 4*metric(i0, i1)*metric(i2, i3) + 4*metric(i0, i3)*metric(i1, i2)
gamma_trace(GammaMatrix(i0) * GammaMatrix(-i0))
#> 16
simplify_gamma_expression(GammaMatrix(i0) * GammaMatrix(-i0))
#> Matrix(4, 4)
#> [4, 0, 0, 0]
#> [0, 4, 0, 0]
#> [0, 0, 4, 0]
#> [0, 0, 0, 4]
pm <- TensorHead("p", list(LorentzIndex)); qm <- TensorHead("q", list(LorentzIndex))
simplify_gpgp(pm(i0) * GammaMatrix(-i0) * pm(i1) * GammaMatrix(-i1))
#> 4*p(L_0)*p(-L_0)
extract_type_tens(pm(i0) * GammaMatrix(-i0) * qm(i1), GammaMatrix)
#> [[1]]
#> GammaMatrix(-L_0) 
#> 
#> [[2]]
#> p(L_0)*q(i1)
```

## Solvable models

The hydrogen atom (`sympy.physics.hydrogen`), the three dimensional
isotropic harmonic oscillator (`sympy.physics.sho`), the one dimensional
harmonic oscillator (`sympy.physics.qho_1d`) and the particle on a ring
(`sympy.physics.pring`). The functions carry the module name as a
prefix.

``` r

r <- Symbol("r"); theta <- Symbol("theta"); phi <- Symbol("phi")
hydrogen_R_nl(1, 0, r)
#> 2*exp(-r)
hydrogen_R_nl(2, 1, r, 1)
#> sqrt(6)*r*exp(-r/2)/12
hydrogen_E_nl(2)
#> -1/8
hydrogen_E_nl_dirac(2, 0)
#> -0.125002080189006
hydrogen_E_nl_dirac(2, 0, Z = 1, c = 1)
#> -1 + sqrt(2)/2
hydrogen_Psi_nlm(2, 1, 1, r, phi, theta)
#> -r*exp(-r/2)*exp(I*phi)*sin(theta)/(8*sqrt(pi))
hydrogen_Y_lm(1, 0, theta, phi)
#> sqrt(3)*cos(phi)/(2*sqrt(pi))
hydrogen_Z_lm(1, 1, theta, phi)
#> -sqrt(3)*exp(I*theta)*sin(phi)/(4*sqrt(pi)) - sqrt(3)*exp(-I*theta)*sin(phi)/(4*sqrt(pi))
sho_R_nl(0, 0, 1, r)
#> 2*2^(3/4)*exp(-r^2)/pi^(1/4)
sho_E_nl(1, 2, Symbol("hw"))
#> 11*hw/2
qho_psi_n(0, x, 1, 1)
#> exp(-x^2/(2*hbar))/(hbar^(1/4)*pi^(1/4))
qho_psi_n(1, x, Symbol("m"), Symbol("omega"))
#> sqrt(2)*x*(m*omega)^(3/4)*exp(-m*omega*x^2/(2*hbar))/(hbar^(3/4)*pi^(1/4))
qho_E_n(1, Symbol("omega"))
#> 3*hbar*omega/2
qho_coherent_state(2, 0.5)
#> 0.110312112823074*sqrt(2)
pring_wavefunction(1, x)
#> sqrt(2)*exp(I*x)/(2*sqrt(pi))
pring_energy(1, Symbol("m"), Symbol("r"))
#> hbar^2/(2*m*r^2)
```

## Wigner symbols and Gaunt coefficients

`sympy.physics.wigner` evaluates the angular momentum coupling
coefficients exactly.

``` r

wigner_3j(2, 6, 4, 0, 0, 0)
#> sqrt(715)/143
clebsch_gordan(Rational(3, 2), S_Half, 2, Rational(3, 2), S_Half, 2)
#> 1
wigner_6j(3, 3, 3, 3, 3, 3)
#> -1/14
wigner_9j(1, 1, 1, 1, 1, 1, 1, 1, 0)
#> 1/18
wigner_9j(1, 1, 1, 1, 1, 1, 1, 1, 0, prec = 20)
#> 0.055555555555555555556
racah(3, 3, 3, 3, 3, 3)
#> -1/14
gaunt(1, 0, 1, 1, 0, -1)
#> -1/(2*sqrt(pi))
real_gaunt(2, 2, 2, 0, 2, 2)
#> -sqrt(5)/(7*sqrt(pi))
doit(dot_rot_grad_Ynm(3, 2, 2, 0, theta, phi))
#> 3*sqrt(55)*Ynm(5, 2, theta, phi)/(11*sqrt(pi))
wigner_d_small(S_Half, S_Pi / 2)
#> Matrix(2, 2)
#> [sqrt(2)/2 , sqrt(2)/2]
#> [-sqrt(2)/2, sqrt(2)/2]
wigner_d(S_Half, Symbol("a"), Symbol("b"), Symbol("c"))
#> Matrix(2, 2)
#> [exp(I*a/2)*exp(I*c/2)*cos(b/2)  , exp(I*a/2)*exp(-I*c/2)*sin(b/2) ]
#> [-exp(-I*a/2)*exp(I*c/2)*sin(b/2), exp(-I*a/2)*exp(-I*c/2)*cos(b/2)]
```

## Physics matrices

`sympy.physics.matrices` has the Pauli matrices, the Dirac gamma
matrices, the Minkowski metric, the parallel axis theorem matrix, and
the (deprecated) discrete Fourier transform matrix
[`mdft()`](https://adamrsophie.github.io/symcompR/reference/msigma.md),
which now points to
[`DFT()`](https://adamrsophie.github.io/symcompR/reference/structured-matrices.md).

``` r

msigma(1)
#> Matrix(2, 2)
#> [0, 1]
#> [1, 0]
msigma(1) %*% msigma(2)
#> Matrix(2, 2)
#> [I, 0 ]
#> [0, -I]
mgamma(0)
#> Matrix(4, 4)
#> [1, 0, 0 , 0 ]
#> [0, 1, 0 , 0 ]
#> [0, 0, -1, 0 ]
#> [0, 0, 0 , -1]
mgamma(5)
#> Matrix(4, 4)
#> [0, 0, 1, 0]
#> [0, 0, 0, 1]
#> [1, 0, 0, 0]
#> [0, 1, 0, 0]
mgamma(1, TRUE)
#> Matrix(4, 4)
#> [0, 0, 0 , -1]
#> [0, 0, -1, 0 ]
#> [0, 1, 0 , 0 ]
#> [1, 0, 0 , 0 ]
minkowski_tensor
#> Matrix(4, 4)
#> [1, 0 , 0 , 0 ]
#> [0, -1, 0 , 0 ]
#> [0, 0 , -1, 0 ]
#> [0, 0 , 0 , -1]
pat_matrix(2, 1, 0, 0)
#> Matrix(3, 3)
#> [0, 0, 0]
#> [0, 2, 0]
#> [0, 0, 2]
suppressWarnings(mdft(2))
#> Matrix(2, 2)
#> [sqrt(2)/2, sqrt(2)/2 ]
#> [sqrt(2)/2, -sqrt(2)/2]
```
