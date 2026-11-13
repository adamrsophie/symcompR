# Expected values are SymPy outputs (sympy.physics.quantum and paulialgebra).

test_that("states, operators and their products print as SymPy", {
  psi <- Ket("psi"); phi <- Bra("phi")
  A <- Symbol("A", commutative = FALSE); B_ <- Symbol("B", commutative = FALSE)
  expect_equal(sstr(psi), "|psi>")
  expect_equal(sstr(phi), "<phi|")
  expect_equal(sstr(phi * psi), "<phi|psi>")
  expect_equal(sstr(Dagger(A)), "Dagger(A)")
  expect_equal(sstr(Commutator(A, B_)), "[A,B]")
  expect_equal(sstr(AntiCommutator(A, B_)), "{A,B}")
  expect_equal(srepr(psi), "Ket(Symbol('psi'))")
  expect_equal(latex(phi * psi), "\\left\\langle \\phi \\right. {\\left|\\psi\\right\\rangle }")
})

test_that("spin operators act on spin states", {
  expect_equal(sstr(qapply(Jz * JzKet(1, 1))), "hbar*|1,1>")
  expect_equal(sstr(represent(Jz, basis = Jz)), "Matrix([[hbar/2, 0], [0, -hbar/2]])")
  h <- S(1) / 2
  expect_equal(sstr(couple(TensorProduct(JzKet(h, h), JzKet(h, h)))), "|1,1,j1=1/2,j2=1/2>")
  expect_equal(sstr(doit(CG(1, 0, 1, 0, 2, 0))), "sqrt(6)/3")
})

test_that("gates act on qubits", {
  expect_equal(sstr(qapply(HadamardGate(0) * Qubit("0"))), "sqrt(2)*|0>/2 + sqrt(2)*|1>/2")
  expect_equal(sstr(qapply(CNOT(1, 0) * Qubit("10"))), "|11>")
})

test_that("bosons, Pauli matrices and traces", {
  a <- BosonOp("a")
  expect_equal(sstr(doit(Commutator(a, Dagger(a)))), "1")
  expect_equal(sstr(evaluate_pauli_product(Pauli(1) * Pauli(2))), "I*sigma3")
  expect_equal(sstr(Pauli(3)^3), "sigma3")
  expect_equal(sstr(qsimplify_pauli(SigmaX() * SigmaY())), "I*SigmaZ()")
  A <- Symbol("A", commutative = FALSE); B_ <- Symbol("B", commutative = FALSE)
  expect_true(identical(Tr(A * B_)$key, Tr(B_ * A)$key))
})

test_that("parse_latex reads bras, kets and inner products", {
  expect_equal(sstr(parse_latex("\\langle x | y \\rangle", backend = "lark")), "<x|y>")
  expect_equal(sstr(parse_latex("| x \\rangle", backend = "lark")), "|x>")
  expect_equal(sstr(parse_latex("\\langle x |", backend = "lark")), "<x|")
})

test_that("fidelity of states with Float weights", {
  h <- S(1) / 2
  up <- JzKet(h, h); down <- JzKet(h, -h)
  updown <- up / sqrt(S(2)) + down / sqrt(S(2))
  f <- fidelity(Density(list(updown, Float(1.0))), Density(list(down, Float(1.0))))
  expect_lt(abs(Re(as.complex(N(f))) - 1 / sqrt(2)), 1e-3)
})

test_that("Density and entropy are shared with sympy.stats", {
  d <- Density(list(Ket(0), 0.5), list(Ket(1), 0.5))
  expect_equal(d$type, "QDensity")
  expect_equal(qclass(d), "Density")
  expect_equal(sstr(doit(d)), "0.5*|0><0| + 0.5*|1><1|")
  h <- S(1) / 2
  expect_equal(sstr(entropy(Density(list(JzKet(h, h), h), list(JzKet(h, -h), h)))), "log(2)/2")
  expect_equal(Density(Die("X", 6))$type, "Density")
})

test_that("an integer seed reproduces SymPy's random reduction", {
  x <- XGate(0); y <- YGate(0); z <- ZGate(0); h <- HadamardGate(0)
  ids <- bfs_identity_search(list(x, y, z), 1)
  out <- random_reduce(list(x, y, z, x, y, h), ids, seed = 3)
  expect_equal(vapply(out, sstr, ""), c("X(0)", "Y(0)", "H(0)"))
  expect_equal(Density(Tuple(Ket(0), S(1) / 2))$type, "QDensity")
})
