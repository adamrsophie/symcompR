# Expected values are SymPy outputs (sympy.physics.secondquant).

test_that("second quantization follows Wick's theorem", {
  p <- Symbol("p"); q <- Symbol("q")
  i <- Symbol("i", below_fermi = TRUE); a <- Symbol("a", above_fermi = TRUE)
  expect_equal(sstr(wicks(sym_F(p) * Fd(q))),
               "KroneckerDelta(_a, q)*KroneckerDelta(p, q) - NO(CreateFermion(q)*AnnihilateFermion(p))")
  expect_equal(sstr(NO(sym_F(i) * Fd(a))), "-NO(CreateFermion(a)*AnnihilateFermion(i))")
  expect_equal(sstr(contraction(Fd(i), sym_F(i))), "1")
  expect_equal(sstr(apply_operators(Bd(0) * BKet(list(1, 2)))), "sqrt(2)*FockStateBosonKet((2, 2))")
  expect_equal(sstr(sq_Commutator(B(0), Bd(0))), "1")
  expect_equal(sstr(AntiSymmetricTensor("v", list(a, i), list(i, a))), "-AntiSymmetricTensor(v, (a, i), (a, i))")
  expect_equal(latex(Fd(p)), "{a^\\dagger_{p}}")
})

test_that("nested commutators, Fock kets in LaTeX, bases and doit hints", {
  i <- Symbol("i", below_fermi = TRUE); j <- Symbol("j", below_fermi = TRUE)
  a <- Symbol("a", above_fermi = TRUE)
  X <- Symbol("X", commutative = FALSE); Y <- Symbol("Y", commutative = FALSE)
  Z <- Symbol("Z", commutative = FALSE)
  expect_equal(sstr(sq_Commutator(Y, sq_Commutator(X, Z))), "[Y,[X,Z]]")
  expect_equal(latex(FockStateFermionKet(list(a, i), 1)), "\\left|\\left( a, \\  i\\right)\\right\\rangle")
  expect_equal(sstr(doit(sq_Commutator(Fd(j), NO(Fd(a) * sym_F(i))), wicks = TRUE)),
               "-KroneckerDelta(i, j)*CreateFermion(a)")
  expect_equal(sstr(FixedBosonicBasis(2, 2)), "[FockState((2, 0)), FockState((1, 1)), FockState((0, 2))]")
})
