# Cases from sympy/tensor/tests/test_tensor.py (test_tensor_matching,
# test_TensAdd_matching, test_TensMul_matching, test_TensMul_subs); the
# expected strings are SymPy's output.

.wt_setup <- function() {
  R3 <- TensorIndexType("R3", dim = 3)
  ii <- tensor_indices("p q r s t", R3)
  env <- list2env(ii)
  for (n in c("a", "b", "c", "d", "e", "f")) assign(n, WildTensorIndex(n, R3, ignore_updown = TRUE), envir = env)
  env$R3 <- R3
  env$g <- WildTensorIndex("g", R3)
  env$delta <- R3$delta; env$eps <- R3$epsilon
  env$K <- TensorHead("K", list(R3)); env$V <- TensorHead("V", list(R3)); env$A <- TensorHead("A", list(R3, R3))
  env$W <- WildTensorHead("W", unordered_indices = TRUE); env$U <- WildTensorHead("U")
  env$wi <- Wild("wi"); env$k <- Symbol("K")
  env
}
.wt_fmt <- function(m) {
  if (is.null(m)) return(NULL)
  paste(names(m), vapply(m, function(v) if (inherits(v, "sym")) sstr(v) else format(v), ""), sep = ": ", collapse = ", ")
}
.wt_s <- function(x) if (inherits(x, "sym")) sstr(x) else format(x)

test_that("wild indices and wild tensors match as in SymPy", {
  with(.wt_setup(), {
    expect_equal(.wt_fmt(sym_match(q, a)), "a: q")
    expect_equal(.wt_fmt(sym_match(-q, a)), "a: -q")
    expect_null(sym_match(-q, g))
    expect_equal(.wt_fmt(sym_match(q, g)), "g: q")
    expect_null(sym_match(eps(p, q, r), eps(p, -a, a)))
    expect_equal(.wt_fmt(sym_match(eps(p, q, r), eps(p, -b, a))), "-b: q, a: r")
    expect_equal(.wt_s(sym_replace(eps(p, -q, r), eps(a, b, c), 1)), "1")
    expect_equal(.wt_fmt(sym_match(K(p) * V(q), W())), "W: K(p)*V(q)")
    expect_equal(.wt_fmt(sym_match(K(p), W(a))), "a: p, W(R3): _WildTensExpr(K(p))")
    expect_equal(.wt_fmt(sym_match(K(p) * V(q), W(a, p))), "a: q, W(R3,R3): _WildTensExpr(K(p)*V(q))")
    expect_equal(.wt_fmt(sym_match(K(p) * V(q), W(p, q))), "W(R3,R3): _WildTensExpr(K(p)*V(q))")
    expect_equal(.wt_fmt(sym_match(A(q, p), W(p, q))), "W(R3,R3): _WildTensExpr(A(q, p))")
    expect_null(sym_match(A(q, p), U(p, q)))
    expect_equal(.wt_s(sym_replace(K(q) * K(p), W(q, p), 1)), "1")
    expect_equal(length(sym_match(delta(q, p), delta(p, q))), 0L)
    expect_null(sym_match(eps(q, p, r), eps(p, q, r)))
    expect_false(is.null(sym_match(eps(q, r, p), eps(p, q, r))))
    expect_equal(format(W(a)), "W(a)")
    expect_equal(format(W()), "W")
    expect_equal(format(-a), "-a")
    expect_s3_class(-a, "WildTensorIndex")
  })
})

test_that("sums of tensors match as in SymPy", {
  with(.wt_setup(), {
    expect_equal(length(sym_match(K(p) * K(q) + V(p) * V(q), K(p) * K(q) + V(p) * V(q))), 0L)
    expect_false(is.null(sym_match(K(p) * K(q) + V(p) * V(q), K(p) * K(q) + V(p) * V(q))))
    expect_null(sym_match(K(p) * K(q) + V(p) * V(q) + K(p) * V(q) + V(p) * K(q), K(p) * K(q) + V(p) * V(q)))
    r1 <- sym_replace(K(p) * K(q) + V(p) * V(q) + K(p) * V(q) + K(q) * V(p),
                      W(p, q) + K(p) * K(q) + V(p) * V(q), W(p, q) + 3 * K(p) * V(q))
    expect_equal(format(.tn_doit(r1)), "K(q)*V(p) + 4*K(p)*V(q)")
  })
})

test_that("products of tensors match and replace as in SymPy", {
  with(.wt_setup(), {
    expect_equal(.wt_fmt(sym_match(K(p), wi * K(p))), "wi: 1")
    expect_equal(.wt_fmt(sym_match(eps(p, r, q), wi * eps(p, q, r))), "wi: -1")
    expect_equal(.wt_s(sym_replace(K(p) * V(-p), W(a) * V(-a), 1)), "1")
    expect_equal(.wt_s(sym_replace(K(q) * K(p) * V(-p), W(q, a) * V(-a), 1)), "1")
    expect_equal(.wt_s(sym_replace(K(p) * V(-p), K(-a) * V(a), 1)), "1")
    expect_equal(.wt_s(sym_replace(K(q) * K(p) * V(-p), W(q) * U(p) * V(-p), 1)), "1")
    expect_equal(format(.tn_doit(sym_replace(K(p) * V(q), W() * K(p) * V(q), W() * V(p) * V(q)))), "V(p)*V(q)")
    r1 <- sym_replace(eps(r, p, q) * eps(-r, -s, -t), eps(e, a, b) * eps(-e, c, d),
                      delta(a, c) * delta(b, d) - delta(a, d) * delta(b, c))
    expect_equal(format(canon_bp(.tn_doit(r1))), "KD(p, -s)*KD(q, -t) - KD(p, -t)*KD(q, -s)")
    r2 <- sym_replace(eps(r, p, q) * eps(-r, -p, -q), eps(c, a, b) * eps(-c, d, f),
                      delta(a, d) * delta(b, f) - delta(a, f) * delta(b, d))
    expect_equal(.wt_s(.tn_doit(contract_delta(r2, delta))), "6")
    expect_equal(format(.tn_doit(sym_replace(V(-p) * V(q) * V(-q), wi * W() * V(a) * V(-a), wi * W()))), "V(-p)")
    expect_equal(.wt_s(.tn_doit(sym_replace(k^4 * K(r) * K(-r), wi * W() * K(a) * K(-a), wi * W() * k^2))), "K^6")
    expect_equal(format(sym_replace(K(p) * V(q), W(q) * K(p), W(p) * W(q))), "V(p)*V(q)")
    expect_equal(format(sym_replace(K(p) * V(q) * V(r), W(q, r) * K(p), W(p, r) * W(q, s) * V(-s))),
                 "V(p)*V(r)*V(q)*V(R_0)*V(-R_0)")
    D <- tensor_indices("R_0 R_1 R_2 R_3", R3)
    m <- sym_match(delta(-D[[1]], -D[[2]]) * K(D[[3]]) * K(D[[4]]) * K(-D[[4]]), W() * K(a) * K(-a))
    expect_equal(.wt_fmt(m), "R_0: R_3, -R_0: -R_3, W: KD(-R_0, -R_1)*K(R_2)")
    expect_equal(format(sym_replace(K(p) * V(-p), K(a), V(a) * V(q) * V(-q))), "V(R_0)*V(R_1)*V(-R_1)*V(-R_0)")
  })
})
