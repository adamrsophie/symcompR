# Further _eval_rewrite_as_<target> methods of SymPy: the erfc, sign,
# hyper, SingularityFunction, asin and tractable targets. rewrite() merges
# these with its own table. Each rule takes the arguments of the node and
# returns NULL to leave the node alone.

.fs_hyp <- function(ap, bq, z) hyper(ap, bq, z)

.fs_erf_hyper <- function(z) .mul(list(S_Two, z, .inv(.sqrt(S_Pi)), .fs_hyp(list(S_Half), list(.rat(3, 2)), .neg(.sq(z)))))

# Heaviside._eval_rewrite_as_sign
.fs_heaviside_sign <- function(arg, H0 = S_Half) {
  if (!isTRUE(ask_fact(arg, "extended_real"))) return(NULL)
  sg <- .div(.add(list(fcall("sign", list(arg)), S_One)), S_Two)
  h0 <- fcall("Heaviside", list(S_Zero, H0))
  pw1 <- Piecewise(list(sg, Ne(arg, S_Zero)), list(h0, TRUE))
  Piecewise(list(sg, Eq(h0, S_Half)), list(pw1, TRUE))
}

# DiracDelta and Heaviside _eval_rewrite_as_SingularityFunction
.fs_singularity <- function(name, args) {
  free <- free_symbols(fnode(name, args))
  if (name == "DiracDelta" && !length(free)) {
    if (length(args) == 1L && is_zero_num(args[[1L]])) return(fcall("SingularityFunction", list(S_Zero, S_Zero, S_NegativeOne)))
    if (length(args) == 2L && is_zero_num(args[[1L]]) && is_one_num(args[[2L]]))
      return(fcall("SingularityFunction", list(S_Zero, S_Zero, .int(-2L))))
  }
  if (length(free) != 1L)
    stop(.sympy_condition("TypeError", "\nrewrite(SingularityFunction) does not support arguments with more that\none variable."))
  x <- free[[1L]]
  a <- .fs_solve_first(args[[1L]], x)
  n <- if (name == "Heaviside") S_Zero else if (length(args) == 1L) S_NegativeOne else .sub(.neg(args[[2L]]), S_One)
  fcall("SingularityFunction", list(x, a, n))
}

# solve(f, x)[0]
.fs_solve_first <- function(f, x) {
  s <- solve(f, x)
  if (is_sym(s)) s else s[[1L]]
}

.fs_tractable_rule <- function(nm) {
  force(nm)
  function(...) .gr_tractable1(fnode(nm, list(...)), NULL)
}

# built once per session
.fs_memo <- new.env(parent = emptyenv())
.fs_rewrite_rules <- function() {
  if (is.null(.fs_memo$rewrite)) .fs_memo$rewrite <- .fs_make_rewrite_rules()
  .fs_memo$rewrite
}

.fs_make_rewrite_rules <- function() {
  hp <- .mul(list(S_Half, S_Pi))
  tr <- c("gamma", "uppergamma", "erf", "erfc", "erfi", "Ei", "li", "Li", "fibonacci", "harmonic",
          "RisingFactorial", "FallingFactorial", "binomial", "besseli", "besselk", "sinh", "cosh",
          "tanh", "coth", "csch", "sech", .GR_LOGREW)
  list(
    erfc = list(
      erf = function(z) .sub(S_One, fcall("erfc", list(z))),
      erfi = function(z) .sub(.mul(list(S_I, fcall("erfc", list(.mul(list(S_I, z)))))), S_I),
      erf2 = function(x, y) .sub(fcall("erfc", list(x)), fcall("erfc", list(y)))
    ),
    sign = list(
      Heaviside = .fs_heaviside_sign,
      Abs = function(a) .div(a, fcall("sign", list(a)))
    ),
    hyper = list(
      erf = .fs_erf_hyper,
      erfc = function(z) .sub(S_One, .fs_erf_hyper(z)),
      erfi = function(z) .mul(list(S_Two, z, .inv(.sqrt(S_Pi)), .fs_hyp(list(S_Half), list(.rat(3, 2)), .sq(z)))),
      erf2 = function(x, y) .sub(.fs_erf_hyper(y), .fs_erf_hyper(x)),
      elliptic_k = function(m) .mul(list(S_Pi, S_Half, .fs_hyp(list(S_Half, S_Half), list(S_One), m))),
      elliptic_e = function(...) {
        a <- list(...)
        if (length(a) != 1L) return(NULL)
        .mul(list(hp, .fs_hyp(list(.rat(-1, 2), S_Half), list(S_One), a[[1L]])))
      },
      li = function(z) {
        lz <- fcall("log", list(z))
        .add(list(.mul(list(lz, .fs_hyp(list(S_One, S_One), list(S_Two, S_Two), lz))),
                  .mul(list(S_Half, .sub(fcall("log", list(lz)), fcall("log", list(.inv(lz)))))),
                  S_EulerGamma))
      },
      fresnels = function(z) .mul(list(S_Pi, .pow(z, .int(3L)), .rat(1, 6),
                                       .fs_hyp(list(.rat(3, 4)), list(.rat(3, 2), .rat(7, 4)),
                                               .neg(.div(.mul(list(.sq(S_Pi), .pow(z, .int(4L)))), .int(16L)))))),
      fresnelc = function(z) .mul(list(z, .fs_hyp(list(.rat(1, 4)), list(S_Half, .rat(5, 4)),
                                                  .neg(.div(.mul(list(.sq(S_Pi), .pow(z, .int(4L)))), .int(16L)))))),
      airyai = function(z) {
        pf1 <- .inv(.mul(list(.pow(.int(3L), .rat(2, 3)), fcall("gamma", list(.rat(2, 3))))))
        pf2 <- .div(z, .mul(list(.pow(.int(3L), .rat(1, 3)), fcall("gamma", list(.rat(1, 3))))))
        .sub(.mul(list(pf1, .fs_hyp(list(), list(.rat(2, 3)), .div(.pow(z, .int(3L)), .int(9L))))),
             .mul(list(pf2, .fs_hyp(list(), list(.rat(4, 3)), .div(.pow(z, .int(3L)), .int(9L))))))
      },
      airybi = function(z) {
        pf1 <- .inv(.mul(list(.pow(.int(3L), .rat(1, 6)), fcall("gamma", list(.rat(2, 3))))))
        pf2 <- .div(.mul(list(z, .pow(.int(3L), .rat(1, 6)))), fcall("gamma", list(.rat(1, 3))))
        .add(list(.mul(list(pf1, .fs_hyp(list(), list(.rat(2, 3)), .div(.pow(z, .int(3L)), .int(9L))))),
                  .mul(list(pf2, .fs_hyp(list(), list(.rat(4, 3)), .div(.pow(z, .int(3L)), .int(9L)))))))
      },
      airyaiprime = function(z) {
        pf1 <- .div(.sq(z), .mul(list(S_Two, .pow(.int(3L), .rat(2, 3)), fcall("gamma", list(.rat(2, 3))))))
        pf2 <- .inv(.mul(list(.pow(.int(3L), .rat(1, 3)), fcall("gamma", list(.rat(1, 3))))))
        .sub(.mul(list(pf1, .fs_hyp(list(), list(.rat(5, 3)), .div(.pow(z, .int(3L)), .int(9L))))),
             .mul(list(pf2, .fs_hyp(list(), list(.rat(1, 3)), .div(.pow(z, .int(3L)), .int(9L))))))
      },
      airybiprime = function(z) {
        pf1 <- .div(.sq(z), .mul(list(S_Two, .pow(.int(3L), .rat(1, 6)), fcall("gamma", list(.rat(2, 3))))))
        pf2 <- .div(.pow(.int(3L), .rat(1, 6)), fcall("gamma", list(.rat(1, 3))))
        .add(list(.mul(list(pf1, .fs_hyp(list(), list(.rat(5, 3)), .div(.pow(z, .int(3L)), .int(9L))))),
                  .mul(list(pf2, .fs_hyp(list(), list(.rat(1, 3)), .div(.pow(z, .int(3L)), .int(9L)))))))
      },
      catalan = function(n) .fs_hyp(list(.sub(S_One, n), .neg(n)), list(S_Two), S_One),
      betainc = function(a, b, x1, x2) {
        h <- function(x) .mul(list(.pow(x, a), .fs_hyp(list(a, .sub(S_One, b)), list(.add(list(a, S_One))), x)))
        .div(.sub(h(x2), h(x1)), a)
      },
      betainc_regularized = function(a, b, x1, x2) {
        h <- function(x) .mul(list(.pow(x, a), .fs_hyp(list(a, .sub(S_One, b)), list(.add(list(a, S_One))), x)))
        .div(.div(.sub(h(x2), h(x1)), a), fcall("beta", list(a, b)))
      }
    ),
    SingularityFunction = list(
      DiracDelta = function(...) .fs_singularity("DiracDelta", list(...)),
      Heaviside = function(...) .fs_singularity("Heaviside", list(...)[1L])
    ),
    asin = list(
      acos = function(x) .sub(hp, fcall("asin", list(x))),
      atan = function(a) .mul(list(.sqrt(.sq(a)), .inv(a), .sub(hp, fcall("asin", list(.inv(.sqrt(.add(list(S_One, .sq(a)))))))))),
      acot = function(a) .mul(list(a, .sqrt(.inv(.sq(a))),
                                   .sub(hp, fcall("asin", list(.div(.sqrt(.neg(.sq(a))), .sqrt(.sub(.neg(.sq(a)), S_One)))))))),
      asec = function(a) .sub(hp, fcall("asin", list(.inv(a)))),
      acsc = function(a) fcall("asin", list(.inv(a))),
      asinh = function(x) .mul(list(.neg(S_I), fnode("asin", list(.mul(list(S_I, x)))))),
      acosh = function(x) .mul(list(.sqrt(.sub(x, S_One)), .inv(.sqrt(.sub(S_One, x))), .sub(hp, fcall("asin", list(x)))))
    ),
    tractable = stats::setNames(lapply(tr, .fs_tractable_rule), tr)
  )
}
