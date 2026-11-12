# Rewrites between the Bessel functions, after the _eval_rewrite_as_*
# methods of sympy/functions/special/bessel.py.

local({
  f <- function(name, ...) fcall(name, list(...))
  sq <- function(z) .sqrt(z)
  hlf <- S_Half
  notint <- function(nu) isFALSE(ask_fact(nu, "integer"))
  jtoi <- function(nu, z) .mul(list(fcall("exp", list(.mul(list(S_I, S_Pi, nu, hlf)))), f("besseli", nu, .mul(list(fcall("polar_lift", list(.neg(S_I))), z)))))
  itoj <- function(nu, z) .mul(list(fcall("exp", list(.mul(list(.neg(S_I), S_Pi, nu, hlf)))), f("besselj", nu, .mul(list(fcall("polar_lift", list(S_I)), z)))))
  ytoj <- function(nu, z) if (notint(nu)) .mul(list(fcall("csc", list(.mul(list(S_Pi, nu)))), .sub(.mul(list(fcall("cos", list(.mul(list(S_Pi, nu)))), f("besselj", nu, z))), f("besselj", .neg(nu), z)))) else NULL
  ktoi <- function(nu, z) if (notint(nu)) .mul(list(S_Pi, fcall("csc", list(.mul(list(S_Pi, nu)))), .sub(f("besseli", .neg(nu), z), f("besseli", nu, z)), hlf)) else NULL
  # assume_integer_order: the rewrite applies to integer orders only
  intord <- function(g) function(nu, z) if (isTRUE(ask_fact(nu, "integer"))) g(nu, z) else NULL
  then <- function(g, target) function(nu, z) { r <- g(nu, z); if (is.null(r)) NULL else rewrite(r, target) }
  sqp2z <- function(z) sq(.div(S_Pi, .mul(list(S_Two, z))))
  rules <- list(
    besseli = list(
      besselj = jtoi,
      bessely = then(ytoj, "besseli"),
      besselk = ktoi),
    besselj = list(
      besseli = itoj,
      bessely = ytoj,
      besselk = then(ktoi, "besselj"),
      jn = function(nu, z) .mul(list(sqp2z(z), f("besselj", .add(list(nu, hlf)), z))),
      yn = intord(function(nu, z) .mul(list(.pow(S_NegativeOne, .add(list(nu, S_One))), sqp2z(z), f("besselj", .sub(.neg(nu), hlf), z)))),
      hn1 = intord(function(nu, z) .mul(list(sqp2z(z), .add(list(f("besselj", .add(list(nu, hlf)), z),
                                     .mul(list(S_I, .pow(S_NegativeOne, .add(list(nu, S_One))), f("besselj", .sub(.neg(nu), hlf), z)))))))),
      hn2 = intord(function(nu, z) .mul(list(sqp2z(z), .add(list(f("besselj", .add(list(nu, hlf)), z),
                                     .mul(list(.neg(S_I), .pow(S_NegativeOne, .add(list(nu, S_One))), f("besselj", .sub(.neg(nu), hlf), z)))))))),
      hankel1 = NULL),
    bessely = list(
      besselj = function(nu, z) if (notint(nu)) .sub(.mul(list(fcall("csc", list(.mul(list(S_Pi, nu)))), f("bessely", .neg(nu), z))),
                                                       .mul(list(fcall("cot", list(.mul(list(S_Pi, nu)))), f("bessely", nu, z)))) else NULL,
      besseli = then(itoj, "bessely"),
      besselk = then(then(ktoi, "besselj"), "bessely"),
      jn = function(nu, z) .mul(list(.pow(S_NegativeOne, nu), sqp2z(z), f("bessely", .sub(.neg(nu), hlf), z))),
      yn = intord(function(nu, z) .mul(list(sqp2z(z), f("bessely", .add(list(nu, hlf)), z)))),
      hn1 = intord(function(nu, z) .mul(list(sqp2z(z), .add(list(.mul(list(.pow(S_NegativeOne, nu), f("bessely", .sub(.neg(nu), hlf), z))),
                                     .mul(list(S_I, f("bessely", .add(list(nu, hlf)), z)))))))),
      hn2 = intord(function(nu, z) .mul(list(sqp2z(z), .add(list(.mul(list(.pow(S_NegativeOne, nu), f("bessely", .sub(.neg(nu), hlf), z))),
                                     .mul(list(.neg(S_I), f("bessely", .add(list(nu, hlf)), z))))))))),
    jn = list(
      besselj = function(nu, z) .mul(list(sq(.div(.mul(list(S_Two, z)), S_Pi)), f("jn", .sub(nu, hlf), z))),
      besseli = function(nu, z) rewrite(itoj(nu, z), "jn"),
      yn = function(nu, z) .mul(list(.pow(S_NegativeOne, .add(list(nu, S_One))), f("jn", .sub(.neg(nu), S_One), z))),
      hn1 = function(nu, z) .add(list(f("jn", nu, z), .mul(list(S_I, rewrite(f("yn", nu, z), "jn"))))),
      hn2 = function(nu, z) .add(list(f("jn", nu, z), .mul(list(.neg(S_I), rewrite(f("yn", nu, z), "jn")))))),
    yn = list(
      bessely = function(nu, z) .mul(list(sq(.div(.mul(list(S_Two, z)), S_Pi)), f("yn", .sub(nu, hlf), z))),
      besselk = then(then(then(ktoi, "besselj"), "bessely"), "yn"),
      jn = function(nu, z) .mul(list(.pow(S_NegativeOne, nu), f("yn", .sub(.neg(nu), S_One), z))),
      hn1 = function(nu, z) .add(list(rewrite(f("jn", nu, z), "yn"), .mul(list(S_I, f("yn", nu, z))))),
      hn2 = function(nu, z) .add(list(rewrite(f("jn", nu, z), "yn"), .mul(list(.neg(S_I), f("yn", nu, z)))))),
    hankel1 = list(hn1 = intord(function(nu, z) .mul(list(sqp2z(z), f("hankel1", nu, z))))),
    hankel2 = list(hn2 = intord(function(nu, z) .mul(list(sqp2z(z), f("hankel2", nu, z))))))
  for (tg in names(rules)) {
    cur <- .rewrite_rules[[tg]]
    if (is.null(cur)) cur <- list()
    for (nm in names(rules[[tg]])) if (!is.null(rules[[tg]][[nm]])) cur[[nm]] <- rules[[tg]][[nm]]
    .rewrite_rules[[tg]] <<- cur
  }
})
