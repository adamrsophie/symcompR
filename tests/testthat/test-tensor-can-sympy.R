# Expected values are SymPy's outputs for the same calls (Tensor canonicalization).

skip_unless_full_tests()

skip_on_cran()

.bt <- new.env(parent = asNamespace("symcompR"))
eval(parse(text = r"---(tup <- function(x) structure(x, py_tuple = TRUE)
.perm_pystr <- function(p, repr) {
  a <- p$array; n <- length(a)
  if (!n) return(if (repr) "Permutation()" else "()")
  cf <- cyclic_form(p)
  sep <- if (repr) ", " else " "
  s <- paste0(vapply(cf, function(c) paste0("(", paste(c, collapse = sep), ")"), ""), collapse = "")
  if (!((n - 1L) %in% unlist(cf))) s <- paste0("(", n - 1L, ")", s)
  if (repr) paste0("Permutation", s) else s
}
pystr <- function(x, top = TRUE) {
  if (is.null(x)) return("None")
  if (inherits(x, "Permutation")) return(.perm_pystr(x, !top))
  if (inherits(x, "sym")) return(sstr(x))
  if (inherits(x, "TensExpr")) return(.tn_str(x))
  if (inherits(x, c("TensorIndex", "TensorIndexType", "TensorHead", "TensorSymmetry"))) return(format(x))
  if (inherits(x, "SymArray")) {
    nest <- function(data, shape) {
      if (length(shape) == 1L) return(paste0("[", paste(vapply(data, sstr, ""), collapse = ", "), "]"))
      m <- length(data) / shape[1L]
      paste0("[", paste(vapply(seq_len(shape[1L]), function(k) nest(data[(k - 1L) * m + seq_len(m)], shape[-1L]), ""), collapse = ", "), "]")
    }
    return(if (length(x$shape)) nest(x$data, x$shape) else sstr(x$data[[1L]]))
  }
  if (isTRUE(attr(x, "py_tuple"))) {
    el <- vapply(x, pystr, "", top = FALSE)
    return(paste0("(", paste(el, collapse = ", "), if (length(el) == 1L) "," else "", ")"))
  }
  if (is.list(x)) {
    el <- vapply(x, pystr, "", top = FALSE)
    if (!is.null(names(x))) return(paste0("{", paste(paste0(names(x), ": ", el), collapse = ", "), "}"))
    return(paste0("[", paste(el, collapse = ", "), "]"))
  }
  if (is.integer(x)) return(paste0("[", paste(x, collapse = ", "), "]"))
  if (is.logical(x) && length(x) == 1L) return(if (x) "True" else "False")
  if (is.double(x) && length(x) == 1L) return(format(x, scientific = FALSE))
  if (is.character(x) && length(x) == 1L) return(if (top) x else paste0("'", x, "'"))
  paste(format(x), collapse = " ")
}
.errcls <- function(e) { cl <- setdiff(class(e), c("simpleError", "error", "condition", "sympy_error")); if (length(cl)) cl[1L] else "Error" }
run_battery <- function(cases, out, env, tl = 90) {
  lines <- readLines(cases); lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#")]
  con <- file(out, "w")
  for (ln in lines) {
    rx <- strsplit(ln, " ||| ", fixed = TRUE)[[1L]][2L]
    if (startsWith(ln, "@@")) {
      r <- tryCatch({ eval(parse(text = rx), envir = env); NULL }, error = function(e) conditionMessage(e))
      if (!is.null(r)) cat("SETUP ERROR:", ln, "::", r, "\n")
      next
    }
    o <- tryCatch({ setTimeLimit(cpu = tl, elapsed = 10 * tl, transient = TRUE); r <- pystr(eval(parse(text = rx), envir = env)); setTimeLimit(cpu = Inf, elapsed = Inf); r },
                  error = function(e) { setTimeLimit(cpu = Inf, elapsed = Inf); if (grepl("time limit", conditionMessage(e))) "TIMEOUT" else paste0("ERR ", .errcls(e), ": ", conditionMessage(e)) })
    writeLines(gsub("\n", " ", o), con); flush(con)
  }
  close(con)
}


  P <- function(...) Permutation(c(...))
  s <- get_symmetric_group_sgs(1); base1 <- s[[1]]; gens1 <- s[[2]]
  s <- get_symmetric_group_sgs(1, 1); base1a <- s[[1]]; gens1a <- s[[2]]
  s <- get_symmetric_group_sgs(2); base2 <- s[[1]]; gens2 <- s[[2]]
  s <- get_symmetric_group_sgs(3); base3 <- s[[1]]; gens3 <- s[[2]]
  s <- get_symmetric_group_sgs(2, 1); base2a <- s[[1]]; gens2a <- s[[2]]
  s <- get_symmetric_group_sgs(3, 1); base3a <- s[[1]]; gens3a <- s[[2]]
  baser <- riemann_bsgs[[1]]; gensr <- riemann_bsgs[[2]]
  s <- bsgs_direct_product(base1, gens1, base2a, gens2a); base_f <- s[[1]]; gens_f <- s[[2]]
  s <- bsgs_direct_product(base1, gens1, base1, gens1); base_A <- s[[1]]; gens_A <- s[[2]]
  s <- bsgs_direct_product(base1, gens1, base2, gens2); base_s <- s[[1]]; gens_s <- s[[2]]
  tn <- list(integer(), list(P(0, 1, 2, 3)), 1, 0)
  dgens <- list(P(2,1,0,3,4,5,7,6), P(4,1,2,3,0,5,7,6))
  nonmin <- list(c(2L, 0L), list(P(1,0,2,3,5,4), P(2,3,0,1,4,5)))
)---"), envir = .bt)
.run <- function(code) {
  e <- .bt
  NULL
  gsub("\n", " ", get("pystr", envir = .bt)(eval(parse(text = code), envir = e)), fixed = TRUE)
}
.setup <- function(code) invisible(eval(parse(text = code), envir = .bt))
.err <- function(code) tryCatch({ .run(code); "no error" }, error = function(e) gsub("\n", " ", conditionMessage(e), fixed = TRUE))

test_that("Tensor canonicalization match SymPy", {
  expect_equal(.run("perm_af_direct_product(list(c(1,0,2,3), c(0,1,3,2)), list(c(1,0)), 0)"), "[[1, 0, 2, 3, 4, 5], [0, 1, 3, 2, 4, 5], [0, 1, 2, 3, 5, 4]]")
  expect_equal(.run("perm_af_direct_product(list(c(1,0,2,3,5,4), c(0,1,3,2,4,5)), list(c(1,0,2,3)), TRUE)"), "[[1, 0, 2, 3, 4, 5, 7, 6], [0, 1, 3, 2, 4, 5, 6, 7], [0, 1, 2, 3, 5, 4, 6, 7]]")
  expect_equal(.run("dummy_sgs(c(1,2), 0, 4)"), "[[0, 2, 1, 3, 4, 5]]")
  expect_equal(.run("dummy_sgs(c(2,3,4,5), 0, 8)"), "[[0, 1, 3, 2, 4, 5, 6, 7, 8, 9], [0, 1, 2, 3, 5, 4, 6, 7, 8, 9], [0, 1, 4, 5, 2, 3, 6, 7, 8, 9]]")
  expect_equal(.run("dummy_sgs(c(2,3,4,5), 1, 8)"), "[[0, 1, 3, 2, 4, 5, 6, 7, 9, 8], [0, 1, 2, 3, 5, 4, 6, 7, 9, 8], [0, 1, 4, 5, 2, 3, 6, 7, 8, 9]]")
  expect_equal(.run("dummy_sgs(2:7, 0, 8)"), "[[0, 1, 3, 2, 4, 5, 6, 7, 8, 9], [0, 1, 2, 3, 5, 4, 6, 7, 8, 9], [0, 1, 2, 3, 4, 5, 7, 6, 8, 9], [0, 1, 4, 5, 2, 3, 6, 7, 8, 9], [0, 1, 2, 3, 6, 7, 4, 5, 8, 9]]")
  expect_equal(.run("dummy_sgs(2:7, NULL, 8)"), "[[0, 1, 4, 5, 2, 3, 6, 7, 8, 9], [0, 1, 2, 3, 6, 7, 4, 5, 8, 9]]")
  expect_equal(.err("dummy_sgs(0:9, 0, 8)"), "List too large")
  expect_equal(.run("tup(get_symmetric_group_sgs(0))"), "([], [Permutation(1)])")
  expect_equal(.run("tup(get_symmetric_group_sgs(0, 1))"), "([], [Permutation(1)])")
  expect_equal(.run("tup(get_symmetric_group_sgs(1))"), "([], [Permutation(2)])")
  expect_equal(.run("tup(get_symmetric_group_sgs(1, 1))"), "([], [Permutation(2)])")
  expect_equal(.run("tup(get_symmetric_group_sgs(2))"), "([0], [Permutation(3)(0, 1)])")
  expect_equal(.run("tup(get_symmetric_group_sgs(2, 1))"), "([0], [Permutation(0, 1)(2, 3)])")
  expect_equal(.run("tup(get_symmetric_group_sgs(3))"), "([0, 1], [Permutation(4)(0, 1), Permutation(4)(1, 2)])")
  expect_equal(.run("tup(get_symmetric_group_sgs(3, 1))"), "([0, 1], [Permutation(0, 1)(3, 4), Permutation(1, 2)(3, 4)])")
  expect_equal(.run("tup(get_symmetric_group_sgs(4))"), "([0, 1, 2], [Permutation(5)(0, 1), Permutation(5)(1, 2), Permutation(5)(2, 3)])")
  expect_equal(.run("tup(get_symmetric_group_sgs(4, 1))"), "([0, 1, 2], [Permutation(0, 1)(4, 5), Permutation(1, 2)(4, 5), Permutation(2, 3)(4, 5)])")
  expect_equal(.run("tup(do.call(bsgs_direct_product, c(get_symmetric_group_sgs(0), get_symmetric_group_sgs(0))))"), "([], [Permutation(1)])")
  expect_equal(.run("tup(do.call(bsgs_direct_product, c(get_symmetric_group_sgs(0), get_symmetric_group_sgs(2, FALSE))))"), "([0], [Permutation(3)(0, 1)])")
  expect_equal(.run("tup(do.call(bsgs_direct_product, c(get_symmetric_group_sgs(2, FALSE), get_symmetric_group_sgs(0))))"), "([0], [Permutation(3)(0, 1)])")
  expect_equal(.run("tup(do.call(bsgs_direct_product, c(get_symmetric_group_sgs(0), get_symmetric_group_sgs(2, TRUE))))"), "([0], [Permutation(0, 1)(2, 3)])")
  expect_equal(.run("tup(do.call(bsgs_direct_product, c(get_symmetric_group_sgs(2, TRUE), get_symmetric_group_sgs(0))))"), "([0], [Permutation(0, 1)(2, 3)])")
  expect_equal(.run("tup(bsgs_direct_product(base1, gens1, base2, gens2))"), "([1], [Permutation(4)(1, 2)])")
  expect_equal(.run("tup(bsgs_direct_product(base2a, gens2a, baser, gensr))"), "([0, 2, 4], [Permutation(0, 1)(6, 7), Permutation(2, 3)(6, 7), Permutation(4, 5)(6, 7), Permutation(7)(2, 4)(3, 5)])")
  expect_equal(.run("tup(bsgs_direct_product(base2, gens2, base3a, gens3a, FALSE))"), "([0, 4, 5], [Permutation(8)(0, 1), Permutation(4, 5)(7, 8), Permutation(5, 6)(7, 8)])")
  expect_equal(.run(".tc_min_dummies(list(2:7), list(0), 0:9)"), "[0, 1, 2, 2, 2, 2, 2, 2, 8, 9]")
  expect_equal(.run(".tc_min_dummies(list(2:3, 4:7), list(0, NULL), 0:9)"), "[0, 1, 2, 2, 4, 4, 4, 4, 8, 9]")
  expect_equal(.run("get_transversals(baser, gensr)"), "[{0: [0, 1, 2, 3, 4, 5], 1: [1, 0, 2, 3, 5, 4], 2: [2, 3, 0, 1, 4, 5], 3: [3, 2, 0, 1, 5, 4]}, {2: [0, 1, 2, 3, 4, 5], 3: [0, 1, 3, 2, 5, 4]}]")
  expect_equal(.run("get_transversals(base3, gens3)"), "[{0: [0, 1, 2, 3, 4], 1: [1, 0, 2, 3, 4], 2: [2, 0, 1, 3, 4]}, {1: [0, 1, 2, 3, 4], 2: [0, 2, 1, 3, 4]}]")
  expect_equal(.run("get_transversals(base_f, gens_f)"), "[{1: [0, 1, 2, 3, 4], 2: [0, 2, 1, 4, 3]}]")
  expect_equal(.run("get_transversals(integer(), gens2)"), "[]")
  expect_equal(.run("transversal2coset(6, baser, get_transversals(baser, gensr))"), "[[[0, 1, 2, 3, 4, 5], [1, 0, 2, 3, 5, 4], [2, 3, 0, 1, 4, 5], [3, 2, 0, 1, 5, 4]], [[0, 1, 2, 3, 4, 5]], [[0, 1, 2, 3, 4, 5], [0, 1, 3, 2, 5, 4]]]")
  expect_equal(.run("transversal2coset(5, base3, get_transversals(base3, gens3))"), "[[[0, 1, 2, 3, 4], [1, 0, 2, 3, 4], [2, 0, 1, 3, 4]], [[0, 1, 2, 3, 4], [0, 2, 1, 3, 4]]]")
  expect_equal(.run("double_coset_can_rep(list(0:5), list(0), c(0, 2), dgens, get_transversals(c(0, 2), dgens), P(4, 2, 0, 1, 3, 5, 6, 7))"), "[0, 1, 2, 3, 4, 5, 7, 6]")
  expect_equal(.run("double_coset_can_rep(list(0:5), list(0), c(0, 2), dgens, get_transversals(c(0, 2), dgens), P(4, 1, 3, 0, 5, 2, 6, 7))"), "0")
  expect_equal(.run("double_coset_can_rep(list(0:5), list(1), c(0, 2), dgens, get_transversals(c(0, 2), dgens), P(4, 2, 0, 1, 3, 5, 6, 7))"), "[0, 1, 2, 3, 4, 5, 6, 7]")
  expect_equal(.run("double_coset_can_rep(list(0:5), list(NULL), c(0, 2), dgens, get_transversals(c(0, 2), dgens), P(4, 2, 0, 1, 3, 5, 6, 7))"), "[0, 2, 4, 1, 3, 5, 7, 6]")
  expect_equal(.run("canonical_free(c(0, 2), list(P(1,0,2,3,5,4), P(2,3,0,1,4,5), P(0,1,3,2,5,4)), P(2,1,0,3,4,5), 4)"), "[0, 3, 1, 2, 5, 4]")
  expect_equal(.run("{tg <- tensor_gens(baser, gensr, list(integer(), integer()), 0); canonical_free(tg[[2]], tg[[3]], P(0,3,4,6,7,5,2,1,8,9), 2)}"), "[0, 3, 4, 6, 1, 2, 7, 5, 9, 8]")
  expect_equal(.run("canonical_free(integer(), gens2, P(1,0,2,3), 2)"), "[1, 0, 2, 3]")
  expect_equal(.run(".tc_get_map_slots(8L, c(1L, 4L))"), "[0, 1, 1, 2, 4, 3, 4, 5]")
  expect_equal(.run(".tc_lift_sgens(8L, c(5L, 1L), c(0L, 1L), c(0L, 2L, 1L, 3L, 5L, 4L))"), "[2, 1, 4, 3, 5, 0, 7, 6]")
  expect_equal(.run("tup(tensor_gens(base3, gens3, list(integer(), integer())))"), "(8, [0, 1, 3, 4], [Permutation(7)(0, 1), Permutation(7)(1, 2), Permutation(7)(3, 4), Permutation(7)(4, 5), Permutation(7)(0, 3)(1, 4)(2, 5)])")
  expect_equal(.run("tup(tensor_gens(base3, gens3, list(1, 0)))"), "(8, [0, 4], [Permutation(7)(0, 2), Permutation(7)(4, 5)])")
  expect_equal(.run("tup(tensor_gens(base3, gens3, list(integer(), integer(), 0, c(2, 1))))"), "(14, [0, 1, 3, 4, 7], [Permutation(13)(0, 1), Permutation(13)(1, 2), Permutation(13)(3, 4), Permutation(13)(4, 5), Permutation(13)(7, 8), Permutation(13)(0, 3)(1, 4)(2, 5)])")
  expect_equal(.run("tup(tensor_gens(base1, gens1, list(integer(), integer())))"), "(4, [0], [Permutation(3), Permutation(3)(0, 1)])")
  expect_equal(.run("tup(tensor_gens(base1, gens1, list(integer(), integer(), integer()), 1))"), "(5, [0, 1], [Permutation(4), Permutation(0, 1)(3, 4), Permutation(1, 2)(3, 4)])")
  expect_equal(.run("tup(tensor_gens(baser, gensr, list(integer(), integer(), 0), 0))"), "(14, [0, 2, 4, 6, 10], [Permutation(0, 1)(12, 13), Permutation(2, 3)(12, 13), Permutation(13)(0, 2)(1, 3), Permutation(4, 5)(12, 13), Permutation(6, 7)(12, 13), Permutation(13)(4, 6)(5, 7), Permutation(10, 11)(12, 13), Permutation(13)(0, 4)(1, 5)(2, 6)(3, 7)])")
  expect_equal(.run("tup(tensor_gens(baser, gensr, list(c(1, 3), integer()), 1))"), "(10, [4, 6], [Permutation(4, 5)(8, 9), Permutation(6, 7)(8, 9), Permutation(9)(4, 6)(5, 7)])")
  expect_equal(.run("tup(tensor_gens(base2, gens2, list(integer(), integer()), NULL))"), "(6, [0, 2], [Permutation(5)(0, 1), Permutation(5)(2, 3)])")
  expect_equal(.run("tup(tensor_gens(integer(), list(P(0,1,2,3)), list(integer(), integer()), 0))"), "(6, [0], [Permutation(5), Permutation(5)(0, 2)(1, 3)])")
  expect_equal(.run("tup(tensor_gens(integer(), list(P(0,1,2,3)), list(integer(), 1), 0))"), "(6, [], [Permutation(5)])")
  expect_equal(.run("tup(tensor_gens(base2a, gens2a, list(0, 1), 0))"), "(6, [], [Permutation(5)])")
  expect_equal(.run("tup(gens_products(list(base2, gens2, list(integer(), integer()), 0)))"), "(6, [0, 2], [Permutation(5)(0, 1), Permutation(5)(2, 3), Permutation(5)(0, 2)(1, 3)])")
  expect_equal(.run("tup(gens_products(list(base2, gens2, list(1, integer()), 0)))"), "(6, [2], [Permutation(5)(2, 3)])")
  expect_equal(.run("tup(gens_products(list(integer(), list(P(0,1,2,3)), list(integer()), 0), list(integer(), list(P(0,1,2,3)), list(integer()), 0)))"), "(6, [], [Permutation(5)])")
  expect_equal(.run("tup(gens_products(list(baser, gensr, list(integer(), integer()), 0), list(base2, gens2, list(integer(), integer(), integer()), 0)))"), "(16, [0, 2, 4, 6, 8, 10, 12], [Permutation(0, 1)(14, 15), Permutation(2, 3)(14, 15), Permutation(15)(0, 2)(1, 3), Permutation(4, 5)(14, 15), Permutation(6, 7)(14, 15), Permutation(15)(4, 6)(5, 7), Permutation(15)(0, 4)(1, 5)(2, 6)(3, 7), Permutation(15)(8, 9), Permutation(15)(10, 11), Permutation(15)(12, 13), Permutation(15)(8, 10)(9, 11), Permutation(15)(10, 12)(11, 13)])")
  expect_equal(.run("tup(gens_products(list(base3a, gens3a, list(1), 1), list(base1, gens1, list(integer(), integer()), NULL)))"), "(7, [0], [Permutation(0, 2)(5, 6)])")
  expect_equal(.run("tup(get_minimal_bsgs(c(2, 0), list(P(1,0,2,3,5,4), P(2,3,0,1,4,5))))"), "([0, 2], [Permutation(0, 1)(4, 5), Permutation(5)(0, 2)(1, 3), Permutation(2, 3)(4, 5)])")
  expect_equal(.run("tup(get_minimal_bsgs(baser, gensr))"), "([0, 2], [Permutation(0, 1)(4, 5), Permutation(2, 3)(4, 5), Permutation(5)(0, 2)(1, 3)])")
  expect_equal(.run("tup(get_minimal_bsgs(base3, gens3))"), "([0, 1], [Permutation(4)(0, 1), Permutation(4)(1, 2)])")
  expect_equal(.run("get_minimal_bsgs(c(1, 0), list(P(0,2,1,3,4), P(1,0,2,3,4)))"), "None")
  expect_equal(.run("tup(get_minimal_bsgs(integer(), list(P(1,2,0,3,4,5,6), P(0,1,2,4,5,3,6))))"), "([0, 3], [Permutation(6)(0, 1, 2), Permutation(6)(3, 4, 5)])")
  expect_equal(.run("tup(get_minimal_bsgs(integer(), list(P(1,2,3,4,0,5,6), P(4,3,2,1,0,6,5))))"), "([0, 1], [Permutation(6)(0, 1, 2, 3, 4), Permutation(0, 4)(1, 3)(5, 6), Permutation(1, 4)(2, 3)(5, 6)])")
  expect_equal(.run("tup(get_minimal_bsgs(integer(), list(P(0,1,2,3))))"), "([], [Permutation(3)])")
  expect_equal(.run(".tc_is_minimal_bsgs(baser, lapply(gensr, .tc_af))"), "True")
  expect_equal(.run(".tc_is_minimal_bsgs(c(2L, 0L), list(c(1L,0L,2L,3L,5L,4L), c(2L,3L,0L,1L,4L,5L)))"), "False")
  expect_equal(.run("canonicalize(P(1,0,2,3), c(0, 1), 0, list(base1,gens1,1,0), list(base1,gens1,1,0))"), "[0, 1, 2, 3]")
  expect_equal(.run("canonicalize(P(1,0,2,3), c(0, 1), 0, list(base1, gens1, 2, NULL))"), "[0, 1, 2, 3]")
  expect_equal(.run("canonicalize(P(1,0,2,3), c(0, 1), 1, list(base1,gens1,1,0), list(base1,gens1,1,0))"), "[0, 1, 3, 2]")
  expect_equal(.run("canonicalize(P(0,1,2,3), integer(), 0, list(base1,gens1,1,0), list(base1,gens1,1,0))"), "[0, 1, 2, 3]")
  expect_equal(.run("canonicalize(P(1,0,2,3), integer(), 0, list(base1,gens1,1,0), list(base1,gens1,1,0))"), "[1, 0, 2, 3]")
  expect_equal(.run("canonicalize(P(1,3,2,0,4,5), c(2,3), 0, list(base2, gens2, 2, 0))"), "[0, 2, 1, 3, 4, 5]")
  expect_equal(.run("canonicalize(P(1,3,2,0,4,5), c(2,3), 1, list(base2, gens2, 2, 0))"), "[0, 2, 1, 3, 4, 5]")
  expect_equal(.run("canonicalize(P(0,3,2,1,4,5), c(2,3), 1, list(base2, gens2, 2, 0))"), "[0, 2, 1, 3, 5, 4]")
  expect_equal(.run("canonicalize(P(1,3,2,0,4,5), c(2,3), 0, list(base2,gens2,1,0), list(base2,gens2,1,0))"), "[1, 2, 0, 3, 4, 5]")
  expect_equal(.run("canonicalize(P(1,3,2,0,4,5), c(2,3), 1, list(base2,gens2,1,0), list(base2,gens2,1,0))"), "[1, 2, 0, 3, 5, 4]")
  expect_equal(.run("canonicalize(P(2,1,0,3,4,5), 0:3, 0, list(base2, gens2, 1, 0), list(base1, gens1, 1, 0), list(base1, gens1, 1, 0))"), "[0, 2, 1, 3, 4, 5]")
  expect_equal(.run("canonicalize(P(2,1,0,3,4,5), 0:3, 0, tn, list(base1,gens1,1,0), list(base1,gens1,1,0))"), "[0, 2, 3, 1, 4, 5]")
  expect_equal(.run("canonicalize(P(2,1,3,0,4,5), 0:3, 0, tn, tn)"), "[0, 2, 1, 3, 4, 5]")
  expect_equal(.run("canonicalize(P(1,2,3,0,4,5), 0:3, 0, tn, tn)"), "[0, 2, 3, 1, 4, 5]")
  expect_equal(.run("canonicalize(P(4,2,0,3,5,1,6,7), 2:5, 0, tn, tn, tn)"), "[2, 4, 0, 5, 3, 1, 6, 7]")
  expect_equal(.run("canonicalize(P(4,2,0,3,5,1,6,7), 2:5, 0, list(base2,gens2,1,0), tn, tn)"), "[2, 4, 0, 3, 5, 1, 6, 7]")
  expect_equal(.run("canonicalize(P(4,2,0,3,5,1,6,7), 2:5, 0, list(base2,gens2,1,0), tn, list(base2,gens2,1,0))"), "[2, 4, 0, 3, 1, 5, 6, 7]")
  expect_equal(.run("canonicalize(P(4,2,0,3,5,1,6,7), 2:5, 0, list(base2,gens2,1,0), tn, list(base2a,gens2a,1,0))"), "[2, 4, 0, 3, 1, 5, 7, 6]")
  expect_equal(.run("canonicalize(P(2,1,0,3,4), integer(), 0, list(base1, gens1, 3, 0))"), "[0, 1, 2, 3, 4]")
  expect_equal(.run("canonicalize(P(2,1,0,3,4), integer(), 0, list(base1, gens1, 3, 1))"), "[0, 1, 2, 4, 3]")
  expect_equal(.run("canonicalize(P(1,3,2,0,4,5), integer(), 0, list(base2, gens2, 2, 0))"), "[0, 2, 1, 3, 4, 5]")
  expect_equal(.run("canonicalize(P(1,3,2,0,4,5), integer(), 0, list(base2, gens2, 2, 1))"), "[0, 2, 1, 3, 5, 4]")
  expect_equal(.run("canonicalize(P(2,0,1,3,4,5), integer(), 0, list(base2, gens2, 2, 1))"), "[0, 2, 1, 3, 4, 5]")
  expect_equal(.run("canonicalize(P(2,1,0,3,4,5), 0:3, NULL, list(integer(), list(P(0,1,2,3)), 2, 0))"), "[0, 3, 2, 1, 4, 5]")
  expect_equal(.run("canonicalize(P(2,5,0,7,4,3,6,1,8,9), 0:7, NULL, list(integer(), list(P(0,1,2,3)), 4, 0))"), "[0, 3, 2, 1, 4, 7, 6, 5, 8, 9]")
  expect_equal(.run("canonicalize(P(0,5,2,7,6,1,4,3,8,9), 0:7, NULL, list(integer(), list(P(0,1,2,3)), 4, 0))"), "[0, 3, 2, 5, 4, 7, 6, 1, 8, 9]")
  expect_equal(.run("canonicalize(P(12,7,10,3,14,13,4,11,6,1,2,9,0,15,8,5,16,17), 0:15, NULL, list(integer(), list(P(0,1,2,3)), 8, 0))"), "[0, 3, 2, 5, 4, 7, 6, 1, 8, 11, 10, 13, 12, 15, 14, 9, 16, 17]")
  expect_equal(.run("canonicalize(P(2,1,3,0,4,5), list(c(2,3)), list(NULL), list(integer(), list(P(0,1,2,3)), 2, 0))"), "[3, 0, 2, 1, 4, 5]")
  expect_equal(.run("canonicalize(P(1,0,2,3), c(0, 1), 0, list(base1, gens1, 2, 0))"), "[0, 1, 2, 3]")
  expect_equal(.run("canonicalize(P(1,3,5,4,2,0,6,7), 0:5, 0, list(base1, gens1, 6, 0))"), "[0, 1, 2, 3, 4, 5, 6, 7]")
  expect_equal(.run("canonicalize(P(1,3,5,4,2,0,6,7), 0:5, 0, list(base1, gens1, 6, 1))"), "0")
  expect_equal(.run(".tc_canonicalize_naive(P(1,3,5,4,2,0,6,7), 0:5, 0, list(base1, gens1, 6, 1))"), "0")
  expect_equal(.run("canonicalize(P(2,1,0,5,4,3,6,7), 2:5, 0, list(base2, gens2, 3, 0))"), "[0, 2, 1, 4, 3, 5, 6, 7]")
  expect_equal(.run("canonicalize(P(2,1,4,3,0,5,6,7), 2:5, 0, list(base2,gens2,2,0), list(base2,gens2,1,0))"), "[1, 2, 3, 4, 0, 5, 6, 7]")
  expect_equal(.run("canonicalize(P(4,2,1,0,5,3,6,7), 2:5, 0, list(base3, gens3, 2, 0))"), "[0, 2, 4, 1, 3, 5, 6, 7]")
  expect_equal(.run("canonicalize(P(10,4,8, 0,7,9, 6,11,1, 2,3,5, 12,13), 4:11, 0, list(base3, gens3, 4, 0))"), "[0, 4, 6, 1, 5, 8, 2, 3, 10, 7, 9, 11, 12, 13]")
  expect_equal(.run("canonicalize(P(0,2,4,5,7,3,1,6,8,9), 0:7, 0, list(base3, gens3,2,0), list(base2a,gens2a,1,0))"), "0")
  expect_equal(.run("canonicalize(P(0,2,4,5,7,3,1,6,8,9), 0:7, 0, list(base3, gens3,2,1), list(base2a,gens2a,1,0))"), "[0, 2, 4, 1, 3, 6, 5, 7, 8, 9]")
  expect_equal(.run("canonicalize(P(0,2,4,5,7,3,1,6,8,9), 0:7, 1, list(base3, gens3,2,1), list(base2a,gens2a,1,0))"), "[0, 2, 4, 1, 3, 6, 5, 7, 9, 8]")
  expect_equal(.run("canonicalize(P(0,2,4,5,7,3,1,6,8,9), 0:7, NULL, list(base3, gens3,2,1), list(base2a,gens2a,1,0))"), "[0, 2, 4, 1, 3, 7, 5, 6, 8, 9]")
  expect_equal(.run("canonicalize(P(3,5,1,4,2,0,6,7), 2:5, 0, list(base2a, gens2a, 1, NULL), list(base1, gens1, 1, NULL), list(base3a, gens3a, 1, NULL))"), "[2, 4, 1, 0, 3, 5, 7, 6]")
  expect_equal(.run("canonicalize(P(5,7,2,1,3,6,4,0,8,9), 4:7, 0, list(base2a, gens2a, 2, NULL), list(base1, gens1, 1, NULL), list(base3a, gens3a, 1, NULL))"), "[4, 6, 1, 2, 3, 0, 5, 7, 8, 9]")
  expect_equal(.run("canonicalize(P(8,11,5, 9,13,7, 1,10, 3,4, 2,12, 0,6, 14,15), list(0:3, 4:13), list(0, 0), list(base_f, gens_f, 2, 0), list(base_A, gens_A, 4, 0))"), "[4, 6, 8, 5, 10, 12, 0, 7, 1, 11, 2, 9, 3, 13, 15, 14]")
  expect_equal(.run("canonicalize(P(1,3,0,5,4,2,6,7), 0:5, 0, list(base2a, gens2a, 1, 0), list(base2a, gens2a, 2, 0))"), "0")
  expect_equal(.run("canonicalize(P(1,3,0,5,4,2,6,7), 0:5, 0, list(base2a, gens2a, 1, 0), list(base2a, gens2a, 2, 1))"), "[0, 2, 1, 4, 3, 5, 7, 6]")
  expect_equal(.run("canonicalize(P(0,7,3,1,9,5,11,6,10,4,13,2,12,8,14,15), list(2:9, 10:13), list(0, 0), list(base3a, gens3a, 2, 0), list(base_A, gens_A, 4, 0))"), "[0, 2, 4, 1, 6, 8, 10, 3, 11, 7, 12, 5, 13, 9, 15, 14]")
  expect_equal(.run("canonicalize(P(0,2,3,1,4,5), 2:3, 0, list(baser, gensr, 1, 0))"), "[0, 2, 1, 3, 5, 4]")
  expect_equal(.run("canonicalize(P(0,2,3,1,4,5), 2:3, 0, list(nonmin[[1]], nonmin[[2]], 1, 0))"), "[0, 2, 1, 3, 5, 4]")
  expect_equal(.run("canonicalize(P(23,2,1,10,12,8,0,11,15,5,17,19,21,7,13,9,4,14,22,3,16,18,6,20,24,25), 0:23, 0, list(baser, gensr, 6, 0))"), "[0, 2, 4, 6, 1, 3, 8, 10, 5, 7, 12, 14, 9, 11, 16, 18, 13, 15, 20, 22, 17, 19, 21, 23, 24, 25]")
  expect_equal(.run("canonicalize(P(23,2,1,10,12,8,0,11,15,5,17,19,21,7,13,9,4,14,22,3,16,18,6,20,24,25), 0:23, 0, list(nonmin[[1]], nonmin[[2]], 6, 0))"), "[0, 2, 4, 6, 1, 3, 8, 10, 5, 7, 12, 14, 9, 11, 16, 18, 13, 15, 20, 22, 17, 19, 21, 23, 24, 25]")
  expect_equal(.run("canonicalize(P(0,2,5,7,4,6,9,11,8,10,13,15,12,14,17,19,16,18,21,23,20,22,25,27,24,26,29,31,28,30,33,35,32,34,37,39,36,38,1,3,40,41), 0:39, 0, list(baser, gensr, 10, 0))"), "[0, 2, 4, 6, 1, 3, 8, 10, 5, 7, 12, 14, 9, 11, 16, 18, 13, 15, 20, 22, 17, 19, 24, 26, 21, 23, 28, 30, 25, 27, 32, 34, 29, 31, 36, 38, 33, 35, 37, 39, 40, 41]")
  expect_equal(.run("canonicalize(P(0,1,2,3,4,5), 2:3, 0, list(baser, gensr, 1, 0))"), "0")
  expect_equal(.run("canonicalize(P(2,1,0,3,4,5), 2:3, 0, list(baser, gensr, 1, 0))"), "[0, 2, 1, 3, 5, 4]")
  expect_equal(.run("canonicalize(P(4,7,1,3,2,0,5,6,8,9), 2:7, 0, list(baser, gensr, 2, 0))"), "[0, 2, 4, 6, 1, 3, 5, 7, 9, 8]")
  expect_equal(.run(".tc_canonicalize_naive(P(4,7,1,3,2,0,5,6,8,9), 2:7, 0, list(baser, gensr, 2, 0))"), "[0, 2, 4, 6, 1, 3, 5, 7, 9, 8]")
  expect_equal(.run("canonicalize(P(12,10,5,2,8,0,4,6,13,1,7,3,9,11,14,15), 0:13, 0, list(baser,gensr,2,0), list(base2,gens2,3,0))"), "[0, 2, 4, 6, 1, 8, 10, 12, 3, 9, 5, 11, 7, 13, 15, 14]")
  expect_equal(.run("canonicalize(P(10,0,2,6,8,11,1,3,4,5,7,9,12,13), 6:11, 0, list(baser, gensr, 3, 0))"), "[0, 6, 2, 8, 1, 3, 7, 10, 4, 5, 9, 11, 12, 13]")
  expect_equal(.run("canonicalize(P(0,4,3,1,2,5,6,7), 4:5, 0, list(base_f, gens_f, 2, 0))"), "[0, 3, 4, 1, 2, 5, 7, 6]")
  expect_equal(.run("canonicalize(P(4,2,10, 0,11,8, 1,9,6, 5,7,3, 12,13), 4:11, 0, list(base_s, gens_s, 4, 0))"), "[0, 4, 6, 1, 5, 8, 10, 2, 7, 11, 3, 9, 12, 13]")
  expect_equal(.run("canonicalize(P(4,2,10, 0,11,8, 1,9,6, 5,7,3, 12,13), list(4:5, 6:11), list(0, 0), list(base_s, gens_s, 4, 0))"), "[0, 6, 8, 1, 7, 10, 4, 2, 9, 5, 3, 11, 12, 13]")
  expect_equal(.run("canonicalize(P(4,2,10, 0,11,8, 1,9,6, 5,7,3, 12,13), list(4:5, 6:11), list(0, 1), list(base_s, gens_s, 4, 0))"), "[0, 6, 8, 1, 7, 10, 4, 2, 9, 5, 3, 11, 13, 12]")
  expect_equal(.run("canonicalize(P(4,2,10, 0,11,8, 1,9,6, 5,7,3, 12,13), list(4:5, 6:11), list(NULL, 1), list(base_s, gens_s, 4, 0))"), "[0, 6, 8, 1, 7, 10, 4, 2, 9, 5, 3, 11, 13, 12]")
  expect_equal(.err("canonicalize(P(1,0,2,3), c(0,1), 2, list(base1,gens1,2,0))"), "msym must be 0, 1 or None")
  expect_equal(.err("canonicalize(P(1,0,2,3), list(c(0,1)), list(0, 2), list(base1,gens1,2,0))"), "msym entries must be 0, 1 or None")
  expect_equal(.err("canonicalize(P(1,0,2,3), list(c(0,1)), list(0, 0), list(base1,gens1,2,0))"), "dummies and msym must have the same number of elements")
  expect_equal(.err("canonicalize(P(1,0,2,3), c(1,0), 0, list(base1,gens1,2,0))"), "dummies is not valid")
  expect_equal(.err("canonicalize(P(1,0,2,3,4), c(0,1), 0, list(base1,gens1,2,0))"), "g has size 5, generators have size 4")
  expect_equal(.err("canonicalize(P(1,0,2,3,4), integer(), 0, list(c(1, 0), list(P(0,2,1,3,4), P(1,0,2,3,4)), 1, 0))"), "list index out of range")
  expect_equal(.run("canonicalize(P(2,1,0,3,4), c(1, 2), 0, list(c(1, 0), list(P(0,2,1,3,4), P(1,0,2,3,4)), 1, 0))"), "[0, 1, 2, 3, 4]")
  expect_equal(.run("canonicalize(P(2,1,0,3,4), c(1, 2), 1, list(c(1, 0), list(P(0,2,1,4,3), P(1,0,2,4,3)), 1, 0))"), "[0, 1, 2, 4, 3]")
  expect_equal(.run("canonicalize(P(3,0,2,1,5,4,6,7), 2:5, 0, list(c(1, 0), list(P(0,2,1,3,4), P(1,0,2,3,4)), 1, 0), list(base3, gens3, 1, 0))"), "[0, 2, 3, 1, 4, 5, 6, 7]")
})
