# Cables under point or distributed loads, after
# sympy/physics/continuum_mechanics/cable.py.

#' Cables
#'
#' `Cable(support_1, support_2)` is a cable hung between two supports, each
#' given as `list(label, x, y)`. `$apply_load(-1, list(label, x, y,
#' magnitude, direction))` adds a point load with the direction in degrees
#' and `$apply_load(0, list(label, magnitude))` a load distributed along the
#' horizontal; `$remove_loads()` removes loads by label. `$solve()` finds
#' the tension and the reaction loads; with a distributed load it needs the
#' abscissa of the lowest point of the cable. `$tension_at(x)` evaluates the
#' tension of a cable under a distributed load. The properties are
#' `$supports`, `$left_support`, `$right_support`, `$loads`,
#' `$loads_position`, `$length`, `$reaction_loads` and `$tension`.
#' `$change_support()` and `$apply_length()` change the supports and the
#' length, and `$draw()` and `$plot_tension()` draw base R plots. After
#' `sympy.physics.continuum_mechanics.cable`.
#'
#' @param support_1,support_2 The supports, as `list(label, x, y)`.
#' @return A `Cable` object.
#' @examples
#' c <- Cable(list("A", 0, 10), list("B", 10, 10))
#' c$apply_load(-1, list("X", 5, 5, 10, 270))
#' c$solve()
#' c$reaction_loads
#' @export
Cable <- function(support_1, support_2) {
  self <- .pyenv_new("Cable")
  for (f in c("_left_support", "_right_support", "_support_labels")) .cm_set(self, f, list())
  for (f in c("_supports", "_loads_position", "_reaction_loads", "_tension")) .cm_set(self, f, .pd_empty())
  .cm_set(self, "_loads", .pydict(list("distributed", "point_load"), list(.pd_empty(), .pd_empty())))
  .cm_set(self, "_length", 0L)
  .cm_set(self, "_lowest_x_global", S_Zero)
  .cm_set(self, "_lowest_y_global", S_Zero)
  .cm_set(self, "_cable_eqn", NULL)
  .cm_set(self, "_tension_func", NULL)
  if (.cm_same_key(support_1[[1L]], support_2[[1L]])) .physics_error("ValueError", "Supports can not have the same label")
  if (.cm_same(support_1[[2L]], support_2[[2L]])) .physics_error("ValueError", "Supports can not be at the same location")
  x1 <- as_sym(support_1[[2L]]); y1 <- as_sym(support_1[[3L]])
  x2 <- as_sym(support_2[[2L]]); y2 <- as_sym(support_2[[3L]])
  sup <- .pd_set(.pd_empty(), support_1[[1L]], list(x1, y1))
  .cm_set(self, "_supports", .pd_set(sup, support_2[[1L]], list(x2, y2)))
  if (isTRUE(x1 < x2)) {
    .cm_set(self, "_left_support", list(x1, y1)); .cm_set(self, "_right_support", list(x2, y2))
    .cm_set(self, "_support_labels", list(support_1[[1L]], support_2[[1L]]))
  } else {
    .cm_set(self, "_left_support", list(x2, y2)); .cm_set(self, "_right_support", list(x1, y1))
    .cm_set(self, "_support_labels", list(support_2[[1L]], support_1[[1L]]))
  }
  .cable_zero_reactions(self)
  self
}

.cable_zero_reactions <- function(self) {
  rl <- .cm_get(self, "_reaction_loads")
  for (i in .cm_get(self, "_support_labels")) {
    rl <- .pd_set(rl, Symbol(paste0("R_", i, "_x")), 0L)
    rl <- .pd_set(rl, Symbol(paste0("R_", i, "_y")), 0L)
  }
  .cm_set(self, "_reaction_loads", rl)
}

.cable_loads <- function(self, kind) .pd_get(.cm_get(self, "_loads"), kind)
.cable_set_loads <- function(self, kind, d) .cm_set(self, "_loads", .pd_set(.cm_get(self, "_loads"), kind, d))

.cable_tension_at <- function(self, x) {
  tn <- .cm_get(self, "_tension")
  if (!.pd_has(tn, "distributed")) .physics_error("ValueError", "No distributed load added or solve method not called")
  x <- as_sym(x)
  if (isTRUE(x > .cm_get(self, "_right_support")[[1L]]) || isTRUE(x < .cm_get(self, "_left_support")[[1L]]))
    .physics_error("ValueError", "The value of x should be between the two supports")
  subs(.pd_get(tn, "distributed"), Symbol("X"), x - .cm_get(self, "_lowest_x_global"))
}

.cable_apply_length <- function(self, length) {
  l <- .cm_get(self, "_left_support"); r <- .cm_get(self, "_right_support")
  dist <- ((l[[1L]] - r[[1L]])^2 - (l[[2L]] - r[[2L]])^2)^as_sym(0.5)
  if (isTRUE(as_sym(length) < dist))
    .physics_error("ValueError", "length should not be less than the distance between the supports")
  .cm_set(self, "_length", length)
  invisible(NULL)
}

.cable_change_support <- function(self, label, new_support) {
  sup <- .cm_get(self, "_supports")
  if (!.pd_has(sup, label)) .physics_error("ValueError", "No support exists with the given label")
  labels <- .cm_get(self, "_support_labels")
  i <- .pl_index(labels, label)
  rem <- labels[[i %% 2L + 1L]]
  x1 <- .pd_get(sup, rem)[[1L]]; y1 <- .pd_get(sup, rem)[[2L]]
  x <- as_sym(new_support[[2L]]); y <- as_sym(new_support[[3L]])
  # SymPy compares the first character of each load label with the abscissae
  if (length(.cm_get(self, "_loads_position")$keys))
    .physics_error("TypeError", "cannot determine truth value of Relational")
  .cm_set(self, "_supports", .pd_pop(sup, label))
  .cm_set(self, "_reaction_loads", .pd_empty())
  labels <- .pl_remove(labels, label)
  .cm_set(self, "_supports", .pd_set(.cm_get(self, "_supports"), new_support[[1L]], list(x, y)))
  if (isTRUE(x1 < x)) {
    .cm_set(self, "_left_support", list(x1, y1)); .cm_set(self, "_right_support", list(x, y))
    labels <- c(labels, list(new_support[[1L]]))
  } else {
    .cm_set(self, "_left_support", list(x, y)); .cm_set(self, "_right_support", list(x1, y1))
    labels <- c(list(new_support[[1L]]), labels)
  }
  .cm_set(self, "_support_labels", labels)
  .cable_zero_reactions(self)
  invisible(NULL)
}

.cable_apply_load <- function(self, order, load) {
  if (isTRUE(as_sym(order) == -1L) || identical(as.numeric(order), -1)) {
    if (length(.cable_loads(self, "distributed")$keys)) .physics_error("ValueError", "Distributed load already exists")
    label <- load[[1L]]
    pl <- .cable_loads(self, "point_load")
    if (.pd_has(pl, label)) .physics_error("ValueError", "Label already exists")
    x <- as_sym(load[[2L]]); y <- as_sym(load[[3L]])
    if (isTRUE(x > .cm_get(self, "_right_support")[[1L]]) || isTRUE(x < .cm_get(self, "_left_support")[[1L]]))
      .physics_error("ValueError", "The load should be positioned between the supports")
    .cable_set_loads(self, "point_load", .pd_set(pl, label, list(as_sym(load[[4L]]), as_sym(load[[5L]]))))
    .cm_set(self, "_loads_position", .pd_set(.cm_get(self, "_loads_position"), label, list(x, y)))
  } else if (identical(as.numeric(order), 0)) {
    if (length(.cm_get(self, "_loads_position")$keys)) .physics_error("ValueError", "Point load(s) already exist")
    label <- load[[1L]]
    dl <- .cable_loads(self, "distributed")
    if (.pd_has(dl, label)) .physics_error("ValueError", "Label already exists")
    .cable_set_loads(self, "distributed", .pd_set(dl, label, as_sym(load[[2L]])))
  } else {
    .physics_error("ValueError", "Order should be either -1 or 0")
  }
  invisible(NULL)
}

.cable_remove_loads <- function(self, ...) {
  for (i in list(...)) {
    if (!length(.cm_get(self, "_loads_position")$keys)) {
      if (!.pd_has(.cable_loads(self, "distributed"), i))
        .physics_error("ValueError", paste0("Error removing load ", i, ": no such load exists"))
      # SymPy looks the load up under a misspelled key
      .physics_error("KeyError", "'disrtibuted'")
    }
    pl <- .cable_loads(self, "point_load")
    if (!.pd_has(pl, i)) .physics_error("ValueError", paste0("Error removing load ", i, ": no such load exists"))
    .cable_set_loads(self, "point_load", .pd_pop(pl, i))
    .cm_set(self, "_loads_position", .pd_pop(.cm_get(self, "_loads_position"), i))
  }
  invisible(NULL)
}

# the first character of a support label: SymPy indexes the label string
.cable_key <- function(e) if (is.character(e)) substr(e, 1L, 1L) else e[[1L]]

.cable_solve_point <- function(self) {
  lp <- .cm_get(self, "_loads_position"); pl <- .cable_loads(self, "point_load")
  L <- .cm_get(self, "_left_support"); R <- .cm_get(self, "_right_support")
  labels <- .cm_get(self, "_support_labels")
  o <- order(vapply(lp$values, function(p) as.numeric(N(p[[1L]])), 0))
  sp <- c(list(labels[[1L]]), lapply(o, function(k) list(lp$keys[[k]], lp$values[[k]])), list(labels[[2L]]))
  pos <- function(k) .pd_get(lp, sp[[k]][[1L]])
  mag <- function(k) .pd_get(pl, sp[[k]][[1L]])[[1L]]
  ang <- function(k) S_Pi * .pd_get(pl, sp[[k]][[1L]])[[2L]] / 180L
  tension <- .pd_empty()
  mom_l <- S_Zero; mom_r <- S_Zero; F_x <- S_Zero; F_y <- S_Zero
  len <- S_Zero
  tfun <- list()
  x <- Symbol("x")
  n <- length(sp)
  for (i in seq(2L, n - 1L)) {
    p <- pos(i)
    if (i == 2L) len <- len + sqrt((L[[1L]] - p[[1L]])^2 + (L[[2L]] - p[[2L]])^2)
    else len <- len + sqrt((pos(i - 1L)[[1L]] - p[[1L]])^2 + (pos(i - 1L)[[2L]] - p[[2L]])^2)
    if (i == n - 1L) len <- len + sqrt((R[[1L]] - p[[1L]])^2 + (R[[2L]] - p[[2L]])^2)
    mom_l <- mom_l + mag(i) * cos(ang(i)) * Abs(L[[2L]] - p[[2L]])
    mom_l <- mom_l + mag(i) * sin(ang(i)) * Abs(L[[1L]] - p[[1L]])
    F_x <- F_x + mag(i) * cos(ang(i))
    F_y <- F_y + mag(i) * sin(ang(i))
    label <- Symbol(paste0(.cable_key(sp[[i]]), "_", .cable_key(sp[[i + 1L]])))
    y2 <- p[[2L]]; x2 <- p[[1L]]
    if (i == n - 1L) { x1 <- R[[1L]]; y1 <- R[[2L]] } else { x1 <- pos(i + 1L)[[1L]]; y1 <- pos(i + 1L)[[2L]] }
    a <- atan((y1 - y2) / (x1 - x2))
    t <- -mom_l / (Abs(L[[2L]] - p[[2L]]) * cos(a) + Abs(L[[1L]] - p[[1L]]) * sin(a))
    tension <- .pd_set(tension, label, t)
    tfun[[length(tfun) + 1L]] <- list(t, x <= x1)
    mom_r <- mom_r + mag(i) * cos(ang(i)) * Abs(R[[2L]] - p[[2L]])
    mom_r <- mom_r + mag(i) * sin(ang(i)) * Abs(R[[1L]] - p[[1L]])
  }
  label <- Symbol(paste0(.cable_key(sp[[1L]]), "_", .cable_key(sp[[2L]])))
  p <- pos(2L)
  a <- -atan((p[[2L]] - L[[2L]]) / (p[[1L]] - L[[1L]]))
  t <- -mom_r / (Abs(R[[2L]] - p[[2L]]) * cos(a) + Abs(R[[1L]] - p[[1L]]) * sin(a))
  tension <- .pd_set(tension, label, t)
  .cm_set(self, "_tension", tension)
  .cm_set(self, "_length", len)
  tfun <- c(list(list(t, x <= p[[1L]])), tfun)
  .cm_set(self, "_tension_func", do.call(Piecewise, tfun))
  a <- S_Pi / 2L - a
  rl <- .cm_get(self, "_reaction_loads")
  rl <- .pd_set(rl, Symbol(paste0("R_", labels[[1L]], "_x")), -sin(a) * t)
  F_x <- F_x - sin(a) * t
  rl <- .pd_set(rl, Symbol(paste0("R_", labels[[1L]], "_y")), cos(a) * t)
  F_y <- F_y + cos(a) * t
  rl <- .pd_set(rl, Symbol(paste0("R_", labels[[2L]], "_x")), -F_x)
  rl <- .pd_set(rl, Symbol(paste0("R_", labels[[2L]], "_y")), -F_y)
  .cm_set(self, "_reaction_loads", rl)
}

.cable_solve_distributed <- function(self, lowest_x) {
  L <- .cm_get(self, "_left_support"); R <- .cm_get(self, "_right_support")
  lowest_x <- as_sym(lowest_x)
  .cm_set(self, "_lowest_x_global", lowest_x)
  a <- Symbol("a", positive = TRUE); c <- Symbol("c")
  sol <- linsolve(list((L[[1L]] - lowest_x)^2 * a + c - L[[2L]], (R[[1L]] - lowest_x)^2 * a + c - R[[2L]]), list(a, c))
  if (is.null(sol) || is.null(sol[["a"]]) || is_zero_num(sol[["a"]]))
    .physics_error("ValueError", "The lowest point is inconsistent with the supports")
  sa <- sol[["a"]]; sc <- sol[["c"]] %||% c
  A <- sa; C <- sc + sa * lowest_x^2; B <- -2L * sa * lowest_x
  .cm_set(self, "_lowest_y_global", sc)
  lowest_y <- sc
  X <- Symbol("X")
  Y <- A * (X + lowest_x)^2 + B * (X + lowest_x) + C - lowest_y
  force <- .cable_loads(self, "distributed")$values[[1L]]
  H <- (force * (R[[1L]] - lowest_x)^2) / (2L * (R[[2L]] - lowest_y))
  slope <- diff(Y, X)
  .cm_set(self, "_tension", .pd_set(.pd_empty(), "distributed", H / cos(atan(slope))))
  labels <- .cm_get(self, "_support_labels")
  t0 <- self$tension_at(L[[1L]])
  rl <- .cm_get(self, "_reaction_loads")
  rl <- .pd_set(rl, Symbol(paste0("R_", labels[[1L]], "_x")), t0 * cos(atan(subs(slope, X, L[[1L]] - lowest_x))))
  rl <- .pd_set(rl, Symbol(paste0("R_", labels[[1L]], "_y")), t0 * sin(atan(subs(slope, X, L[[1L]] - lowest_x))))
  rl <- .pd_set(rl, Symbol(paste0("R_", labels[[2L]], "_x")), t0 * cos(atan(subs(slope, X, R[[1L]] - lowest_x))))
  rl <- .pd_set(rl, Symbol(paste0("R_", labels[[2L]], "_y")), t0 * sin(atan(subs(slope, X, R[[1L]] - lowest_x))))
  .cm_set(self, "_reaction_loads", rl)
}

.cable_solve <- function(self, ...) {
  args <- list(...)
  if (length(.cm_get(self, "_loads_position")$keys)) {
    .cable_solve_point(self)
  } else if (length(.cable_loads(self, "distributed")$keys)) {
    if (!length(args)) .physics_error("ValueError", "Provide the lowest point of the cable")
    .cable_solve_distributed(self, args[[1L]])
  }
  invisible(NULL)
}

# a base R drawing of the cable, its supports and its loads
.cable_draw <- function(self) {
  L <- .cm_get(self, "_left_support"); R <- .cm_get(self, "_right_support")
  num <- function(e) as.numeric(N(as_sym(e)))
  lp <- .cm_get(self, "_loads_position")
  if (length(lp$keys)) {
    o <- order(vapply(lp$values, function(p) num(p[[1L]]), 0))
    xs <- c(num(L[[1L]]), vapply(lp$values[o], function(p) num(p[[1L]]), 0), num(R[[1L]]))
    ys <- c(num(L[[2L]]), vapply(lp$values[o], function(p) num(p[[2L]]), 0), num(R[[2L]]))
  } else if (length(.cable_loads(self, "distributed")$keys)) {
    if (!.pd_has(.cm_get(self, "_tension"), "distributed"))
      .physics_error("ValueError", "solve method not called and/or values provided for loads and supports not adequate")
    x0 <- num(.cm_get(self, "_lowest_x_global")); y0 <- num(.cm_get(self, "_lowest_y_global"))
    k <- (num(L[[2L]]) - y0) / (num(L[[1L]]) - x0)^2
    xs <- seq(num(L[[1L]]), num(R[[1L]]), length.out = 200L)
    ys <- k * (xs - x0)^2 + y0
  } else {
    .physics_error("ValueError", "solve method not called and/or values provided for loads and supports not adequate")
  }
  graphics::plot(xs, ys, type = "l", col = "brown", lwd = 2, asp = 1, axes = FALSE, xlab = "", ylab = "")
  graphics::points(c(xs[1L], xs[length(xs)]), c(ys[1L], ys[length(ys)]), pch = 15, cex = 2, col = "brown")
  if (length(lp$keys)) {
    pl <- .cable_loads(self, "point_load")
    for (k in seq_along(pl$keys)) {
      p <- .pd_get(lp, pl$keys[[k]]); a <- num(pl$values[[k]][[2L]]) * pi / 180
      len <- 0.1 * diff(range(xs))
      graphics::arrows(num(p[[1L]]), num(p[[2L]]), num(p[[1L]]) + len * cos(a), num(p[[2L]]) + len * sin(a), length = 0.1)
      graphics::text(num(p[[1L]]) + 1.6 * len * cos(a), num(p[[2L]]) + 1.6 * len * sin(a), paste0(.pyfmt(pl$values[[k]][[1L]]), "N"))
    }
  }
  invisible(NULL)
}

.cable_plot_tension <- function(self) {
  L <- .cm_get(self, "_left_support"); R <- .cm_get(self, "_right_support")
  if (length(.cm_get(self, "_loads_position")$keys))
    return(.cm_plot(.cm_get(self, "_tension_func"), Symbol("x"), R[[1L]], "Tension", "T", "blue", lower = L[[1L]]))
  .cm_plot(.pd_get(.cm_get(self, "_tension"), "distributed"), Symbol("X"), R[[1L]], "Tension", "T", "blue", lower = L[[1L]])
}

local({
  field <- function(f) function(self) .cm_get(self, f)
  .pyclass_def("Cable",
    props = list(supports = field("_supports"), left_support = field("_left_support"),
                 right_support = field("_right_support"), loads = field("_loads"),
                 loads_position = field("_loads_position"), length = field("_length"),
                 reaction_loads = field("_reaction_loads"), tension = field("_tension")),
    methods = list(tension_at = .cable_tension_at, apply_length = .cable_apply_length,
                   change_support = .cable_change_support, apply_load = .cable_apply_load,
                   remove_loads = .cable_remove_loads, solve = .cable_solve, draw = .cable_draw,
                   plot_tension = .cable_plot_tension))
})
