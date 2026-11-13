# Three hinged parabolic arches, after
# sympy/physics/continuum_mechanics/arch.py.

#' Arches
#'
#' `Arch(left_support, right_support, crown_x = , crown_y = )` is a three
#' hinged parabolic arch between two supports, each given as `list(x, y)`.
#' The crown abscissa is needed; without `crown_y` the parabola is fitted
#' through the supports. `$apply_load(order, label, start, mag, end, angle)`
#' adds a point load (`order = -1`, with the direction `angle` in degrees)
#' or a vertical load distributed from `start` to `end` (`order = 0`), and
#' `$remove_load()` removes one. `$change_support_type()` makes a support a
#' `"roller"` or a `"hinge"`, `$add_member(y)` adds a tie at height `y`,
#' and `$solve()` finds the reaction forces, after which
#' `$bending_moment_at()`, `$shear_force_at()` and `$axial_force_at()`
#' evaluate the internal forces. The properties are `$get_shape_eqn`,
#' `$get_loads`, `$supports`, `$left_support`, `$right_support` and
#' `$reaction_force`. `$draw()` draws the arch with base R graphics. After
#' `sympy.physics.continuum_mechanics.arch`.
#'
#' @param left_support,right_support The supports, as `list(x, y)`.
#' @param ... `crown_x` and `crown_y`, the position of the crown hinge.
#' @return An `Arch` object.
#' @examples
#' a <- Arch(list(0, 0), list(10, 0), crown_x = 5, crown_y = 5)
#' a$get_shape_eqn
#' a$apply_load(-1, "C", 3, 10, angle = 270)
#' a$solve()
#' a$reaction_force
#' @export
Arch <- function(left_support, right_support, ...) {
  kw <- list(...)
  self <- .pyenv_new("Arch")
  .cm_set(self, "_shape_eqn", NULL)
  .cm_set(self, "_left_support", .pytuple(list(as_sym(left_support[[1L]]), as_sym(left_support[[2L]]))))
  .cm_set(self, "_right_support", .pytuple(list(as_sym(right_support[[1L]]), as_sym(right_support[[2L]]))))
  .cm_set(self, "_crown_x", if (!is.null(kw$crown_x)) as_sym(kw$crown_x))
  .cm_set(self, "_crown_y", if (!is.null(kw$crown_y)) as_sym(kw$crown_y))
  .cm_set(self, "_shape_eqn", self$get_shape_eqn)
  for (f in c("_conc_loads", "_distributed_loads", "_loads_applied", "_moment_x", "_moment_y", "_load_x", "_load_y"))
    .cm_set(self, f, .pd_empty())
  .cm_set(self, "_supports", .pydict(list("left", "right"), list("hinge", "hinge")))
  .cm_set(self, "_member", NULL)
  .cm_set(self, "_member_force", NULL)
  .cm_set(self, "_reaction_force", .pydict(lapply(c("R_A_x", "R_A_y", "R_B_x", "R_B_y"), Symbol), rep(list(0L), 4L)))
  .cm_set(self, "_points_disc_x", list())
  .cm_set(self, "_points_disc_y", list())
  for (f in c("_moment_x_func", "_moment_y_func", "_load_x_func", "_load_y_func"))
    .cm_set(self, f, Piecewise(list(0L, S_true)))
  for (f in c("_bending_moment", "_shear_force", "_axial_force")) .cm_set(self, f, NULL)
  self
}

# Python's truth value of an optional number
.arch_truthy <- function(v) !is.null(v) && !is_zero_num(v)

.arch_shape <- function(self) {
  shape <- .cm_get(self, "_shape_eqn")
  if (.arch_truthy(shape)) return(shape)
  x <- Symbol("x"); y <- Symbol("y"); c <- Symbol("c")
  a <- Symbol("a", positive = FALSE)
  L <- .cm_get(self, "_left_support"); R <- .cm_get(self, "_right_support")
  x0 <- .cm_get(self, "_crown_x"); y0 <- .cm_get(self, "_crown_y")
  if (.arch_truthy(x0) && .arch_truthy(y0)) {
    sol <- linsolve(list(a * (L[[1L]] - x0)^2 + y0 - L[[2L]]), list(a))
    if (is.null(sol) || is.null(sol[["a"]])) .physics_error("IndexError", "list index out of range")
    eqn <- sol[["a"]] * (x - x0)^2 + y0
    if (!.cm_same(subs(eqn, x, R[[1L]]), R[[2L]]))
      .physics_error("ValueError", "provided coordinates of crown and supports are not consistent with parabolic arch")
  } else if (.arch_truthy(x0)) {
    sol <- linsolve(list(a * (L[[1L]] - x0)^2 + c - L[[2L]], a * (R[[1L]] - x0)^2 + c - R[[2L]]), list(a, c))
    if (is.null(sol) || length(sol) < 2L || is_zero_num(sol[["a"]]))
      .physics_error("ValueError", "parabolic arch cannot be constructed with the provided coordinates, try providing crown_y")
    eqn <- sol[["a"]] * (x - x0)^2 + sol[["c"]]
    .cm_set(self, "_crown_y", sol[["c"]])
  } else {
    .physics_error("KeyError", "'please provide crown_x to construct arch'")
  }
  eqn
}

# a set of numbers, as a list without duplicates
.arch_set_add <- function(l, v) if (is.na(.pl_index(l, v))) c(l, list(v)) else l

.arch_dict_add <- function(self, f, k, v) {
  d <- .cm_get(self, f)
  .cm_set(self, f, .pd_set(d, k, if (.pd_has(d, k)) .pd_get(d, k) + v else v))
}

.arch_apply_load <- function(self, order, label, start, mag, end = NULL, angle = NULL) {
  y <- Symbol("y"); x <- Symbol("x"); x0 <- Symbol("x0")
  order <- as_sym(order); mag <- as_sym(mag)
  if (!is.null(angle)) angle <- as_sym(angle)
  if (.pd_has(.cm_get(self, "_loads_applied"), label)) .physics_error("ValueError", "load with the given label already exists")
  if (label %in% c("A", "B")) .physics_error("ValueError", "cannot use the given label, reserved for supports")
  start <- as_sym(start)
  if (is_zero_num(order)) {
    if (is.null(end) || isTRUE(as_sym(end) < start)) .physics_error("KeyError", "'provide end greater than start'")
    end <- as_sym(end)
    .cm_set(self, "_distributed_loads", .pd_set(.cm_get(self, "_distributed_loads"), label,
                                                .pydict(list("start", "end", "f_y"), list(start, end, mag))))
    .cm_set(self, "_points_disc_y", .arch_set_add(.cm_get(self, "_points_disc_y"), start))
    m <- Min(x, end)
    .arch_dict_add(self, "_moment_y", start, -mag * (m - start) * (x0 - (start + m) / 2L))
    .arch_dict_add(self, "_load_y", start, mag * (Min(end, x) - start))
    .cm_set(self, "_loads_applied", .pd_set(.cm_get(self, "_loads_applied"), label, "distributed"))
  }
  if (isTRUE(order == -1L) || identical(order$key, as_sym(-1L)$key)) {
    if (is.null(angle)) .physics_error("TypeError", "please provide direction of force")
    height <- subs(.cm_get(self, "_shape_eqn"), x, start)
    f_x <- mag * cos(rad(angle)); f_y <- mag * sin(rad(angle))
    load <- .pydict(list("x", "y", "f_x", "f_y", "mag", "angle"), list(start, height, f_x, f_y, mag, angle))
    .cm_set(self, "_conc_loads", .pd_set(.cm_get(self, "_conc_loads"), label, load))
    .cm_set(self, "_points_disc_x", .arch_set_add(.cm_get(self, "_points_disc_x"), start))
    .cm_set(self, "_points_disc_y", .arch_set_add(.cm_get(self, "_points_disc_y"), start))
    .arch_dict_add(self, "_moment_x", start, f_x * (y - height))
    .arch_dict_add(self, "_load_x", start, f_x)
    .arch_dict_add(self, "_moment_y", start, -f_y * (x0 - start))
    .arch_dict_add(self, "_load_y", start, f_y)
    .cm_set(self, "_loads_applied", .pd_set(.cm_get(self, "_loads_applied"), label, "concentrated"))
  }
  invisible(NULL)
}

.arch_remove_load <- function(self, label) {
  y <- Symbol("y"); x <- Symbol("x"); x0 <- Symbol("x0")
  dl <- .cm_get(self, "_distributed_loads"); cl <- .cm_get(self, "_conc_loads")
  set_remove <- function(f, v) {
    l <- .cm_get(self, f)
    if (is.na(.pl_index(l, v))) .physics_error("KeyError", .pyfmt(v))
    .cm_set(self, f, .pl_remove(l, v))
  }
  if (.pd_has(dl, label)) {
    .cm_set(self, "_loads_applied", .pd_pop(.cm_get(self, "_loads_applied"), label))
    ld <- .pd_get(dl, label)
    start <- .pd_get(ld, "start"); end <- .pd_get(ld, "end"); mag <- .pd_get(ld, "f_y")
    set_remove("_points_disc_y", start)
    .arch_dict_add(self, "_load_y", start, -mag * (Min(x, end) - start))
    .arch_dict_add(self, "_moment_y", start, mag * (Min(x, end) - start) * (x0 - (start + Min(x, end)) / 2L))
    .cm_set(self, "_distributed_loads", .pd_pop(dl, label))
    cat(sprintf("removed load %s: %s\n", label, .arch_fmt(ld)))
  } else if (.pd_has(cl, label)) {
    .cm_set(self, "_loads_applied", .pd_pop(.cm_get(self, "_loads_applied"), label))
    ld <- .pd_get(cl, label)
    start <- .pd_get(ld, "x")
    set_remove("_points_disc_x", start)
    set_remove("_points_disc_y", start)
    .arch_dict_add(self, "_moment_y", start, .pd_get(ld, "f_y") * (x0 - start))
    .arch_dict_add(self, "_moment_x", start, -.pd_get(ld, "f_x") * (y - .pd_get(ld, "y")))
    .arch_dict_add(self, "_load_x", start, -.pd_get(ld, "f_x"))
    .arch_dict_add(self, "_load_y", start, -.pd_get(ld, "f_y"))
    .cm_set(self, "_conc_loads", .pd_pop(cl, label))
    cat(sprintf("removed load %s: %s\n", label, .arch_fmt(ld)))
  } else {
    .physics_error("ValueError", "label not found")
  }
  invisible(NULL)
}

# repr() of a dictionary with string keys
.arch_fmt <- function(d) paste0("{", paste(vapply(seq_along(d$keys), function(i)
  paste0("'", d$keys[[i]], "': ", .pyfmt(d$values[[i]])), ""), collapse = ", "), "}")

.arch_add_member <- function(self, y) {
  L <- .cm_get(self, "_left_support"); R <- .cm_get(self, "_right_support")
  cy <- .cm_get(self, "_crown_y"); cx <- .cm_get(self, "_crown_x")
  ys <- as_sym(y)
  lo <- if (isTRUE(R[[2L]] < L[[2L]])) R[[2L]] else L[[2L]]
  if (isTRUE(ys > cy) || isTRUE(ys < lo))
    .physics_error("ValueError", sprintf("position of support must be between y=%s and y=%s", .pyfmt(lo), .pyfmt(cy)))
  x <- Symbol("x")
  a <- subs(diff(.cm_get(self, "_shape_eqn"), x), x, cx + 1L) / 2L
  xd <- sqrt((ys - cy) / a)
  .cm_set(self, "_member", .pytuple(list(cx + xd, cx - xd, y)))
  invisible(NULL)
}

.arch_at <- function(field, var) function(self, pos = NULL, ...) {
  kw <- list(...)
  e <- .cm_get(self, field)
  if (is.null(pos)) return(e)
  if (!is.null(kw$dir)) return(limit(e, Symbol(var), pos, dir = kw$dir))
  subs(e, Symbol(var), pos)
}

.arch_solve <- function(self) {
  y <- Symbol("y"); x <- Symbol("x"); x0 <- Symbol("x0")
  L <- .cm_get(self, "_left_support"); Rs <- .cm_get(self, "_right_support")
  cx <- .cm_get(self, "_crown_x"); cy <- .cm_get(self, "_crown_y")
  srt <- .cm_sorted
  mx <- Piecewise(list(0L, S_true)); my <- mx; lx <- mx; ly <- mx
  ax <- S_Zero; ay <- S_Zero; lxa <- S_Zero; lya <- S_Zero
  for (p in srt(.cm_get(self, "_points_disc_x"))) {
    cond <- x >= p
    lxa <- lxa + .pd_get(.cm_get(self, "_load_x"), p)
    ax <- ax + .pd_get(.cm_get(self, "_moment_x"), p)
    lx <- Piecewise(list(lxa, cond), list(lx, S_true))
    mx <- Piecewise(list(ax, cond), list(mx, S_true))
  }
  for (p in srt(.cm_get(self, "_points_disc_y"))) {
    cond <- x >= p
    ay <- ay + .pd_get(.cm_get(self, "_moment_y"), p)
    lya <- lya + .pd_get(.cm_get(self, "_load_y"), p)
    ly <- Piecewise(list(lya, cond), list(ly, S_true))
    my <- Piecewise(list(ay, cond), list(my, S_true))
  }
  for (f in c("_moment_x_func", "_moment_y_func", "_load_x_func", "_load_y_func"))
    .cm_set(self, f, switch(f, "_moment_x_func" = mx, "_moment_y_func" = my, "_load_x_func" = lx, ly))
  s2 <- function(e, a, va, b, vb) subs(subs(e, a, va), b, vb)
  moment_A <- s2(my, x, Rs[[1L]], x0, L[[1L]]) + s2(mx, x, Rs[[1L]], y, L[[2L]])
  mh_left <- s2(my, x, cx, x0, cx) + s2(mx, x, cx, y, cy)
  mh_right <- s2(my, x, Rs[[1L]], x0, cx) - s2(my, x, cx, x0, cx) + s2(mx, x, Rs[[1L]], y, cy) - s2(mx, x, cx, y, cy)
  net_x <- subs(lx, x, Rs[[1L]]); net_y <- subs(ly, x, Rs[[1L]])
  sup <- .cm_get(self, "_supports")
  left <- .pd_get(sup, "left"); right <- .pd_get(sup, "right")
  member <- .cm_get(self, "_member")
  if ((left == "roller" || right == "roller") && is.null(member)) {
    cat("member must be added if any of the supports is roller\n")
    return(invisible(NULL))
  }
  RAx <- Symbol("R_A_x"); RAy <- Symbol("R_A_y"); RBx <- Symbol("R_B_x"); RBy <- Symbol("R_B_y"); Tt <- Symbol("T")
  ge <- function(a, b) isTRUE(as_sym(a) >= as_sym(b))
  mxs <- function(a, b) if (isTRUE(a < b)) b else a
  eqs <- NULL
  if (!is.null(member)) m2 <- as_sym(member[[3L]])
  span <- Rs[[1L]] - L[[1L]]; rise <- Rs[[2L]] - L[[2L]]
  if (left == "roller" && right == "roller") {
    if (ge(m2, mxs(L[[2L]], Rs[[2L]]))) {
      if (!is_zero_num(net_x)) .physics_error("ValueError", "net force in x direction not possible under the specified conditions")
      eqs <- list(RAx, RBx, RAy + RBy + net_y, RBy * span - RBx * rise + moment_A,
                  mh_right + RBy * (Rs[[1L]] - cx) + Tt * (m2 - cy))
    } else if (ge(m2, L[[2L]])) {
      eqs <- list(RAx, RBx, RAy + RBy + net_y, RBy * span - Tt * (m2 - L[[2L]]) + moment_A, Tt + net_x)
    } else if (ge(m2, Rs[[2L]])) {
      eqs <- list(RAx, RBx, RAy + RBy + net_y, RBy * span + Tt * (m2 - L[[2L]]) + moment_A, Tt - net_x)
    }
  } else if (left == "roller") {
    if (ge(m2, mxs(L[[2L]], Rs[[2L]]))) {
      eqs <- list(RAx, RBx + net_x, RAy + RBy + net_y, RBy * span - RBx * rise + moment_A,
                  mh_left + RAy * (L[[1L]] - cx) - Tt * (m2 - cy))
    } else if (ge(m2, L[[2L]])) {
      eqs <- list(RAx, RBx + Tt + net_x, RAy + RBy + net_y, RBy * span - RBx * rise - Tt * (m2 - L[[1L]]) + moment_A,
                  mh_left + RAy * (L[[1L]] - cx) - Tt * (m2 - cy))
    } else if (ge(m2, Rs[[1L]])) {
      eqs <- list(RAx, RBx - Tt + net_x, RAy + RBy + net_y, mh_left + RAy * (L[[1L]] - cx),
                  moment_A + RBy * span - RBx * rise + Tt * (m2 - L[[2L]]))
    }
  } else if (right == "roller") {
    if (ge(m2, mxs(L[[2L]], Rs[[2L]]))) {
      eqs <- list(RBx, RAx + net_x, RAy + RBy + net_y, mh_right + RBy * (Rs[[1L]] - cx) + Tt * (m2 - cy),
                  moment_A + RBy * span)
    } else if (ge(m2, L[[2L]])) {
      eqs <- list(RBx, RAx + Tt + net_x, RAy + RBy + net_y, mh_right + RBy * (Rs[[1L]] - cx),
                  moment_A - Tt * (m2 - L[[2L]]) + RBy * span)
    } else if (ge(m2, Rs[[2L]])) {
      # SymPy builds this equation with a single argument
      .physics_error("TypeError", "Equality.__new__() missing 1 required positional argument: 'rhs'")
    }
  } else {
    eqs <- list(RAx + RBx + net_x, RAy + RBy + net_y, RBy * span - RBx * rise + moment_A,
                mh_right + RBy * (Rs[[1L]] - cx) - RBx * (Rs[[2L]] - cy))
  }
  if (is.null(eqs))
    .physics_error("UnboundLocalError", "cannot access local variable 'solution' where it is not associated with a value")
  unknowns <- if (length(eqs) == 5L) list(RAx, RAy, RBx, RBy, Tt) else list(RAx, RAy, RBx, RBy)
  sol <- .cm_linsolve(eqs, unknowns)
  names(sol) <- vapply(unknowns, function(u) u$name, "")
  rf <- .cm_get(self, "_reaction_force")
  for (k in rf$keys) rf <- .pd_set(rf, k, sol[[k$name]])
  .cm_set(self, "_reaction_force", rf)
  shape <- .cm_get(self, "_shape_eqn")
  .cm_set(self, "_bending_moment", -(subs(mx, x, x0) + subs(my, x, x0) - sol[["R_A_y"]] * (x0 - L[[1L]]) +
                                       sol[["R_A_x"]] * (subs(shape, x, x0) - L[[2L]])))
  angle <- atan(diff(shape, x))
  fx <- lx + sol[["R_A_x"]]; fy <- ly + sol[["R_A_y"]]
  .cm_set(self, "_axial_force", fx * cos(angle) + fy * sin(angle))
  .cm_set(self, "_shear_force", -fx * sin(angle) + fy * cos(angle))
  invisible(NULL)
}

# a base R drawing of the arch, its supports, its tie and its loads
.arch_draw <- function(self) {
  L <- .cm_get(self, "_left_support"); R <- .cm_get(self, "_right_support")
  num <- function(e) as.numeric(N(as_sym(e)))
  x <- Symbol("x")
  shape <- .cm_get(self, "_shape_eqn")
  xs <- seq(num(L[[1L]]), num(R[[1L]]), length.out = 200L)
  ys <- vapply(xs, function(v) num(subs(shape, x, v)), 0)
  graphics::plot(xs, ys, type = "l", col = "brown", lwd = 3, asp = 1, axes = FALSE, xlab = "", ylab = "")
  graphics::points(c(xs[1L], xs[200L]), c(ys[1L], ys[200L]), pch = 17, cex = 2, col = "black")
  graphics::points(num(.cm_get(self, "_crown_x")), num(.cm_get(self, "_crown_y")), pch = 1)
  m <- .cm_get(self, "_member")
  if (!is.null(m)) graphics::segments(num(m[[2L]]), num(m[[3L]]), num(m[[1L]]), num(m[[3L]]), col = "grey40", lwd = 2)
  cl <- .cm_get(self, "_conc_loads")
  len <- 0.1 * diff(range(xs))
  for (ld in cl$values) {
    a <- num(.pd_get(ld, "angle")) * pi / 180
    px <- num(.pd_get(ld, "x")); py <- num(.pd_get(ld, "y"))
    graphics::arrows(px - len * cos(a), py - len * sin(a), px, py, length = 0.1, col = "orange")
  }
  dl <- .cm_get(self, "_distributed_loads")
  for (ld in dl$values) {
    st <- num(.pd_get(ld, "start")); en <- num(.pd_get(ld, "end"))
    graphics::rect(st, max(ys) + 0.05 * len, en, max(ys) + 0.5 * len, col = "orange", border = NA)
  }
  invisible(NULL)
}

local({
  field <- function(f) function(self) .cm_get(self, f)
  .pyclass_def("Arch",
    props = list(get_shape_eqn = .arch_shape,
                 get_loads = function(self) .pydict(list("concentrated", "distributed"),
                                                    list(.cm_get(self, "_conc_loads"), .cm_get(self, "_distributed_loads"))),
                 supports = field("_supports"), left_support = field("_left_support"),
                 right_support = field("_right_support"), reaction_force = field("_reaction_force")),
    methods = list(
      apply_load = .arch_apply_load, remove_load = .arch_remove_load,
      change_support_position = function(self, left_support = NULL, right_support = NULL) {
        if (!is.null(left_support)) .cm_set(self, "_left_support", .pytuple(list(left_support[[1L]], left_support[[2L]])))
        if (!is.null(right_support)) .cm_set(self, "_right_support", .pytuple(list(right_support[[1L]], right_support[[2L]])))
        .cm_set(self, "_shape_eqn", NULL)
        .cm_set(self, "_shape_eqn", self$get_shape_eqn)
        invisible(NULL)
      },
      change_crown_position = function(self, crown_x = NULL, crown_y = NULL) {
        .cm_set(self, "_crown_x", if (!is.null(crown_x)) as_sym(crown_x))
        .cm_set(self, "_crown_y", if (!is.null(crown_y)) as_sym(crown_y))
        .cm_set(self, "_shape_eqn", NULL)
        .cm_set(self, "_shape_eqn", self$get_shape_eqn)
        invisible(NULL)
      },
      change_support_type = function(self, left_support = NULL, right_support = NULL) {
        sup <- .cm_get(self, "_supports")
        for (side in c("left", "right")) {
          v <- if (side == "left") left_support else right_support
          if (is.null(v) || !nzchar(v)) next
          if (!v %in% c("roller", "hinge")) .physics_error("ValueError", "supports must only be roller or hinge")
          sup <- .pd_set(sup, side, v)
        }
        .cm_set(self, "_supports", sup)
        invisible(NULL)
      },
      add_member = .arch_add_member,
      shear_force_at = .arch_at("_shear_force", "x"),
      bending_moment_at = .arch_at("_bending_moment", "x0"),
      axial_force_at = .arch_at("_axial_force", "x"),
      solve = .arch_solve, draw = .arch_draw))
})
