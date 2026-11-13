# Two dimensional trusses, after
# sympy/physics/continuum_mechanics/truss.py.

#' Trusses
#'
#' `Truss()` is a plane truss of nodes joined by members. `$add_node()`,
#' `$remove_node()`, `$add_member()`, `$remove_member()`,
#' `$change_node_label()` and `$change_member_label()` build it; each takes
#' any number of tuples, such as `list(label, x, y)` for a node and
#' `list(label, start, end)` for a member. `$apply_load(list(node,
#' magnitude, direction))` adds a load with the direction in degrees,
#' `$apply_support(list(node, type))` a `"pinned"` or `"roller"` support,
#' and `$solve()` finds the reaction loads and the internal forces. The
#' properties are `$nodes`, `$node_labels`, `$node_positions`, `$members`,
#' `$member_lengths`, `$supports`, `$loads`, `$reaction_loads` and
#' `$internal_forces`. `$draw()` draws the truss with base R graphics. After
#' `sympy.physics.continuum_mechanics.truss`.
#'
#' @return A `Truss` object.
#' @examples
#' t <- Truss()
#' t$add_node(list("A", 0, 0), list("B", 2, 2), list("C", 3, 0))
#' t$add_member(list("AB", "A", "B"), list("BC", "B", "C"), list("AC", "A", "C"))
#' t$apply_load(list("B", Symbol("P"), -90))
#' t$apply_support(list("A", "pinned"), list("C", "roller"))
#' t$solve()
#' t$reaction_loads
#' @export
Truss <- function() {
  self <- .pyenv_new("Truss")
  for (f in c("_nodes", "_node_labels", "_node_positions", "_node_position_x", "_node_position_y"))
    .cm_set(self, f, list())
  for (f in c("_members", "_loads", "_supports", "_nodes_occupied", "_member_lengths", "_reaction_loads",
              "_internal_forces", "_node_coordinates"))
    .cm_set(self, f, .pd_empty())
  self
}

# the items of a method taking tuples: each must be subscriptable
.truss_items <- function(...) {
  lapply(list(...), function(i) {
    if (is_sym(i) || !is.list(i))
      .physics_error("TypeError", sprintf("'%s' object is not subscriptable", .py_typename(i)))
    i
  })
}

.truss_add_node <- function(self, ...) {
  for (i in .truss_items(...)) {
    label <- i[[1L]]; x <- as_sym(i[[2L]]); y <- as_sym(i[[3L]])
    nc <- .cm_get(self, "_node_coordinates")
    if (.pd_has(nc, label)) .physics_error("ValueError", "Node needs to have a unique label")
    if (any(vapply(nc$values, function(v) .cm_tuple_eq(v, list(x, y)), TRUE)))
      .physics_error("ValueError", "A node already exists at the given position")
    .cm_push(self, "_nodes", .pytuple(list(label, x, y)))
    .cm_push(self, "_node_labels", label)
    .cm_push(self, "_node_positions", .pytuple(list(x, y)))
    .cm_push(self, "_node_position_x", x)
    .cm_push(self, "_node_position_y", y)
    .cm_set(self, "_node_coordinates", .pd_set(nc, label, list(x, y)))
  }
  invisible(NULL)
}

.truss_remove_node <- function(self, ...) {
  for (label in list(...)) {
    if (!.pd_has(.cm_get(self, "_node_coordinates"), label)) .physics_error("ValueError", "No such node exists in the truss")
    labels <- .cm_get(self, "_node_labels")
    k <- .pl_index(labels, label)
    x <- .cm_get(self, "_node_position_x")[[k]]; y <- .cm_get(self, "_node_position_y")[[k]]
    mem <- .cm_get(self, "_members")
    for (m in mem$values)
      if (.cm_same_key(label, m[[1L]]) || .cm_same_key(label, m[[2L]]))
        .physics_error("ValueError", "The given node already has member attached to it")
    .cm_set(self, "_nodes", .pl_remove(.cm_get(self, "_nodes"), list(label, x, y)))
    .cm_set(self, "_node_labels", .pl_remove(labels, label))
    .cm_set(self, "_node_positions", .pl_remove(.cm_get(self, "_node_positions"), list(x, y)))
    .cm_set(self, "_node_position_x", .pl_remove(.cm_get(self, "_node_position_x"), x))
    .cm_set(self, "_node_position_y", .pl_remove(.cm_get(self, "_node_position_y"), y))
    for (f in c("_loads", "_supports")) if (.pd_has(.cm_get(self, f), label)) .cm_set(self, f, .pd_pop(.cm_get(self, f), label))
    .cm_set(self, "_node_coordinates", .pd_pop(.cm_get(self, "_node_coordinates"), label))
  }
  invisible(NULL)
}

.cm_same_key <- function(a, b) identical(.pyobj_k(a), .pyobj_k(b))

.truss_add_member <- function(self, ...) {
  for (i in .truss_items(...)) {
    label <- i[[1L]]; start <- i[[2L]]; end <- i[[3L]]
    nc <- .cm_get(self, "_node_coordinates")
    occ <- .cm_get(self, "_nodes_occupied")
    if (!.pd_has(nc, start) || !.pd_has(nc, end) || .cm_same_key(start, end))
      .physics_error("ValueError", "The start and end points of the member must be unique nodes")
    if (.pd_has(.cm_get(self, "_members"), label))
      .physics_error("ValueError", "A member with the same label already exists for the truss")
    if (.pd_has(occ, .pytuple(list(start, end))))
      .physics_error("ValueError", "A member already exists between the two nodes")
    .cm_set(self, "_members", .pd_set(.cm_get(self, "_members"), label, list(start, end)))
    ps <- .pd_get(nc, start); pe <- .pd_get(nc, end)
    .cm_set(self, "_member_lengths", .pd_set(.cm_get(self, "_member_lengths"), label,
                                             sqrt((pe[[1L]] - ps[[1L]])^2 + (pe[[2L]] - ps[[2L]])^2)))
    occ <- .pd_set(occ, .pytuple(list(start, end)), TRUE)
    .cm_set(self, "_nodes_occupied", .pd_set(occ, .pytuple(list(end, start)), TRUE))
    .cm_set(self, "_internal_forces", .pd_set(.cm_get(self, "_internal_forces"), label, 0L))
  }
  invisible(NULL)
}

.truss_remove_member <- function(self, ...) {
  for (label in list(...)) {
    mem <- .cm_get(self, "_members")
    if (!.pd_has(mem, label)) .physics_error("ValueError", "No such member exists in the Truss")
    m <- .pd_get(mem, label)
    occ <- .pd_pop(.cm_get(self, "_nodes_occupied"), .pytuple(list(m[[1L]], m[[2L]])))
    .cm_set(self, "_nodes_occupied", .pd_pop(occ, .pytuple(list(m[[2L]], m[[1L]]))))
    for (f in c("_members", "_member_lengths", "_internal_forces")) .cm_set(self, f, .pd_pop(.cm_get(self, f), label))
  }
  invisible(NULL)
}

# load[0] -= R for the first load of the node in the given direction
.truss_shift_load <- function(self, new_label, label, direction, R) {
  loads <- .cm_get(self, "_loads")
  l <- .pd_get(loads, new_label)
  for (j in seq_along(l)) {
    if (.cm_same(l[[j]][[2L]], direction)) {
      l[[j]][[1L]] <- l[[j]][[1L]] - R
      if (is_zero_num(l[[j]][[1L]])) l <- l[-j]
      break
    }
  }
  .cm_set(self, "_loads", .pd_set(loads, new_label, l))
}

.truss_change_node_label <- function(self, ...) {
  for (i in .truss_items(...)) {
    label <- i[[1L]]; new_label <- i[[2L]]
    nc <- .cm_get(self, "_node_coordinates")
    if (!.pd_has(nc, label)) .physics_error("ValueError", "No such node exists for the Truss")
    if (.pd_has(nc, new_label)) .physics_error("ValueError", "A node with the given label already exists")
    nodes <- .cm_get(self, "_nodes")
    for (k in seq_along(nodes)) {
      node <- nodes[[k]]
      if (!.cm_same_key(node[[1L]], label)) next
      nodes[[k]] <- .pytuple(list(new_label, node[[2L]], node[[3L]]))
      .cm_set(self, "_nodes", nodes)
      labels <- .cm_get(self, "_node_labels")
      labels[[.pl_index(labels, label)]] <- new_label
      .cm_set(self, "_node_labels", labels)
      nc <- .pd_set(nc, new_label, .pd_get(nc, label))
      .cm_set(self, "_node_coordinates", .pd_pop(nc, label))
      sup <- .cm_get(self, "_supports")
      if (.pd_has(sup, label)) {
        sup <- .pd_pop(.pd_set(sup, new_label, .pd_get(sup, label)), label)
        .cm_set(self, "_supports", sup)
      }
      Rn <- function(l, d) Symbol(paste0("R_", .pyfmt(l), "_", d))
      if (.pd_has(sup, new_label)) {
        type <- .pd_get(sup, new_label)
        loads <- .cm_get(self, "_loads")
        if (type == "pinned") {
          rl <- .cm_get(self, "_reaction_loads")
          kx <- paste0("R_", .pyfmt(label), "_x"); ky <- paste0("R_", .pyfmt(label), "_y")
          if (.pd_has(rl, kx) && .pd_has(rl, ky)) {
            rl <- .pd_set(rl, paste0("R_", .pyfmt(new_label), "_x"), .pd_get(rl, kx))
            rl <- .pd_set(rl, paste0("R_", .pyfmt(new_label), "_y"), .pd_get(rl, ky))
            .cm_set(self, "_reaction_loads", .pd_pop(.pd_pop(rl, kx), ky))
          }
          .cm_set(self, "_loads", .pd_set(loads, new_label, .pd_get(loads, label)))
          .truss_shift_load(self, new_label, label, as_sym(90L), Rn(label, "y"))
          .truss_shift_load(self, new_label, label, S_Zero, Rn(label, "x"))
          # SymPy calls apply_load(new_label, R_x, 0) without a tuple here
          self$apply_load(new_label, Rn(new_label, "x"), 0L)
        } else if (type == "roller") {
          .cm_set(self, "_loads", .pd_set(loads, new_label, .pd_get(loads, label)))
          .truss_shift_load(self, new_label, label, as_sym(90L), Rn(label, "y"))
          self$apply_load(new_label, Rn(new_label, "y"), 90L)
        }
      } else {
        loads <- .cm_get(self, "_loads")
        if (.pd_has(loads, label)) .cm_set(self, "_loads", .pd_pop(.pd_set(loads, new_label, .pd_get(loads, label)), label))
      }
      mem <- .cm_get(self, "_members")
      occ <- .cm_get(self, "_nodes_occupied")
      for (j in seq_along(mem$values)) {
        m <- mem$values[[j]]
        if (.cm_same_key(m[[1L]], label)) {
          mem$values[[j]][[1L]] <- new_label
          occ <- .pd_set(.pd_set(occ, .pytuple(list(new_label, m[[2L]])), TRUE), .pytuple(list(m[[2L]], new_label)), TRUE)
          occ <- .pd_pop(.pd_pop(occ, .pytuple(list(label, m[[2L]]))), .pytuple(list(m[[2L]], label)))
        } else if (.cm_same_key(m[[2L]], label)) {
          mem$values[[j]][[2L]] <- new_label
          occ <- .pd_set(.pd_set(occ, .pytuple(list(m[[1L]], new_label)), TRUE), .pytuple(list(new_label, m[[1L]])), TRUE)
          occ <- .pd_pop(.pd_pop(occ, .pytuple(list(m[[1L]], label))), .pytuple(list(label, m[[1L]])))
        }
      }
      .cm_set(self, "_members", mem)
      .cm_set(self, "_nodes_occupied", occ)
    }
  }
  invisible(NULL)
}

.truss_change_member_label <- function(self, ...) {
  for (i in .truss_items(...)) {
    label <- i[[1L]]; new_label <- i[[2L]]
    if (!.pd_has(.cm_get(self, "_members"), label)) .physics_error("ValueError", "No such member exists for the Truss")
    for (f in c("_members", "_member_lengths", "_internal_forces")) {
      d <- .cm_get(self, f)
      .cm_set(self, f, .pd_pop(.pd_set(d, new_label, .pd_get(d, label)), label))
    }
  }
  invisible(NULL)
}

.truss_apply_load <- function(self, ...) {
  for (i in .truss_items(...)) {
    location <- i[[1L]]; magnitude <- as_sym(i[[2L]]); direction <- as_sym(i[[3L]])
    if (!.pd_has(.cm_get(self, "_node_coordinates"), location))
      .physics_error("ValueError", "Load must be applied at a known node")
    loads <- .cm_get(self, "_loads")
    l <- if (.pd_has(loads, location)) .pd_get(loads, location) else list()
    .cm_set(self, "_loads", .pd_set(loads, location, c(l, list(list(magnitude, direction)))))
  }
  invisible(NULL)
}

.truss_remove_load <- function(self, ...) {
  for (i in .truss_items(...)) {
    location <- i[[1L]]; magnitude <- as_sym(i[[2L]]); direction <- as_sym(i[[3L]])
    if (!.pd_has(.cm_get(self, "_node_coordinates"), location))
      .physics_error("ValueError", "Load must be removed from a known node")
    loads <- .cm_get(self, "_loads")
    l <- .pd_get(loads, location)
    if (is.na(.pl_index(l, list(magnitude, direction))))
      .physics_error("ValueError", "No load of this magnitude and direction has been applied at this node")
    l <- .pl_remove(l, list(magnitude, direction))
    loads <- .pd_set(loads, location, l)
    if (!length(l)) loads <- .pd_pop(loads, location)
    .cm_set(self, "_loads", loads)
  }
  invisible(NULL)
}

.truss_apply_support <- function(self, ...) {
  for (i in .truss_items(...)) {
    location <- i[[1L]]; type <- i[[2L]]
    if (!.pd_has(.cm_get(self, "_node_coordinates"), location))
      .physics_error("ValueError", "Support must be added on a known node")
    sup <- .cm_get(self, "_supports")
    R <- function(d) Symbol(paste0("R_", .pyfmt(location), "_", d))
    if (!.pd_has(sup, location)) {
      if (identical(type, "pinned")) {
        self$apply_load(list(location, R("x"), 0L))
        self$apply_load(list(location, R("y"), 90L))
      } else if (identical(type, "roller")) {
        self$apply_load(list(location, R("y"), 90L))
      }
    } else if (identical(.pd_get(sup, location), "pinned")) {
      if (identical(type, "roller")) self$remove_load(list(location, R("x"), 0L))
    } else if (identical(.pd_get(sup, location), "roller")) {
      if (identical(type, "pinned")) self$apply_load(list(location, R("x"), 0L))
    }
    .cm_set(self, "_supports", .pd_set(.cm_get(self, "_supports"), location, type))
  }
  invisible(NULL)
}

.truss_remove_support <- function(self, ...) {
  for (location in list(...)) {
    if (!.pd_has(.cm_get(self, "_node_coordinates"), location)) .physics_error("ValueError", "No such node exists in the Truss")
    sup <- .cm_get(self, "_supports")
    if (!.pd_has(sup, location)) .physics_error("ValueError", "No support has been added to the given node")
    R <- function(d) Symbol(paste0("R_", .pyfmt(location), "_", d))
    type <- .pd_get(sup, location)
    if (identical(type, "pinned")) {
      self$remove_load(list(location, R("x"), 0L))
      self$remove_load(list(location, R("y"), 90L))
    } else if (identical(type, "roller")) {
      self$remove_load(list(location, R("y"), 90L))
    }
    .cm_set(self, "_supports", .pd_pop(.cm_get(self, "_supports"), location))
  }
  invisible(NULL)
}

.truss_solve <- function(self) {
  nodes <- .cm_get(self, "_nodes"); sup <- .cm_get(self, "_supports"); loads <- .cm_get(self, "_loads")
  nc <- .cm_get(self, "_node_coordinates"); mem <- .cm_get(self, "_members")
  labels <- .cm_get(self, "_node_labels")
  stype <- function(l) if (.pd_has(sup, l)) .pd_get(sup, l) else ""
  count <- 0L
  for (node in nodes) count <- count + switch(stype(node[[1L]]), pinned = 2L, roller = 1L, 0L)
  n2 <- 2L * length(nodes)
  if (n2 != length(mem$keys) + count) .physics_error("ValueError", "The given truss cannot be solved")
  coef <- rep(list(S_Zero), n2 * n2)
  at <- function(r, c) r * n2 + c + 1L
  add <- function(r, c, v) coef[[at(r, c)]] <<- coef[[at(r, c)]] + v
  load_m <- rep(list(S_Zero), n2)
  row <- 0L
  for (node in nodes) {
    if (.pd_has(loads, node[[1L]])) {
      rx <- Symbol(paste0("R_", .pyfmt(node[[1L]]), "_x")); ry <- Symbol(paste0("R_", .pyfmt(node[[1L]]), "_y"))
      for (ld in .pd_get(loads, node[[1L]])) {
        if (!.cm_same(ld[[1L]], rx) && !.cm_same(ld[[1L]], ry)) {
          load_m[[row + 1L]] <- load_m[[row + 1L]] - ld[[1L]] * cos(S_Pi * ld[[2L]] / 180L)
          load_m[[row + 2L]] <- load_m[[row + 2L]] - ld[[1L]] * sin(S_Pi * ld[[2L]] / 180L)
        }
      }
    }
    row <- row + 2L
  }
  cols <- 0L; row <- 0L
  for (node in nodes) {
    t <- stype(node[[1L]])
    if (t == "pinned") {
      add(row, cols, 1L); add(row + 1L, cols + 1L, 1L); cols <- cols + 2L
    } else if (t == "roller") {
      add(row + 1L, cols, 1L); cols <- cols + 1L
    }
    row <- row + 2L
  }
  for (m in mem$values) {
    ps <- .pd_get(nc, m[[1L]]); pe <- .pd_get(nc, m[[2L]])
    len <- sqrt((ps[[1L]] - pe[[1L]])^2 + (ps[[2L]] - pe[[2L]])^2)
    si <- .pl_index(labels, m[[1L]]) - 1L; ei <- .pl_index(labels, m[[2L]]) - 1L
    add(si * 2L, cols, (pe[[1L]] - ps[[1L]]) / len)
    add(si * 2L + 1L, cols, (pe[[2L]] - ps[[2L]]) / len)
    add(ei * 2L, cols, (ps[[1L]] - pe[[1L]]) / len)
    add(ei * 2L + 1L, cols, (ps[[2L]] - pe[[2L]]) / len)
    cols <- cols + 1L
  }
  forces <- (inverse(.new_mat(coef, n2, n2)) %*% .new_mat(load_m, n2, 1L))$data
  symbolic <- function(e) e$type %in% c("Symbol", "Mul", "Add")
  min_load <- NULL
  for (node in nodes) if (.pd_has(loads, node[[1L]]))
    for (ld in .pd_get(loads, node[[1L]])) if (!symbolic(ld[[1L]]))
      if (is.null(min_load) || isTRUE(ld[[1L]] < min_load)) min_load <- ld[[1L]]
  for (j in seq_along(forces)) {
    if (symbolic(forces[[j]])) next
    small <- if (is.null(min_load)) TRUE else isTRUE(abs(as.numeric(N(forces[[j]] / min_load))) < 1e-10)
    if (small) forces[[j]] <- S_Zero
  }
  rl <- .pd_empty(); i <- 1L
  for (node in nodes) {
    t <- stype(node[[1L]])
    if (t == "pinned") {
      rl <- .pd_set(rl, paste0("R_", .pyfmt(node[[1L]]), "_x"), forces[[i]])
      rl <- .pd_set(rl, paste0("R_", .pyfmt(node[[1L]]), "_y"), forces[[i + 1L]])
      i <- i + 2L
    } else if (t == "roller") {
      rl <- .pd_set(rl, paste0("R_", .pyfmt(node[[1L]]), "_y"), forces[[i]])
      i <- i + 1L
    }
  }
  .cm_set(self, "_reaction_loads", rl)
  inf <- .cm_get(self, "_internal_forces")
  for (k in mem$keys) { inf <- .pd_set(inf, k, forces[[i]]); i <- i + 1L }
  .cm_set(self, "_internal_forces", inf)
  invisible(NULL)
}

# a base R drawing of the nodes, members, supports and loads
.truss_draw <- function(self, subs_dict = NULL) {
  nc <- .cm_get(self, "_node_coordinates")
  val <- function(e) {
    if (!is.null(subs_dict)) e <- .cm_dict_subs(e, subs_dict)
    v <- suppressWarnings(as.numeric(N(e)))
    if (is.na(v)) .physics_error("ValueError", "provided substituted dictionary is not adequate")
    v
  }
  xs <- vapply(nc$values, function(p) val(p[[1L]]), 0)
  ys <- vapply(nc$values, function(p) val(p[[2L]]), 0)
  lim <- max(diff(range(xs)), diff(range(ys)), 1)
  graphics::plot(xs, ys, pch = 19, asp = 1, axes = FALSE, xlab = "", ylab = "",
                 xlim = range(xs) + c(-0.1, 0.1) * lim, ylim = range(ys) + c(-0.1, 0.1) * lim)
  mem <- .cm_get(self, "_members")
  for (m in mem$values) {
    i <- .pd_idx(nc, m[[1L]]); j <- .pd_idx(nc, m[[2L]])
    graphics::segments(xs[i], ys[i], xs[j], ys[j], col = "brown", lwd = 3)
  }
  sup <- .cm_get(self, "_supports")
  for (k in seq_along(sup$keys)) {
    i <- .pd_idx(nc, sup$keys[[k]])
    graphics::points(xs[i], ys[i] - 0.03 * lim, pch = if (sup$values[[k]] == "pinned") 17 else 1, cex = 2, col = "blue")
  }
  loads <- .cm_get(self, "_loads")
  for (k in seq_along(loads$keys)) {
    i <- .pd_idx(nc, loads$keys[[k]])
    for (ld in loads$values[[k]]) {
      if (grepl("^R_", .str(ld[[1L]]))) next
      a <- val(ld[[2L]]) * pi / 180
      graphics::arrows(xs[i] - 0.15 * lim * cos(a), ys[i] - 0.15 * lim * sin(a), xs[i], ys[i], col = "orange", length = 0.1)
    }
  }
  graphics::text(xs, ys, vapply(nc$keys, .pyfmt, ""), pos = 3)
  invisible(NULL)
}

local({
  field <- function(f) function(self) .cm_get(self, f)
  .pyclass_def("Truss",
    props = list(nodes = field("_nodes"), node_labels = field("_node_labels"), node_positions = field("_node_positions"),
                 members = field("_members"), member_lengths = field("_member_lengths"), supports = field("_supports"),
                 loads = field("_loads"), reaction_loads = field("_reaction_loads"), internal_forces = field("_internal_forces")),
    methods = list(add_node = .truss_add_node, remove_node = .truss_remove_node, add_member = .truss_add_member,
                   remove_member = .truss_remove_member, change_node_label = .truss_change_node_label,
                   change_member_label = .truss_change_member_label, apply_load = .truss_apply_load,
                   remove_load = .truss_remove_load, apply_support = .truss_apply_support,
                   remove_support = .truss_remove_support, solve = .truss_solve, draw = .truss_draw))
})
