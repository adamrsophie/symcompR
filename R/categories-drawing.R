# Diagram layout and Xy-pic output, after sympy/categories/diagram_drawing.py.
# Grids are 0 indexed as SymPy's; objects and edges are compared by keys.

#' Drawing diagrams
#'
#' `DiagramGrid(diagram, groups, layout, transpose)` lays the objects of a
#' diagram out on a grid as SymPy's `DiagramGrid`; `grid_get(g, i, j)` reads a
#' cell (from 0) and `g$width`, `g$height`, `g$morphisms` describe it.
#' `xypic_draw_diagram(diagram, masked, diagram_format, groups, ...)` returns
#' the Xy-pic code of a diagram. `ArrowStringDescription()` and
#' `XypicDiagramDrawer()` are SymPy's building blocks.
#'
#' @param diagram A diagram.
#' @param groups Groups of objects laid out together.
#' @param layout `"sequential"` for the sequential layout.
#' @param transpose Transpose the grid.
#' @param g A grid.
#' @param i,j Indices from 0.
#' @param masked Morphisms to leave out.
#' @param diagram_format Xy-pic format options.
#' @param ... Hints for `DiagramGrid()`.
#' @param unit,curving,curving_amount,looping_start,looping_end,horizontal_direction,vertical_direction,label_position,label
#'   Fields of an arrow description.
#' @return A grid, a string or a drawer.
#' @examples
#' A <- Object("A"); B <- Object("B"); C <- Object("C")
#' f <- NamedMorphism(A, B, "f"); g <- NamedMorphism(B, C, "g")
#' cat(xypic_draw_diagram(Diagram(list(f, g))))
#' @name diagram_drawing
NULL

.gg_new <- function(width, height) {
  g <- new.env()
  g$width <- width; g$height <- height
  g$rows <- lapply(seq_len(height), function(i) vector("list", width))
  g
}
.gg_get <- function(g, i, j) { r <- g$rows[[i + 1L]]; if (j + 1L > length(r)) NULL else r[[j + 1L]] }
.gg_set <- function(g, i, j, v) { r <- g$rows[[i + 1L]]; if (is.null(v)) r[j + 1L] <- list(NULL) else r[[j + 1L]] <- v; g$rows[[i + 1L]] <- r }
.gg_append_row <- function(g) { g$height <- g$height + 1L; g$rows[[g$height]] <- vector("list", g$width) }
.gg_append_column <- function(g) { g$width <- g$width + 1L; for (i in seq_len(g$height)) { r <- g$rows[[i]]; r[g$width] <- list(NULL); g$rows[[i]] <- r } }
.gg_prepend_row <- function(g) { g$height <- g$height + 1L; g$rows <- c(list(vector("list", g$width)), g$rows) }
.gg_prepend_column <- function(g) { g$width <- g$width + 1L; for (i in seq_len(g$height)) g$rows[[i]] <- c(list(NULL), g$rows[[i]]) }

.okey <- function(o) o$key
.edge <- function(a, b) { ks <- sort(c(a$key, b$key), method = "radix"); objs <- list(a, b)[order(c(a$key, b$key), method = "radix")]
                          list(key = paste(ks, collapse = "||"), objs = objs) }
.mdom <- function(m) morphism_domain(m)
.mcod <- function(m) codomain(m)

# morphisms of a Dict node as an ordered list of (morphism, props)
.dict_items <- function(d) lapply(d$args, function(t) list(m = t$args[[1L]], p = t$args[[2L]]))

.simplify_morphisms <- function(items) Filter(function(it) {
  !((it$m$type == "CompositeMorphism" && it$p$type == "EmptySet") || it$m$type == "IdentityMorphism")
}, items)
.merge_items <- function(a, b) {
  out <- a; keys <- vapply(a, function(it) it$m$key, "")
  for (it in b) { k <- it$m$key; w <- which(keys == k); if (length(w)) out[[w]] <- it else { out[[length(out) + 1L]] <- it; keys <- c(keys, k) } }
  out
}

.juxtapose <- function(e1, e2) {
  k1 <- vapply(e1$objs, .okey, ""); k2 <- vapply(e2$objs, .okey, "")
  inter <- intersect(k1, k2)
  if (length(inter) != 1L) return(NULL)
  rest <- c(e1$objs[k1 != inter], e2$objs[k2 != inter])
  if (length(rest) != 2L || identical(rest[[1L]]$key, rest[[2L]]$key)) {
    if (length(rest) == 1L) return(NULL)
  }
  .edge(rest[[1L]], rest[[2L]])
}

.build_skeleton <- function(morphisms) {
  edges <- list(); eobjs <- list()
  for (m in morphisms) {
    e <- .edge(.mdom(m), .mcod(m))
    if (is.null(edges[[e$key]])) { edges[[e$key]] <- list(m); eobjs[[e$key]] <- e } else edges[[e$key]] <- c(edges[[e$key]], list(m))
  }
  e1 <- eobjs
  for (w in e1) for (v in e1) {
    wv <- .juxtapose(w, v)
    if (!is.null(wv) && is.null(eobjs[[wv$key]])) { edges[[wv$key]] <- list(); eobjs[[wv$key]] <- wv }
  }
  list(morphs = edges, edges = eobjs)
}

.list_triangles <- function(sk) {
  tris <- list()
  for (w in sk$edges) for (v in sk$edges) {
    wv <- .juxtapose(w, v)
    if (!is.null(wv) && !is.null(sk$edges[[wv$key]])) {
      ks <- sort(c(w$key, v$key, wv$key), method = "radix")
      tk <- paste(ks, collapse = "##")
      if (is.null(tris[[tk]])) tris[[tk]] <- list(key = tk, edges = list(sk$edges[[ks[1L]]], sk$edges[[ks[2L]]], sk$edges[[ks[3L]]]))
    }
  }
  unname(tris)
}
.tri_objects <- function(tri) { os <- unlist(lapply(tri$edges, `[[`, "objs"), recursive = FALSE); os[!duplicated(vapply(os, .okey, ""))] }
.mlen <- function(m) if (m$type == "CompositeMorphism") length(components(m)) else 1L

.triangle_key <- function(tri, sk) {
  size <- 0
  for (e in tri$edges) { ms <- sk$morphs[[e$key]]; if (length(ms)) size <- size + max(vapply(ms, .mlen, 0L)) }
  objs <- .key_sort(lapply(.tri_objects(tri), .sort_key), .tri_objects(tri))
  list(size, lapply(objs, .sort_key))
}

.empty_point <- function(pt, g) pt[1L] < 0 || pt[2L] < 0 || pt[1L] >= g$height || pt[2L] >= g$width || is.null(.gg_get(g, pt[1L], pt[2L]))

.put_object <- function(coords, obj, g, fringe) {
  i <- coords[1L]; j <- coords[2L]; offset <- c(0L, 0L)
  if (i == -1L) {
    .gg_prepend_row(g); i <- 0L; offset <- c(1L, 0L)
    for (k in seq_along(fringe$l)) fringe$l[[k]] <- lapply(fringe$l[[k]], function(p) p + c(1L, 0L))
  } else if (i == g$height) .gg_append_row(g)
  if (j == -1L) {
    j <- 0L; offset <- c(offset[1L], 1L); .gg_prepend_column(g)
    for (k in seq_along(fringe$l)) fringe$l[[k]] <- lapply(fringe$l[[k]], function(p) p + c(0L, 1L))
  } else if (j == g$width) .gg_append_column(g)
  .gg_set(g, i, j, obj)
  offset
}

.choose_target_cell <- function(pt1, pt2, edge, obj, sk, g) {
  e1 <- .empty_point(pt1, g); e2 <- .empty_point(pt2, g)
  if (e1 && e2) {
    A <- .gg_get(g, edge[[1L]][1L], edge[[1L]][2L])
    k <- .edge(A, obj)$key
    if (length(sk$morphs[[k]])) return(pt1) else return(pt2)
  }
  if (e1) return(pt1)
  if (e2) return(pt2)
  NULL
}

.find_triangle_to_weld <- function(tris, fringe, g) {
  for (tri in tris) for (fe in fringe$l) {
    k <- .edge(.gg_get(g, fe[[1L]][1L], fe[[1L]][2L]), .gg_get(g, fe[[2L]][1L], fe[[2L]][2L]))$key
    if (k %in% vapply(tri$edges, `[[`, "", "key")) return(list(tri = tri, edge = fe))
  }
  NULL
}

.fringe_remove <- function(fringe, fe) {
  keys <- vapply(fringe$l, function(z) paste(c(z[[1L]], z[[2L]]), collapse = ","), "")
  w <- which(keys == paste(c(fe[[1L]], fe[[2L]]), collapse = ","))[1L]
  fringe$l <- fringe$l[-w]
}

.weld_triangle <- function(tri, we, fringe, g, sk) {
  a <- we[[1L]]; b <- we[[2L]]
  ga <- .gg_get(g, a[1L], a[2L]); gb <- .gg_get(g, b[1L], b[2L])
  obj <- Filter(function(o) !(o$key %in% c(ga$key, gb$key)), .tri_objects(tri))[[1L]]
  target <- NULL
  if (abs(a[1L] - b[1L]) == 1L && abs(a[2L] - b[2L]) == 1L) {
    target <- c(a[1L], b[2L])
    if (!is.null(.gg_get(g, target[1L], target[2L]))) {
      target <- c(b[1L], a[2L])
      if (!is.null(.gg_get(g, target[1L], target[2L]))) { .fringe_remove(fringe, we); return(TRUE) }
    }
  } else if (a[1L] == b[1L]) {
    target <- .choose_target_cell(c(a[1L] + 1L, a[2L]), c(a[1L] + 1L, b[2L]), we, obj, sk, g)
    if (is.null(target)) {
      target <- .choose_target_cell(c(a[1L] - 1L, a[2L]), c(a[1L] - 1L, b[2L]), we, obj, sk, g)
      if (is.null(target)) { .fringe_remove(fringe, we); return(TRUE) }
    }
  } else if (a[2L] == b[2L]) {
    target <- .choose_target_cell(c(a[1L], a[2L] + 1L), c(b[1L], a[2L] + 1L), we, obj, sk, g)
    if (is.null(target)) {
      target <- .choose_target_cell(c(a[1L], a[2L] - 1L), c(b[1L], a[2L] - 1L), we, obj, sk, g)
      if (is.null(target)) { .fringe_remove(fringe, we); return(TRUE) }
    }
  }
  off <- .put_object(target, obj, g, fringe)
  target <- target + off; a <- a + off; b <- b + off
  fringe$l <- c(fringe$l, list(list(a, target), list(b, target)))
  FALSE
}

.pick_root_edge <- function(tri, sk) {
  cand <- Filter(function(e) length(sk$morphs[[e$key]]) > 0L, tri$edges)
  cand <- lapply(cand, function(e) .key_sort(lapply(e$objs, .sort_key), e$objs))
  cand <- .key_sort(lapply(cand, function(c) lapply(c, .sort_key)), cand)
  cand[[1L]]
}

.grow_pseudopod <- function(tris, fringe, g, sk, placed) {
  for (i in seq_len(g$height) - 1L) for (j in seq_len(g$width) - 1L) {
    obj <- .gg_get(g, i, j)
    if (is.null(obj)) next
    good <- Filter(function(tri) {
      ok <- vapply(.tri_objects(tri), .okey, "")
      obj$key %in% ok && !any(setdiff(ok, obj$key) %in% placed)
    }, tris)
    if (!length(good)) next
    tri <- good[[1L]]
    cand <- Filter(function(e) length(sk$morphs[[e$key]]) > 0L, tri$edges)
    cand <- .key_sort(lapply(cand, function(e) .sort_key(FiniteSet(e$objs))), cand)
    edges <- Filter(function(e) obj$key %in% vapply(e$objs, .okey, ""), cand)
    edge <- edges[[1L]]
    other <- Filter(function(o) o$key != obj$key, edge$objs)[[1L]]
    nb <- list(c(i - 1L, j), c(i, j + 1L), c(i + 1L, j), c(i, j - 1L), c(i - 1L, j - 1L), c(i - 1L, j + 1L), c(i + 1L, j - 1L), c(i + 1L, j + 1L))
    for (pt in nb) if (.empty_point(pt, g)) {
      off <- .put_object(pt, other, g, fringe)
      ii <- i + off[1L]; jj <- j + off[2L]; pt <- pt + off
      fringe$l <- c(fringe$l, list(list(c(ii, jj), pt)))
      return(other)
    }
  }
  NULL
}

.generic_layout <- function(diagram, merged) {
  allo <- diagram_objects(diagram)$args
  allk <- vapply(allo, .okey, "")
  if (length(allo) == 1L) { g <- .gg_new(1L, 1L); .gg_set(g, 0L, 0L, allo[[1L]]); return(g) }
  sk <- .build_skeleton(merged)
  g <- .gg_new(2L, 1L)
  if (length(sk$edges) == 1L) {
    objs <- .key_sort(lapply(allo, .sort_key), allo)
    .gg_set(g, 0L, 0L, objs[[1L]]); .gg_set(g, 0L, 1L, objs[[2L]])
    return(g)
  }
  tris <- .list_triangles(sk)
  tris <- Filter(function(tri) sum(vapply(tri$edges, function(e) length(sk$morphs[[e$key]]) > 0L, NA)) >= 2L, tris)
  tris <- .key_sort(lapply(tris, .triangle_key, sk = sk), tris)
  root <- .pick_root_edge(tris[[1L]], sk)
  .gg_set(g, 0L, 0L, root[[1L]]); .gg_set(g, 0L, 1L, root[[2L]])
  fringe <- new.env(); fringe$l <- list(list(c(0L, 0L), c(0L, 1L)))
  placed <- vapply(root, .okey, "")
  while (!setequal(placed, allk)) {
    w <- .find_triangle_to_weld(tris, fringe, g)
    if (!is.null(w)) {
      if (.weld_triangle(w$tri, w$edge, fringe, g, sk)) next
      placed <- union(placed, vapply(.tri_objects(w$tri), .okey, ""))
    } else {
      nobj <- .grow_pseudopod(tris, fringe, g, sk, placed)
      if (is.null(nobj)) {
        rem <- allo[!(allk %in% placed)]
        rg <- DiagramGrid(subdiagram_from_objects(diagram, FiniteSet(rem)))
        fw <- g$width + rg$width; fh <- max(g$height, rg$height)
        fg <- .gg_new(fw, fh)
        # SymPy loops i over the width and j over the height here
        for (i in seq_len(g$width) - 1L) for (j in seq_len(g$height) - 1L)
          if (i < fh && j < fw) .gg_set(fg, i, j, .gg_get(g, i, j))
        for (i in seq_len(rg$height) - 1L) for (j in seq_len(rg$width) - 1L) .gg_set(fg, i, g$width + j, grid_get(rg, i, j))
        return(fg)
      }
      placed <- union(placed, nobj$key)
    }
    tris <- Filter(function(tri) !all(vapply(.tri_objects(tri), .okey, "") %in% placed), tris)
  }
  g
}

.undirected_graph <- function(objects, merged) {
  adj <- stats::setNames(lapply(objects, function(o) list()), vapply(objects, .okey, ""))
  for (m in merged) {
    adj[[.mdom(m)$key]] <- c(adj[[.mdom(m)$key]], list(.mcod(m)))
    adj[[.mcod(m)$key]] <- c(adj[[.mcod(m)$key]], list(.mdom(m)))
  }
  lapply(adj, function(l) .key_sort(lapply(l, .sort_key), l))
}

.sequential_layout <- function(diagram, merged) {
  objects <- diagram_objects(diagram)$args
  so <- .key_sort(lapply(objects, .sort_key), objects)
  adj <- .undirected_graph(objects, merged)
  lens <- vapply(so, function(o) length(adj[[o$key]]), 0L)
  root <- so[[which.min(lens)]]
  g <- .gg_new(1L, 1L); .gg_set(g, 0L, 0L, root)
  placed <- root$key
  dummy <- new.env(); dummy$l <- list()
  place <- function(pt) {
    np <- c(pt[1L], pt[2L] + 1L)
    for (ao in adj[[.gg_get(g, pt[1L], pt[2L])$key]]) {
      if (ao$key %in% placed) next
      .put_object(np, ao, g, dummy)
      placed <<- c(placed, ao$key)
      place(np)
      np <- c(np[1L] + 1L, np[2L])
    }
  }
  place(c(0L, 0L))
  g
}

.connected_components_dg <- function(objects, items) {
  merged <- lapply(items, `[[`, "m")
  adj <- .undirected_graph(objects, merged)
  idx <- stats::setNames(rep(NA_integer_, length(objects)), vapply(objects, .okey, ""))
  objmap <- stats::setNames(objects, names(idx))
  trav <- function(o, ci) { idx[[o$key]] <<- ci; for (n in adj[[o$key]]) if (is.na(idx[[n$key]])) trav(n, ci) }
  ci <- 0L
  for (k in names(adj)) if (is.na(idx[[k]])) { trav(objmap[[k]], ci); ci <- ci + 1L }
  lapply(seq_len(ci) - 1L, function(c) {
    comp <- names(idx)[idx == c]
    cur <- Filter(function(it) .mdom(it$m)$key %in% comp && .mcod(it$m)$key %in% comp, items)
    pairs <- lapply(cur, function(it) list(it$m, it$p))
    if (length(comp) == 1L) pairs[[length(pairs) + 1L]] <- list(IdentityMorphism(objmap[[comp]]), S_EmptySet)
    Diagram(pairs)
  })
}

#' @rdname diagram_drawing
#' @export
DiagramGrid <- function(diagram, groups = NULL, layout = NULL, transpose = FALSE) {
  prem <- .simplify_morphisms(.dict_items(premises(diagram)))
  conc <- .simplify_morphisms(.dict_items(conclusions(diagram)))
  allm <- .merge_items(prem, conc)
  merged <- Filter(function(m) !identical(.mdom(m)$key, .mcod(m)$key), lapply(allm, `[[`, "m"))
  comps <- .connected_components_dg(diagram_objects(diagram)$args, allm)
  same_objs <- is_sym(groups) && identical(groups$key, diagram_objects(diagram)$key)
  if (length(groups) && !same_objs) {
    g <- .handle_groups(diagram, groups, merged, list(layout = layout, transpose = transpose))
  } else if (length(comps) > 1L) {
    comps <- .key_sort(lapply(comps, .sort_key), comps)
    grids <- lapply(comps, function(cm) DiagramGrid(cm, layout = layout, transpose = transpose))
    g <- .gg_new(sum(vapply(grids, function(z) z$width, 0L)), max(vapply(grids, function(z) z$height, 0L)))
    sj <- 0L
    for (gr in grids) {
      for (i in seq_len(gr$height) - 1L) for (j in seq_len(gr$width) - 1L) .gg_set(g, i, sj + j, grid_get(gr, i, j))
      sj <- sj + gr$width
    }
  } else if (!is.null(layout)) {
    g <- if (layout == "sequential") .sequential_layout(diagram, merged) else stop("unknown layout", call. = FALSE)
  } else g <- .generic_layout(diagram, merged)
  if (isTRUE(transpose)) {
    tg <- .gg_new(g$height, g$width)
    for (i in seq_len(g$height) - 1L) for (j in seq_len(g$width) - 1L) .gg_set(tg, j, i, .gg_get(g, i, j))
    g <- tg
  }
  structure(list(grid = g, width = g$width, height = g$height, morphisms = allm), class = "DiagramGrid")
}

#' @rdname diagram_drawing
#' @export
grid_get <- function(g, i, j) .gg_get(g$grid, i, j)

#' @export
format.DiagramGrid <- function(x, ...) {
  rows <- lapply(seq_len(x$height) - 1L, function(i) paste0("[", paste(vapply(seq_len(x$width) - 1L, function(j) {
    o <- grid_get(x, i, j); if (is.null(o)) "None" else .str(o) }, ""), collapse = ", "), "]"))
  paste0("[", paste(rows, collapse = ", "), "]")
}
#' @export
print.DiagramGrid <- function(x, ...) { cat(format(x), "\n"); invisible(x) }

# ---------------------------------------------------------------------------
# Xy-pic

#' @rdname diagram_drawing
#' @export
ArrowStringDescription <- function(unit, curving, curving_amount, looping_start, looping_end,
                                   horizontal_direction, vertical_direction, label_position, label) {
  a <- new.env()
  a$unit <- unit; a$curving <- curving; a$curving_amount <- curving_amount
  a$looping_start <- looping_start; a$looping_end <- looping_end
  a$horizontal_direction <- horizontal_direction; a$vertical_direction <- vertical_direction
  a$label_position <- label_position; a$label <- label
  a$label_displacement <- ""; a$arrow_style <- ""; a$forced_label_position <- FALSE
  class(a) <- "ArrowStringDescription"
  a
}
#' @export
format.ArrowStringDescription <- function(x, ...) {
  # Python truthiness: None and "" are false; None prints as "None"
  tr <- function(v) !is.null(v) && length(v) == 1L && !identical(v, "") && !identical(v, FALSE)
  st <- function(v) if (is.null(v)) "None" else as.character(v)
  cs <- if (tr(x$curving)) sprintf("@/%s%d%s/", x$curving, as.integer(x$curving_amount), st(x$unit)) else ""
  ls <- if (tr(x$looping_start) && tr(x$looping_end)) sprintf("@(%s,%s)", x$looping_start, x$looping_end) else ""
  ss <- if (tr(x$arrow_style)) paste0("@", x$arrow_style) else ""
  sprintf("\\ar%s%s%s[%s%s]%s%s{%s}", cs, ls, ss, st(x$horizontal_direction), st(x$vertical_direction),
          st(x$label_position), st(x$label_displacement), st(x$label))
}

#' @rdname diagram_drawing
#' @export
XypicDiagramDrawer <- function() {
  d <- new.env()
  d$unit <- "mm"; d$default_curving_amount <- 3; d$default_curving_step <- 4
  d$arrow_formatters <- list(); d$default_arrow_formatter <- NULL
  class(d) <- "XypicDiagramDrawer"
  d
}

.xy_coords <- function(oc, o) oc[[o$key]]

.process_loop <- function(i, j, g, info, oc) {
  quad <- c(0, 0, 0, 0); obj <- .gg_get(g, i, j)
  for (it in info) {
    m <- it$m; si <- it$s
    if (identical(.mdom(m)$key, obj$key) && identical(.mcod(m)$key, obj$key)) {
      ls <- paste(si$looping_start, si$looping_end)
      k <- match(ls, c("r u", "u l", "l d", "d r"))
      if (!is.na(k)) quad[k] <- quad[k] + 1
      next
    }
    if (identical(.mdom(m)$key, obj$key)) { e <- .xy_coords(oc, .mcod(m)); out <- TRUE }
    else if (identical(.mcod(m)$key, obj$key)) { e <- .xy_coords(oc, .mdom(m)); out <- FALSE }
    else next
    di <- e[1L] - i; dj <- e[2L] - j; mc <- si$curving
    if (di != 0 && dj != 0) {
      q <- if (di > 0 && dj > 0) 1 else if (di > 0 && dj < 0) 2 else if (di < 0 && dj < 0) 3 else 4
      quad[q] <- quad[q] + 1
    } else if (di == 0) {
      if (dj > 0) { up <- if (out) 1 else 4; lo <- if (out) 4 else 1 } else { up <- if (out) 3 else 2; lo <- if (out) 2 else 3 }
      if (nzchar(mc)) { if (mc == "^") quad[up] <- quad[up] + 1 else if (mc == "_") quad[lo] <- quad[lo] + 1 }
      else { quad[up] <- quad[up] + 1; quad[lo] <- quad[lo] + 1 }
    } else if (dj == 0) {
      if (di < 0) { lf <- if (out) 2 else 1; rt <- if (out) 1 else 2 } else { lf <- if (out) 4 else 3; rt <- if (out) 3 else 4 }
      if (nzchar(mc)) { if (mc == "^") quad[lf] <- quad[lf] + 1 else if (mc == "_") quad[rt] <- quad[rt] + 1 }
      else { quad[lf] <- quad[lf] + 1; quad[rt] <- quad[rt] + 1 }
    }
  }
  fq <- 1
  for (k in 1:4) if (quad[k] < quad[fq]) fq <- k
  lse <- list(c("r", "u"), c("u", "l"), c("l", "d"), c("d", "r"))[[fq]]
  list(curving = "", label_pos = "^", ls = lse[1L], le = lse[2L])
}

.process_straight <- function(i, j, target, g, info, oc, horizontal) {
  backwards <- FALSE
  start <- if (horizontal) j else i; end <- target
  if (end < start) { tmp <- start; start <- end; end <- tmp; backwards <- TRUE }
  side1 <- list(); side2 <- list(); straight <- list()
  if (end - start > 1L) for (k in (start + 1L):(end - 1L)) {
    obj <- if (horizontal) .gg_get(g, i, k) else .gg_get(g, k, j)
    if (is.null(obj)) next
    for (it in info) {
      m <- it$m
      if (identical(.mdom(m)$key, obj$key)) e <- .xy_coords(oc, .mcod(m))
      else if (identical(.mcod(m)$key, obj$key)) e <- .xy_coords(oc, .mdom(m))
      else next
      if (horizontal) {
        if (e[1L] > i) side2 <- c(side2, list(it)) else if (e[1L] < i) side1 <- c(side1, list(it)) else if (!nzchar(it$s$curving)) straight <- c(straight, list(it))
      } else {
        if (e[2L] > j) side2 <- c(side2, list(it)) else if (e[2L] < j) side1 <- c(side1, list(it)) else if (!nzchar(it$s$curving)) straight <- c(straight, list(it))
      }
    }
  }
  setpos <- function(first) for (it in straight) {
    c1 <- .xy_coords(oc, .mdom(it$m)); c2 <- .xy_coords(oc, .mcod(it$m))
    less <- if (horizontal) c1[2L] < c2[2L] else c1[1L] < c2[1L]
    it$s$label_position <- if (less) first else (if (first == "^") "_" else "^")
    it$s$forced_label_position <- TRUE
  }
  if (horizontal) {
    if (length(side1) < length(side2)) { cv <- if (backwards) "_" else "^"; setpos("_") }
    else { cv <- if (backwards) "^" else "_"; setpos("^") }
  } else {
    if (length(side1) < length(side2)) { cv <- if (backwards) "^" else "_"; setpos("^") }
    else { cv <- if (backwards) "_" else "^"; setpos("_") }
  }
  list(curving = cv, label_pos = cv)
}

.process_morphism <- function(d, g, m, oc, info) {
  rep_cond <- function(t, gt, lt) if (t > 0) strrep(gt, t) else strrep(lt, -t)
  c0 <- .xy_coords(oc, .mdom(m)); c1 <- .xy_coords(oc, .mcod(m))
  i <- c0[1L]; j <- c0[2L]; di <- c1[1L] - i; dj <- c1[2L] - j
  vd <- rep_cond(di, "d", "u"); hd <- rep_cond(dj, "r", "l")
  curving <- ""; lpos <- "^"; ls <- ""; le <- ""
  if (di == 0 && dj == 0) { r <- .process_loop(i, j, g, info, oc); curving <- r$curving; lpos <- r$label_pos; ls <- r$ls; le <- r$le }
  else if (di == 0 && abs(dj) > 1) { r <- .process_straight(i, j, c1[2L], g, info, oc, TRUE); curving <- r$curving; lpos <- r$label_pos }
  else if (dj == 0 && abs(di) > 1) { r <- .process_straight(i, j, c1[1L], g, info, oc, FALSE); curving <- r$curving; lpos <- r$label_pos }
  ks <- sort(c(.mdom(m)$key, .mcod(m)$key))
  count <- sum(vapply(info, function(it) identical(sort(unique(c(.mdom(it$m)$key, .mcod(it$m)$key))), sort(unique(ks))), NA))
  amount <- ""
  if (nzchar(curving)) amount <- d$default_curving_amount + count * d$default_curving_step
  else if (count) {
    curving <- "^"
    filt <- sum(vapply(info, function(it) identical(.mdom(it$m)$key, .mdom(m)$key) && identical(.mcod(it$m)$key, .mcod(m)$key) && it$s$curving == curving, NA))
    amount <- d$default_curving_amount + filt * d$default_curving_step
  }
  name <- switch(m$type,
    IdentityMorphism = paste0("id_{%s}", latex(.gg_get(g, i, j))),
    CompositeMorphism = paste(rev(vapply(components(m), function(cp) latex(Symbol(morphism_name(cp))), "")), collapse = "\\circ "),
    NamedMorphism = latex(Symbol(morphism_name(m))), "")
  ArrowStringDescription(d$unit, curving, amount, ls, le, hd, vd, lpos, name)
}

.free_h <- function(di, dj, cj, g) {
  if (dj < cj) { s <- dj; e <- cj; bw <- FALSE } else { s <- cj; e <- dj; bw <- TRUE }
  fu <- if (di == 0) TRUE else all(vapply(s:e, function(j) !is.null(.gg_get(g, di - 1L, j)), NA))
  fd <- if (di == g$height - 1L) TRUE else !any(vapply(s:e, function(j) !is.null(.gg_get(g, di + 1L, j)), NA))
  c(fu, fd, bw)
}
.free_v <- function(di, ci, dj, g) {
  if (di < ci) { s <- di; e <- ci; bw <- FALSE } else { s <- ci; e <- di; bw <- TRUE }
  fl <- if (dj == 0) TRUE else !any(vapply(s:e, function(i) !is.null(.gg_get(g, i, dj - 1L)), NA))
  fr <- if (dj == g$width - 1L) TRUE else !any(vapply(s:e, function(i) !is.null(.gg_get(g, i, dj + 1L)), NA))
  c(fl, fr, bw)
}
.free_d <- function(di, ci, dj, cj, g) {
  rng <- function(a, b) if (a < b) a:b else b:a
  if (di < ci && dj < cj) { si <- di; sj <- dj; ei <- ci; ej <- cj; bw <- FALSE }
  else if (di > ci && dj > cj) { si <- ci; sj <- cj; ei <- di; ej <- dj; bw <- TRUE }
  if (di < ci && dj > cj) { si <- di; sj <- dj; ei <- ci; ej <- cj; bw <- TRUE }
  else if (di > ci && dj < cj) { si <- ci; sj <- cj; ei <- di; ej <- dj; bw <- FALSE }
  alpha <- (ei - si) / (ej - sj)
  fu <- TRUE; fd <- TRUE
  for (i in rng(si, ei)) {
    if (!fu && !fd) break
    for (j in rng(sj, ej)) {
      if (!fu && !fd) break
      if (i == si && j == sj) next
      a1 <- if (j == sj) Inf else (i - si) / (j - sj)
      if (!is.null(.gg_get(g, i, j))) {
        if (is.infinite(a1) || abs(a1) > abs(alpha)) fd <- FALSE
        else if (abs(a1) < abs(alpha)) fu <- FALSE
      }
    }
  }
  c(fu, fd, bw)
}

.push_labels_out <- function(info, g, oc) {
  setp <- function(f1, f2, p1, p2, bw, s) {
    if (bw) { tmp <- p1; p1 <- p2; p2 <- tmp }
    if (f1 && !f2) s$label_position <- p1 else if (f2 && !f1) s$label_position <- p2
  }
  for (it in info) {
    s <- it$s; m <- it$m
    if (nzchar(s$curving) || s$forced_label_position) next
    if (identical(.mdom(m)$key, .mcod(m)$key)) next
    d0 <- .xy_coords(oc, .mdom(m)); c0 <- .xy_coords(oc, .mcod(m))
    if (d0[1L] == c0[1L]) { r <- .free_h(d0[1L], d0[2L], c0[2L], g); setp(r[1L], r[2L], "^", "_", r[3L], s) }
    else if (d0[2L] == c0[2L]) { r <- .free_v(d0[1L], c0[1L], d0[2L], g); setp(r[1L], r[2L], "_", "^", r[3L], s) }
    else { r <- .free_d(d0[1L], c0[1L], d0[2L], c0[2L], g); setp(r[1L], r[2L], "^", "_", r[3L], s) }
  }
}

.morphism_sort_key <- function(m, oc) {
  c0 <- .xy_coords(oc, .mdom(m)); c1 <- .xy_coords(oc, .mcod(m))
  if (identical(.mdom(m)$key, .mcod(m)$key)) return(list(3, 0, .sort_key(m)))
  if (c1[1L] == c0[1L]) return(list(1, abs(c1[2L] - c0[2L]), .sort_key(m)))
  if (c1[2L] == c0[2L]) return(list(1, abs(c1[1L] - c0[1L]), .sort_key(m)))
  list(2, 0, .sort_key(m))
}

.xy_draw <- function(d, diagram, grid, masked = NULL, diagram_format = "") {
  items <- grid$morphisms
  if (length(masked)) { mk <- vapply(masked, `[[`, "", "key"); items <- Filter(function(it) !(it$m$key %in% mk), items) }
  g <- grid$grid
  oc <- list()
  for (i in seq_len(g$height) - 1L) for (j in seq_len(g$width) - 1L) { o <- .gg_get(g, i, j); if (!is.null(o)) oc[[o$key]] <- c(i, j) }
  items <- .key_sort(lapply(items, function(it) .morphism_sort_key(it$m, oc)), items)
  info <- list()
  for (it in items) {
    sd <- .process_morphism(d, g, it$m, oc, info)
    if (!is.null(d$default_arrow_formatter)) d$default_arrow_formatter(sd)
    for (pp in it$p$args) if (!is.null(d$arrow_formatters[[.str(pp)]])) d$arrow_formatters[[.str(pp)]](sd)
    info[[length(info) + 1L]] <- list(m = it$m, s = sd)
  }
  .push_labels_out(info, g, oc)
  objm <- list()
  for (it in info) { k <- .mdom(it$m)$key; objm[[k]] <- c(objm[[k]], list(it)) }
  res <- paste0("\\xymatrix", diagram_format, "{\n")
  for (i in seq_len(g$height) - 1L) {
    for (j in seq_len(g$width) - 1L) {
      o <- .gg_get(g, i, j)
      if (!is.null(o)) {
        res <- paste0(res, latex(o), " ")
        for (it in objm[[o$key]]) res <- paste0(res, format(it$s), " ")
      }
      if (j < g$width - 1L) res <- paste0(res, "& ")
    }
    if (i < g$height - 1L) res <- paste0(res, "\\\\")
    res <- paste0(res, "\n")
  }
  paste0(res, "}\n")
}

#' @rdname diagram_drawing
#' @export
xypic_draw_diagram <- function(diagram, masked = NULL, diagram_format = "", groups = NULL, ...) {
  grid <- DiagramGrid(diagram, groups, ...)
  .xy_draw(XypicDiagramDrawer(), diagram, grid, masked, diagram_format)
}

# DiagramGrid._handle_groups: lay each group out on its own grid, then the
# groups on a grid of groups
.handle_groups <- function(diagram, groups, merged, hints) {
  obj_groups <- list(); grids <- list()
  to_fs <- function(gr) if (is_sym(gr) && gr$type == "FiniteSet") gr else if (is.list(gr) && !is_sym(gr)) FiniteSet(gr) else gr
  lay_out <- function(gr, local) {
    if (is_sym(gr) && gr$type == "FiniteSet") {
      for (o in gr$args) obj_groups[[o$key]] <<- gr
      h <- if (length(local)) local else hints
      grids[[gr$key]] <<- DiagramGrid(subdiagram_from_objects(diagram, gr), layout = h$layout, transpose = isTRUE(h$transpose))
    } else obj_groups[[gr$key]] <<- gr
  }
  dict <- is.list(groups) && length(groups) && all(vapply(groups, function(z) is.list(z) && !is_sym(z) && !is.null(names(z)) && "group" %in% names(z), NA))
  if (dict) for (gh in groups) lay_out(gh$group, gh$hints)
  else for (gr in (if (is_sym(groups)) groups$args else groups)) lay_out(to_fs(gr), NULL)
  nm <- list()
  for (m in merged) {
    dg <- obj_groups[[.mdom(m)$key]]; cg <- obj_groups[[.mcod(m)$key]]
    if (!identical(dg$key, cg$key)) nm[[length(nm) + 1L]] <- NamedMorphism(dg, cg, "dummy")
  }
  top <- DiagramGrid(Diagram(nm))
  gsize <- function(gr) if (!is.null(gr) && !is.null(grids[[gr$key]])) c(grids[[gr$key]]$height, grids[[gr$key]]$width) else c(1L, 1L)
  rh <- vapply(seq_len(top$height) - 1L, function(i) max(vapply(seq_len(top$width) - 1L, function(j) gsize(grid_get(top, i, j))[1L], 0L)), 0L)
  cw <- vapply(seq_len(top$width) - 1L, function(j) max(vapply(seq_len(top$height) - 1L, function(i) gsize(grid_get(top, i, j))[2L], 0L)), 0L)
  g <- .gg_new(sum(cw), sum(rh))
  rr <- 0L
  for (li in seq_len(top$height) - 1L) {
    rc <- 0L
    for (lj in seq_len(top$width) - 1L) {
      o <- grid_get(top, li, lj)
      if (!is.null(o) && !is.null(grids[[o$key]])) {
        lg <- grids[[o$key]]
        for (i in seq_len(lg$height) - 1L) for (j in seq_len(lg$width) - 1L) .gg_set(g, rr + i, rc + j, grid_get(lg, i, j))
      } else .gg_set(g, rr, rc, o)
      rc <- rc + cw[lj + 1L]
    }
    rr <- rr + rh[li + 1L]
  }
  g
}
