# Categories, after sympy/categories/baseclasses.py: objects, morphisms,
# categories and diagrams (with the closure of premises under composition).

#' Categories
#'
#' `Object("A")` is an object; `NamedMorphism(A, B, "f")`,
#' `IdentityMorphism(A)` and `CompositeMorphism(f, g)` (also `g * f`) are
#' morphisms with `morphism_domain()`, `codomain()`, `components()` and
#' `flatten_morphism()`. `Category(name, objects, commutative_diagrams)` is
#' a category. `Diagram(premises, conclusions)` closes its premises under
#' composition and adds identities, as SymPy's `Diagram`; premises and
#' conclusions are lists of morphisms or named lists (a list of pairs
#' `list(morphism, properties)`). `hom()`, `is_subdiagram()` and
#' `subdiagram_from_objects()` are SymPy's methods.
#'
#' @param name A name.
#' @param domain,codomain,A,B Objects.
#' @param ... Morphisms.
#' @param m,x A morphism, category or diagram.
#' @param new_name The name of the flattened morphism.
#' @param objects A set of objects.
#' @param commutative_diagrams Diagrams.
#' @param premises,conclusions Morphisms, or pairs of a morphism and its
#'   properties.
#' @param diagram A diagram.
#' @return A sym object.
#' @examples
#' A <- Object("A"); B <- Object("B"); C <- Object("C")
#' f <- NamedMorphism(A, B, "f"); g <- NamedMorphism(B, C, "g")
#' g * f
#' Diagram(list(f, g))
#' @name categories
NULL

#' @rdname categories
#' @export
Object <- function(name) new_sym("Object", name = as.character(name), key = paste0("Obj:", name))
register_str_printer("Object", function(e) paste0("Object(\"", e$name, "\")"))
register_latex_printer("Object", function(e) .latex(Symbol(e$name)))

.Str <- function(s) new_sym("Str", name = as.character(s), key = paste0("Str:", s))
register_str_printer("Str", function(e) e$name)

.is_morphism <- function(m) is_sym(m) && m$type %in% c("NamedMorphism", "IdentityMorphism", "CompositeMorphism")

#' @rdname categories
#' @export
IdentityMorphism <- function(domain) new_sym("IdentityMorphism", list(domain))
#' @rdname categories
#' @export
NamedMorphism <- function(domain, codomain, name) {
  if (!nzchar(name)) stop("Empty morphism names not allowed.", call. = FALSE)
  new_sym("NamedMorphism", list(domain, codomain, .Str(name)))
}
register_str_printer("IdentityMorphism", function(e) paste0("IdentityMorphism(", .str(e$args[[1L]]), ")"))
register_str_printer("NamedMorphism", function(e) paste0("NamedMorphism(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]),
                                                         ", \"", e$args[[3L]]$name, "\")"))
register_str_printer("CompositeMorphism", function(e) paste0("CompositeMorphism((",
  paste(vapply(e$args[[1L]]$args, .str, ""), collapse = ", "), if (length(e$args[[1L]]$args) == 1L) "," else "", "))"))

#' @rdname categories
#' @export
morphism_domain <- function(m) if (m$type == "CompositeMorphism") morphism_domain(m$args[[1L]]$args[[1L]]) else m$args[[1L]]
#' @rdname categories
#' @export
codomain <- function(m) switch(m$type,
  IdentityMorphism = m$args[[1L]],
  CompositeMorphism = { cs <- m$args[[1L]]$args; codomain(cs[[length(cs)]]) },
  m$args[[2L]])
#' @rdname categories
#' @export
morphism_name <- function(m) m$args[[3L]]$name
#' @rdname categories
#' @export
components <- function(m) if (.is_bd(m)) .vec_components(m) else m$args[[1L]]$args

#' @rdname categories
#' @export
CompositeMorphism <- function(...) {
  cs <- list(...)
  if (length(cs) == 1L && !.is_morphism(cs[[1L]])) cs <- if (is_sym(cs[[1L]])) cs[[1L]]$args else as.list(cs[[1L]])
  add <- function(t, m) {
    if (m$type == "CompositeMorphism") c(t, components(m))
    else if (m$type == "IdentityMorphism") t
    else c(t, list(m))
  }
  norm <- list()
  if (length(cs) > 1L) for (k in seq_len(length(cs) - 1L)) {
    cur <- cs[[k]]; fol <- cs[[k + 1L]]
    if (!.is_morphism(cur) || !.is_morphism(fol)) stop("All components must be morphisms.", call. = FALSE)
    if (!identical(codomain(cur)$key, morphism_domain(fol)$key)) stop("Uncomposable morphisms.", call. = FALSE)
    norm <- add(norm, cur)
  }
  norm <- add(norm, cs[[length(cs)]])
  if (!length(norm)) return(cs[[1L]])
  if (length(norm) == 1L) return(norm[[1L]])
  new_sym("CompositeMorphism", list(do.call(Tuple, norm)))
}

#' @rdname categories
#' @export
flatten_morphism <- function(m, new_name) NamedMorphism(morphism_domain(m), codomain(m), new_name)

#' @rdname categories
#' @export
Category <- function(name, objects = S_EmptySet, commutative_diagrams = list()) {
  if (!nzchar(name)) stop("A Category cannot have an empty name.", call. = FALSE)
  objects <- as_sym(objects)
  cls <- if (objects$type == "Class") objects else new_sym("Class", list(objects))
  new_sym("Category", list(.Str(name), cls, FiniteSet(commutative_diagrams)))
}
register_str_printer("Category", function(e) paste0("Category(\"", e$args[[1L]]$name, "\")"))
#' @rdname categories
#' @export
category_objects <- function(x) x$args[[2L]]
#' @rdname categories
#' @export
commutative_diagrams <- function(x) x$args[[3L]]

# ---------------------------------------------------------------------------
# Diagrams; the morphism dictionaries are environments keyed by sym keys
# that keep insertion order

.dmap <- function() { e <- new.env(); e$keys <- character(); e$m <- list(); e$p <- list(); e }
.dmap_has <- function(d, m) m$key %in% d$keys
.dmap_union <- function(d, m, props) {
  k <- m$key
  if (k %in% d$keys) { d$p[[k]] <- Union(d$p[[k]], props); return(TRUE) }
  d$keys <- c(d$keys, k); d$m[[k]] <- m; d$p[[k]] <- props
  FALSE
}

.add_morphism_closure <- function(ms, m, props, add_identities = TRUE, recurse_composites = TRUE) {
  if (.dmap_union(ms, m, props)) return(invisible())
  if (m$type == "IdentityMorphism") {
    if (props$type != "EmptySet") stop("Instances of IdentityMorphism cannot have properties.", call. = FALSE)
    return(invisible())
  }
  if (add_identities) {
    .dmap_union(ms, IdentityMorphism(morphism_domain(m)), S_EmptySet)
    .dmap_union(ms, IdentityMorphism(codomain(m)), S_EmptySet)
  }
  snap <- lapply(ms$keys, function(k) list(ms$m[[k]], ms$p[[k]]))
  for (ex in snap) {
    em <- ex[[1L]]; np <- Intersection(ex[[2L]], props)
    if (identical(morphism_domain(m)$key, codomain(em)$key)) .dmap_union(ms, CompositeMorphism(em, m), np)
    if (identical(codomain(m)$key, morphism_domain(em)$key)) .dmap_union(ms, CompositeMorphism(m, em), np)
  }
  if (m$type == "CompositeMorphism" && recurse_composites)
    for (cmp in components(m)) .add_morphism_closure(ms, cmp, S_EmptySet, add_identities)
  invisible()
}

.props_set <- function(p) {
  if (is_sym(p) && p$type %in% c("FiniteSet", "EmptySet")) return(p)
  if (is.character(p)) return(FiniteSet(lapply(p, Symbol)))
  if (is.list(p)) return(FiniteSet(lapply(p, function(z) if (is.character(z)) Symbol(z) else as_sym(z))))
  FiniteSet(as_sym(p))
}
.dmap_dict <- function(d) Dict(lapply(d$keys, function(k) d$m[[k]]), lapply(d$keys, function(k) d$p[[k]]))
# a premises argument as a list of (morphism, properties) pairs
.diagram_pairs <- function(a) {
  if (!length(a)) return(list())
  if (is_sym(a) && a$type == "Dict") return(lapply(a$args, function(t) list(t$args[[1L]], t$args[[2L]])))
  if (all(vapply(a, .is_morphism, NA))) return(lapply(a, function(m) list(m, NULL)))
  lapply(a, function(p) list(p[[1L]], p[[2L]]))
}

#' @rdname categories
#' @export
Diagram <- function(premises = list(), conclusions = list()) {
  pm <- .dmap(); cm <- .dmap()
  objects <- S_EmptySet
  for (pp in .diagram_pairs(premises)) {
    m <- pp[[1L]]
    objects <- Union(objects, FiniteSet(morphism_domain(m), codomain(m)))
    .add_morphism_closure(pm, m, if (is.null(pp[[2L]])) S_EmptySet else .props_set(pp[[2L]]))
  }
  for (pp in .diagram_pairs(conclusions)) {
    m <- pp[[1L]]
    if (isTRUE(.rel_to_r(.contains(objects, morphism_domain(m)))) && isTRUE(.rel_to_r(.contains(objects, codomain(m)))))
      .add_morphism_closure(cm, m, if (is.null(pp[[2L]])) S_EmptySet else .props_set(pp[[2L]]),
                            add_identities = FALSE, recurse_composites = FALSE)
  }
  new_sym("Diagram", list(.dmap_dict(pm), .dmap_dict(cm), objects))
}
register_str_printer("Diagram", function(e) paste0("Diagram(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))

#' @rdname categories
#' @export
premises <- function(x) x$args[[1L]]
#' @rdname categories
#' @export
conclusions <- function(x) x$args[[2L]]
#' @rdname categories
#' @export
diagram_objects <- function(x) x$args[[3L]]

.dict_get <- function(d, k) { for (t in d$args) if (identical(t$args[[1L]]$key, k$key)) return(t$args[[2L]]); NULL }

#' @rdname categories
#' @export
hom <- function(x, A, B) {
  pick <- function(d) {
    ms <- Filter(function(m) identical(morphism_domain(m)$key, A$key) && identical(codomain(m)$key, B$key),
                 lapply(d$args, function(t) t$args[[1L]]))
    if (length(ms)) FiniteSet(ms) else S_EmptySet
  }
  if (x$type == "Category") stop("hom-sets are not implemented in Category.", call. = FALSE)
  list(pick(premises(x)), pick(conclusions(x)))
}

#' @rdname categories
#' @export
is_subdiagram <- function(x, diagram) {
  ok <- function(big, small) all(vapply(small$args, function(t) {
    v <- .dict_get(big, t$args[[1L]]); !is.null(v) && identical(v$key, t$args[[2L]]$key)
  }, NA))
  ok(premises(x), premises(diagram)) && ok(conclusions(x), conclusions(diagram))
}

#' @rdname categories
#' @export
subdiagram_from_objects <- function(x, objects) {
  objects <- as_sym(objects)
  if (!isTRUE(is_subset(objects, diagram_objects(x)))) stop("Supplied objects should all belong to the diagram.", call. = FALSE)
  keep <- function(d) Filter(function(t) {
    m <- t$args[[1L]]
    isTRUE(.rel_to_r(.contains(objects, morphism_domain(m)))) && isTRUE(.rel_to_r(.contains(objects, codomain(m))))
  }, d$args)
  np <- lapply(keep(premises(x)), function(t) list(t$args[[1L]], t$args[[2L]]))
  nc <- lapply(keep(conclusions(x)), function(t) list(t$args[[1L]], t$args[[2L]]))
  Diagram(np, nc)
}
register_str_printer("Class", function(e) paste0("Class(", .str(e$args[[1L]]), ")"))
