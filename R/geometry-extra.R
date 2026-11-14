# Coordinate list helpers shared by integrals-intpoly.R and physics-optics-utils.R.
.vsub <- function(a, b) Map(.sub, a, b)
.vdot <- function(a, b) .add(Map(function(u, v) .mul(list(u, v)), a, b))
.vcross <- function(a, b) list(.sub(.mul(list(a[[2L]], b[[3L]])), .mul(list(a[[3L]], b[[2L]]))),
                               .sub(.mul(list(a[[3L]], b[[1L]])), .mul(list(a[[1L]], b[[3L]]))),
                               .sub(.mul(list(a[[1L]], b[[2L]])), .mul(list(a[[2L]], b[[1L]]))))
.vnorm <- function(a) simplify(.sqrt(.add(lapply(a, .sq))))
.geo_zero <- function(v) is_zero_num(simplify(v))
.ldir <- function(l) .vsub(l$p2$coords, l$p1$coords)
.ldim <- function(l) length(l$p1$coords)
