# Ready made models, after sympy/physics/mechanics/models.py.

.pm_syms <- function(prefix, n) lapply(seq_len(n) - 1L, function(i) Symbol(paste0(prefix, i)))

#' Example multibody models
#'
#' `multi_mass_spring_damper()` builds the equations of motion of `n`
#' masses hanging in series from springs and dampers, and
#' `n_link_pendulum_on_cart()` those of an `n` link pendulum on a sliding
#' cart, as SymPy's `sympy.physics.mechanics.models`. Both return a
#' [KanesMethod()] object whose `kanes_equations()` have been formed.
#'
#' @param n Number of masses or links.
#' @param apply_gravity Add gravity to each mass.
#' @param apply_external_forces Add a specified force to each mass.
#' @param cart_force Apply a specified force to the cart.
#' @param joint_torques Apply specified torques at the joints.
#' @return A `KanesMethod` object.
#' @examples
#' kane <- multi_mass_spring_damper(1)
#' kane$mass_matrix
#' kane$forcing
#' @export
multi_mass_spring_damper <- function(n = 1, apply_gravity = FALSE, apply_external_forces = FALSE) {
  n <- as.integer(n)
  mass <- .pm_syms("m", n)
  stiffness <- .pm_syms("k", n)
  damping <- .pm_syms("c", n)
  acceleration_due_to_gravity <- Symbol("g")
  coordinates <- dynamicsymbols(sprintf("x:%d", n))
  speeds <- dynamicsymbols(sprintf("v:%d", n))
  specifieds <- dynamicsymbols(sprintf("f:%d", n))
  ceiling <- ReferenceFrame("N")
  origin <- PhysPoint("origin")
  origin$set_vel(ceiling, 0)
  points <- list(origin)
  kinematic_equations <- list()
  particles <- list()
  forces <- list()
  for (i in seq_len(n)) {
    last <- points[[length(points)]]
    center <- last$locatenew(sprintf("center%d", i - 1L), coordinates[[i]] * ceiling$x)
    center$set_vel(ceiling, .pv_point_vel(last, ceiling) + speeds[[i]] * ceiling$x)
    points[[length(points) + 1L]] <- center
    block <- Particle(sprintf("block%d", i - 1L), center, mass[[i]])
    kinematic_equations[[i]] <- .sub(speeds[[i]], .pv_diff(coordinates[[i]], .pv_t()))
    total_force <- .sub(.neg(.mul(list(stiffness[[i]], coordinates[[i]]))), .mul(list(damping[[i]], speeds[[i]])))
    if (i < n)
      total_force <- .add(list(total_force, .mul(list(stiffness[[i + 1L]], coordinates[[i + 1L]])),
                               .mul(list(damping[[i + 1L]], speeds[[i + 1L]]))))
    if (isTRUE(apply_gravity)) total_force <- .add(list(total_force, .mul(list(mass[[i]], acceleration_due_to_gravity))))
    if (isTRUE(apply_external_forces)) total_force <- .add(list(total_force, specifieds[[i]]))
    forces[[i]] <- .pytuple(list(center, total_force * ceiling$x))
    particles[[i]] <- block
  }
  kane <- KanesMethod(ceiling, q_ind = coordinates, u_ind = speeds, kd_eqs = kinematic_equations)
  kane$kanes_equations(particles, forces)
  kane
}

#' @rdname multi_mass_spring_damper
#' @export
n_link_pendulum_on_cart <- function(n = 1, cart_force = TRUE, joint_torques = FALSE) {
  if (n <= 0) .pv_err("ValueError", "The number of links must be a positive integer.")
  n <- as.integer(n)
  q <- dynamicsymbols(sprintf("q:%d", n + 1L))
  u <- dynamicsymbols(sprintf("u:%d", n + 1L))
  if (isTRUE(joint_torques)) Tq <- dynamicsymbols(sprintf("T1:%d", n + 1L))
  m <- .pm_syms("m", n + 1L)
  l <- .pm_syms("l", n)
  g <- Symbol("g"); t <- Symbol("t")
  I <- ReferenceFrame("I")
  O <- PhysPoint("O")
  O$set_vel(I, 0)
  P0 <- PhysPoint("P0")
  P0$set_pos(O, q[[1L]] * I$x)
  P0$set_vel(I, u[[1L]] * I$x)
  Pa0 <- Particle("Pa0", P0, m[[1L]])
  frames <- list(I)
  points <- list(P0)
  particles <- list(Pa0)
  forces <- list(.pytuple(list(P0, .neg(.mul(list(m[[1L]], g))) * I$y)))
  kindiffs <- list(.sub(.pv_diff(q[[1L]], t), u[[1L]]))
  specified <- if (isTRUE(cart_force) || isTRUE(joint_torques)) list() else NULL
  for (i in seq_len(n)) {
    Bi <- I$orientnew(sprintf("B%d", i - 1L), "Axis", list(q[[i + 1L]], I$z))
    Bi$set_ang_vel(I, u[[i + 1L]] * I$z)
    frames[[length(frames) + 1L]] <- Bi
    Pi <- points[[length(points)]]$locatenew(sprintf("P%d", i), l[[i]] * Bi$y)
    Pi$v2pt_theory(points[[length(points)]], I, Bi)
    points[[length(points) + 1L]] <- Pi
    particles[[length(particles) + 1L]] <- Particle(paste0("Pa", i), Pi, m[[i + 1L]])
    forces[[length(forces) + 1L]] <- .pytuple(list(Pi, .neg(.mul(list(m[[i + 1L]], g))) * I$y))
    if (isTRUE(joint_torques)) {
      specified[[length(specified) + 1L]] <- Tq[[i]]
      if (i == 1L) forces[[length(forces) + 1L]] <- .pytuple(list(I, .neg(Tq[[i]]) * I$z))
      if (i == n) forces[[length(forces) + 1L]] <- .pytuple(list(Bi, Tq[[i]] * I$z))
      else forces[[length(forces) + 1L]] <- .pytuple(list(Bi, Tq[[i]] * I$z - Tq[[i + 1L]] * I$z))
    }
    kindiffs[[length(kindiffs) + 1L]] <- .sub(.pv_diff(q[[i + 1L]], t), u[[i + 1L]])
  }
  if (isTRUE(cart_force)) {
    F <- dynamicsymbols("F")
    forces[[length(forces) + 1L]] <- .pytuple(list(P0, F * I$x))
    specified[[length(specified) + 1L]] <- F
  }
  kane <- KanesMethod(I, q_ind = q, u_ind = u, kd_eqs = kindiffs)
  kane$kanes_equations(particles, forces)
  kane
}
