# Classical mechanics

This vignette covers `sympy.physics.vector` (reference frames, vectors,
dyadics and points), `sympy.physics.mechanics` (bodies, loads, joints,
pathways, actuators and the equations of motion) and
`sympy.physics.biomechanics` (muscle models). The classes keep their
state and are used through `$`, as fields (`N$x`) and methods
(`A$dcm(N)`). `Point`, `Vector` and `Dyadic` are exported as
[`PhysPoint()`](https://adamrsophie.github.io/symcompR/reference/PhysPoint.md),
[`PhysVector()`](https://adamrsophie.github.io/symcompR/reference/PhysVector.md)
and
[`PhysDyadic()`](https://adamrsophie.github.io/symcompR/reference/PhysVector.md)
to avoid clashes with the geometry and vector modules.

``` r

library(symcompR)
x <- Symbol("x"); m <- Symbol("m"); k <- Symbol("k")
r <- Symbol("r"); cc <- Symbol("c")
```

## Reference frames and vectors

[`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md)
creates a frame with unit vectors `$x`, `$y` and `$z`. Indexing a frame
gives its coordinate symbols
([`CoordinateSym()`](https://adamrsophie.github.io/symcompR/reference/CoordinateSym.md)),
with R’s indices starting at 1. Vectors support `+`, scalar `*`, `&`
(dot product), `^` (cross product) and `|` (outer product, giving a
dyadic).

``` r

N <- ReferenceFrame("N")
N$x
#> N.x
N$name
#> [1] "N"
ReferenceFrame("A")[1]
#> A_x
CoordinateSym("A_x", ReferenceFrame("A"), 0)
#> A_x
v <- x * N$x + 2 * N$y
v
#> x*N.x + 2*N.y
v$magnitude()
#> sqrt(x^2 + 4)
N$x ^ N$y
#> N.z
(N$x + N$y) & N$y
#> 1
N$x | N$y
#> (N.x|N.y)
(N$x | N$y) & N$y
#> N.x
PhysVector(0)
#> 0
PhysDyadic(0)
#> 0
```

[`dynamicsymbols()`](https://adamrsophie.github.io/symcompR/reference/dynamicsymbols.md)
creates functions of time `t`, the generalized coordinates and speeds of
mechanics.

``` r

q <- dynamicsymbols("q"); u <- dynamicsymbols("u")
qd <- dynamicsymbols("q", 1); ud <- dynamicsymbols("u", 1)
qd
#> Derivative(q(t), t)
dynamicsymbols("q1:4")
#> [[1]]
#> q1(t)
#> 
#> [[2]]
#> q2(t)
#> 
#> [[3]]
#> q3(t)
(x * N$x + q * N$y)$diff(q, N)
#> N.y
(q * N$x)$dt(N)
#> Derivative(q(t), t)*N.x
```

Frames are oriented relative to each other; the direction cosine matrix
and the angular velocity follow.

``` r

A <- N$orientnew("A", "Axis", list(q, N$z))
A$dcm(N)
#> Matrix(3, 3)
#> [cos(q(t)) , sin(q(t)), 0]
#> [-sin(q(t)), cos(q(t)), 0]
#> [0         , 0        , 1]
A$ang_vel_in(N)
#> Derivative(q(t), t)*N.z
time_derivative(N$orientnew("B", "Axis", list(q, N$x))$y, N)
#> Derivative(q(t), t)*B.z
kinematic_equations(dynamicsymbols("u1 u2 u3"), dynamicsymbols("q1 q2 q3"), "body", "313")
#> [[1]]
#> -(u1(t)*sin(q3(t)) + u2(t)*cos(q3(t)))/sin(q2(t)) + Derivative(q1(t), t)
#> 
#> [[2]]
#> -u1(t)*cos(q3(t)) + u2(t)*sin(q3(t)) + Derivative(q2(t), t)
#> 
#> [[3]]
#> (u1(t)*sin(q3(t)) + u2(t)*cos(q3(t)))*cos(q2(t))/sin(q2(t)) - u3(t) + Derivative(q3(t), t)
get_motion_params(N, position = Symbol("t")^2 * N$x)
#> (2*N.x, 2*t*N.x, t^2*N.x)
partial_velocity(list(dynamicsymbols("u1") * N$x + dynamicsymbols("u2") * N$y),
                 dynamicsymbols("u1 u2"), N)
#> [[1]]
#> [[1]][[1]]
#> N.x
#> 
#> [[1]][[2]]
#> N.y
```

Points have positions and velocities relative to other points.

``` r

O <- PhysPoint("O")
P <- O$locatenew("P", q * N$x)
O$set_vel(N, 0)
P$vel(N)
#> Derivative(q(t), t)*N.x
P$pos_from(O)
#> q(t)*N.x
```

## Printing

The vector printers write time derivatives with primes or dots and drop
the `(t)`.
[`vsprint()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md)
and
[`msprint()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md)
return a string,
[`vsstrrepr()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md)
its representation,
[`vlatex()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md)
and
[`mlatex()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md)
LaTeX, and
[`vprint()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md),
[`mprint()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md),
[`vpprint()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md)
and
[`mpprint()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md)
print directly.

``` r

vsprint(qd * N$x)
#> [1] "q'*N.x"
vsprint(qd)
#> [1] "q'"
vsstrrepr(qd)
#> [1] "q'"
msprint(qd)
#> [1] "q'"
vlatex(dynamicsymbols("q", 2))
#> [1] "\\ddot{q}"
mlatex(qd)
#> [1] "\\dot{q}"
VectorLatexPrinter()$doprint(qd)
#> [1] "\\dot{q}"
vprint(qd * N$x)
#> q'*N.x
mprint(qd)
#> q'
vpprint(qd * N$x)
#> q̇ n_x
mpprint(qd^2)
#>  2
#> q̇
```

[`init_vprinting()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md)
(also
[`mechanics_printing()`](https://adamrsophie.github.io/symcompR/reference/vsprint.md))
makes this the default printing for the rest of the R session, so it is
not run here:

``` r

init_vprinting()
mechanics_printing()
```

## Particles, rigid bodies and loads

[`inertia()`](https://adamrsophie.github.io/symcompR/reference/inertia.md)
builds an inertia dyadic, and
[`Inertia()`](https://adamrsophie.github.io/symcompR/reference/inertia.md)
pairs it with a reference point.

``` r

inertia(N, 1, 2, 3)
#> (N.x|N.x) + 2*(N.y|N.y) + 3*(N.z|N.z)
inertia_of_point_mass(m, x * N$x, N)
#> m*x^2*(N.y|N.y) + m*x^2*(N.z|N.z)
Inertia(inertia(N, 1, 2, 3), PhysPoint("O"))
#> Inertia(dyadic=(N.x|N.x) + 2*(N.y|N.y) + 3*(N.z|N.z), point=O)
Inertia_from_inertia_scalars(PhysPoint("O"), N, 1, 2, 3)
#> Inertia(dyadic=(N.x|N.x) + 2*(N.y|N.y) + 3*(N.z|N.z), point=O)
```

[`Particle()`](https://adamrsophie.github.io/symcompR/reference/Particle.md)
and
[`RigidBody()`](https://adamrsophie.github.io/symcompR/reference/Particle.md)
carry mass; the energy and momentum functions accept any number of them.

``` r

P3 <- PhysPoint("P")
P3$set_vel(N, 3 * N$x)
pa <- Particle("p", P3, 2)
pa$linear_momentum(N)
#> 6*N.x
pa$kinetic_energy(N)
#> 9
linear_momentum(N, pa)
#> 6*N.x
kinetic_energy(N, pa)
#> 9
pa$potential_energy <- m * Symbol("g") * q
potential_energy(pa)
#> g*m*q(t)
Lagrangian(N, pa)
#> -g*m*q(t) + 9
RigidBody("B")$central_inertia
#> B_ixx*(B_frame.x|B_frame.x) + B_ixy*(B_frame.x|B_frame.y) + B_izx*(B_frame.x|B_frame.z) + B_ixy*(B_frame.y|B_frame.x) + B_iyy*(B_frame.y|B_frame.y) + B_iyz*(B_frame.y|B_frame.z) + B_izx*(B_frame.z|B_frame.x) + B_iyz*(B_frame.z|B_frame.y) + B_izz*(B_frame.z|B_frame.z)
RigidBody("B")$mass
#> B_mass
O2 <- PhysPoint("O2")
Pm <- O2$locatenew("Pm", r * N$x)
Pm$set_vel(N, r * u * N$y)
angular_momentum(O2, N, Particle("pm", Pm, m))
#> m*r^2*u(t)*N.z
Oc <- PhysPoint("Oc")
pc1 <- Particle("pc1", Oc$locatenew("Pc1", N$x), 1)
pc2 <- Particle("pc2", Oc$locatenew("Pc2", 3 * N$x), 1)
center_of_mass(Oc, pc1, pc2)
#> 2*N.x
```

Loads are forces at points and torques on frames:

``` r

Force(P3, 3 * N$x)
#> Force(point=P, force=3*N.x)
Torque(N, 2 * N$z)
#> Torque(frame=N, torque=2*N.z)
gravity(-9.81 * N$y, pa)
#> [[1]]
#> Force(point=P, force=- 19.62*N.y)
```

[`find_dynamicsymbols()`](https://adamrsophie.github.io/symcompR/reference/find_dynamicsymbols.md)
lists the dynamic symbols of an expression, and
[`msubs()`](https://adamrsophie.github.io/symcompR/reference/msubs.md)
substitutes into it.

``` r

sapply(find_dynamicsymbols(q^2 + qd), sstr)
#> [1] "Derivative(q(t), t)" "q(t)"
msubs(q^2 + qd, list(list(q, 2)))
#> Derivative(q(t), t) + 4
```

## Pathways and actuators

A pathway connects attachment points: a straight line
([`LinearPathway()`](https://adamrsophie.github.io/symcompR/reference/LinearPathway.md)),
a path through intermediate points
([`ObstacleSetPathway()`](https://adamrsophie.github.io/symcompR/reference/LinearPathway.md)),
or a path wrapping around a geometry
([`WrappingPathway()`](https://adamrsophie.github.io/symcompR/reference/LinearPathway.md)
with
[`WrappingCylinder()`](https://adamrsophie.github.io/symcompR/reference/WrappingSphere.md),
[`WrappingSphere()`](https://adamrsophie.github.io/symcompR/reference/WrappingSphere.md)
or
[`WrappingCone()`](https://adamrsophie.github.io/symcompR/reference/WrappingSphere.md)).

``` r

pA <- PhysPoint("pA"); pB <- PhysPoint("pB")
pB$set_pos(pA, q * N$x)
path <- LinearPathway(pA, pB)
path$length
#> sqrt(q(t)^2)
path$extension_velocity
#> sqrt(q(t)^2)*Derivative(q(t), t)/q(t)
path$to_loads(Symbol("F"))
#> [[1]]
#> Force(point=pA, force=- F*q(t)/sqrt(q(t)^2)*N.x)
#> 
#> [[2]]
#> Force(point=pB, force=F*q(t)/sqrt(q(t)^2)*N.x)
O3 <- PhysPoint("O3"); p1 <- PhysPoint("p1"); p2 <- PhysPoint("p2")
p1$set_pos(O3, r * N$x); p2$set_pos(O3, r * N$y)
WrappingPathway(p1, p2, WrappingCylinder(r, O3, N$z))$length
#> pi*sqrt(r^2)/2
WrappingCylinder(r, PhysPoint("O"), N$x)
#> WrappingCylinder(radius=r, point=O, axis=N.x)
WrappingSphere(r, PhysPoint("O"))
#> WrappingSphere(radius=r, point=O)
WrappingCone(Symbol("alpha"), PhysPoint("O"), N$z)
#> WrappingCone(alpha=alpha, apex=O, axis=N.z)
o1 <- PhysPoint("o1"); o2 <- PhysPoint("o2"); o3 <- PhysPoint("o3")
o2$set_pos(o1, N$x); o3$set_pos(o2, N$y)
ObstacleSetPathway(o1, o2, o3)$length
#> 2
```

Actuators produce loads along a pathway: springs, dampers, friction, and
general forces and torques.

``` r

LinearSpring(k, path)$force
#> -k*sqrt(q(t)^2)
LinearSpring(k, path)$to_loads()
#> [[1]]
#> Force(point=pA, force=k*q(t)*N.x)
#> 
#> [[2]]
#> Force(point=pB, force=- k*q(t)*N.x)
LinearDamper(cc, path)$force
#> -c*sqrt(q(t)^2)*Derivative(q(t), t)/q(t)
DuffingSpring(k, Symbol("beta"), path)$force
#> -beta*(q(t)^2)^(3/2) - k*sqrt(q(t)^2)
CoulombKineticFriction(Symbol("mu_k"), Symbol("f_n"), path)$force
#> -f_n*mu_k*sign(sqrt(q(t)^2)*Derivative(q(t), t)/q(t))
ForceActuator(Symbol("F"), path)$to_loads()
#> [[1]]
#> Force(point=pA, force=- F*q(t)/sqrt(q(t)^2)*N.x)
#> 
#> [[2]]
#> Force(point=pB, force=F*q(t)/sqrt(q(t)^2)*N.x)
TorqueActuator(Symbol("T"), N$z, ReferenceFrame("A"), N)$to_loads()
#> [[1]]
#> Torque(frame=A, torque=T*N.z)
#> 
#> [[2]]
#> Torque(frame=N, torque=- T*N.z)
```

## Joints

Joints connect a parent and a child body and create the coordinates,
speeds and kinematic equations of the connection.

``` r

Pp <- RigidBody("P"); Cp <- RigidBody("C")
J <- PinJoint("J", Pp, Cp)
J$coordinates
#> Matrix(1, 1)
#> [q_J(t)]
J$kdes
#> Matrix(1, 1)
#> [u_J(t) - Derivative(q_J(t), t)]
Cp$frame$dcm(Pp$frame)
#> Matrix(3, 3)
#> [1, 0           , 0          ]
#> [0, cos(q_J(t)) , sin(q_J(t))]
#> [0, -sin(q_J(t)), cos(q_J(t))]
TorqueActuator_at_pin_joint(Symbol("T"), J)$to_loads()
#> [[1]]
#> Torque(frame=C_frame, torque=T*P_frame.x)
#> 
#> [[2]]
#> Torque(frame=P_frame, torque=- T*P_frame.x)
PrismaticJoint("S", RigidBody("P2"), RigidBody("C2"))$coordinates
#> Matrix(1, 1)
#> [q_S(t)]
CylindricalJoint("Y", RigidBody("P3"), RigidBody("C3"))$coordinates
#> Matrix(2, 1)
#> [q0_Y(t)]
#> [q1_Y(t)]
PlanarJoint("L", RigidBody("P4"), RigidBody("C4"))$coordinates
#> Matrix(3, 1)
#> [q0_L(t)]
#> [q1_L(t)]
#> [q2_L(t)]
SphericalJoint("Sp", RigidBody("P5"), RigidBody("C5"))$coordinates
#> Matrix(3, 1)
#> [q0_Sp(t)]
#> [q1_Sp(t)]
#> [q2_Sp(t)]
WeldJoint("W", RigidBody("P6"), RigidBody("C6"))$coordinates
#> Matrix(0, 0)
#> Matrix(0, 0, [])
```

[`Joint()`](https://adamrsophie.github.io/symcompR/reference/Joint.md)
is the abstract base class of the joints. The deprecated
[`Body()`](https://adamrsophie.github.io/symcompR/reference/Body.md) and
[`JointsMethod()`](https://adamrsophie.github.io/symcompR/reference/JointsMethod.md)
still work, with a deprecation warning, but
[`RigidBody()`](https://adamrsophie.github.io/symcompR/reference/Particle.md),
[`Particle()`](https://adamrsophie.github.io/symcompR/reference/Particle.md)
and
[`System()`](https://adamrsophie.github.io/symcompR/reference/System.md)
replace them.

``` r

tryCatch(Joint("J0", Pp, Cp), error = conditionMessage)
#> [1] "Can't instantiate abstract class Joint without an implementation for abstract methods '_generate_coordinates', '_generate_speeds', '_orient_frames', '_set_angular_velocity', '_set_linear_velocity'"
B <- suppressWarnings(Body("B"))
B$apply_force(B$x)
B$loads
#> [[1]]
#> (B_masscenter, B_frame.x)
suppressWarnings(JointsMethod(Pp, J))$q
#> Matrix(1, 1)
#> [q_J(t)]
```

## Equations of motion

Kane’s method for a mass on a spring:

``` r

Pk <- PhysPoint("Pk")
Pk$set_vel(N, u * N$x)
km <- KanesMethod(N, list(q), list(u), kd_eqs = list(qd - u))
invisible(km$kanes_equations(list(Particle("p", Pk, m)), list(Force(Pk, -k * q * N$x))))
km$mass_matrix
#> Matrix(1, 1)
#> [m]
km$forcing
#> Matrix(1, 1)
#> [-k*q(t)]
km$rhs()
#> Matrix(2, 1)
#> [u(t)     ]
#> [-k*q(t)/m]
```

Lagrange’s method takes the Lagrangian:

``` r

LagrangesMethod(m * qd^2 / 2 - k * q^2 / 2, list(q))$form_lagranges_equations()
#> Matrix(1, 1)
#> [k*q(t) + m*Derivative(q(t), (t, 2))]
```

[`Linearizer()`](https://adamrsophie.github.io/symcompR/reference/Linearizer.md)
linearizes the equations, either built from the parts or taken from a
method object:

``` r

km$to_linearizer()$linearize(A_and_B = TRUE)
#> (Matrix([[0, 1], [-k/m, 0]]), Matrix(0, 0, []))
Linearizer(list(qd), list(-u), list(ud), list(k * q), list(0),
           list(), list(), list(), list(q), list(u),
           q_i = list(q), u_i = list(u))$linearize(A_and_B = TRUE)
#> (Matrix([[0, 1], [-k, 0]]), Matrix(0, 0, []))
```

[`System()`](https://adamrsophie.github.io/symcompR/reference/System.md)
collects the coordinates, speeds, bodies, loads and kinematic equations;
[`SymbolicSystem()`](https://adamrsophie.github.io/symcompR/reference/SymbolicSystem.md)
stores equations of motion already derived.

``` r

sy <- System()
sy$add_coordinates(q)
sy$add_speeds(u)
sy$add_kdes(u - qd)
sy$q
#> Matrix(1, 1)
#> [q(t)]
sy$kdes
#> Matrix(1, 1)
#> [u(t) - Derivative(q(t), t)]
System_from_newtonian(RigidBody("Nb"))$frame
#> Nb_frame
ss <- SymbolicSystem(list(q, u), Matrix(list(u, -q)))
ss$states
#> Matrix(2, 1)
#> [q(t)]
#> [u(t)]
ss$comb_explicit_rhs
#> Matrix(2, 1)
#> [u(t) ]
#> [-q(t)]
```

Two ready made models:

``` r

msd <- multi_mass_spring_damper(1)
msd$mass_matrix
#> Matrix(1, 1)
#> [m0]
msd$forcing
#> Matrix(1, 1)
#> [-c0*v0(t) - k0*x0(t)]
cart <- n_link_pendulum_on_cart(1)
cart$mass_matrix_full
#> Matrix(4, 4)
#> [1, 0, 0                , 0                ]
#> [0, 1, 0                , 0                ]
#> [0, 0, m0 + m1          , -l0*m1*cos(q1(t))]
#> [0, 0, -l0*m1*cos(q1(t)), l0^2*m1          ]
cart$forcing_full
#> Matrix(4, 1)
#> [u0(t)                           ]
#> [u1(t)                           ]
#> [-l0*m1*u1(t)^2*sin(q1(t)) + F(t)]
#> [g*l0*m1*sin(q1(t))              ]
```

## Biomechanics

Activation dynamics turn an excitation into a muscle activation.

``` r

act <- FirstOrderActivationDeGroote2016("bicep")
act$x
#> Matrix(1, 1)
#> [a_bicep(t)]
act$constants
#> Matrix(3, 1)
#> [tau_a_bicep]
#> [tau_d_bicep]
#> [b_bicep    ]
act$rhs()
#> Matrix(1, 1)
#> [((1/2 - tanh(b_bicep*(-a_bicep(t) + e_bicep(t)))/2)*(3*a_bicep(t)/2 + 1/2)/tau_d_bicep + (tanh(b_bicep*(-a_bicep(t) + e_bicep(t)))/2 + 1/2)/(tau_a_bicep*(3*a_bicep(t)/2 + 1/2)))*(-a_bicep(t) + e_bicep(t))]
FirstOrderActivationDeGroote2016_with_defaults("bicep")$rhs()
#> Matrix(1, 1)
#> [((1/2 - tanh(10.0*a_bicep(t) - 10.0*e_bicep(t))/2)/(0.0225*a_bicep(t) + 0.0075) + 16.6666666666667*(3*a_bicep(t)/2 + 1/2)*(tanh(10.0*a_bicep(t) - 10.0*e_bicep(t))/2 + 1/2))*(-a_bicep(t) + e_bicep(t))]
ZerothOrderActivation("bicep")$r
#> Matrix(1, 1)
#> [e_bicep(t)]
ZerothOrderActivation_with_defaults("bicep")$order
#> [1] 0
```

The characteristic curves of De Groote et al. (2016) relate the tendon
and fiber forces to their lengths and velocities. Each curve takes its
constants as arguments, and the `_with_defaults` forms use the published
values.

``` r

lT <- Symbol("l_T_tilde")
fT <- TendonForceLengthDeGroote2016_with_defaults(lT)
fT
#> TendonForceLengthDeGroote2016(l_T_tilde, 0.2, 0.995, 0.25, 33.93669377311689)
doit(fT)
#> 4.327167934027389e-16*exp(33.93669377311689*l_T_tilde) - 0.25
diff(fT, lT)
#> 6.787338754623378*exp(33.93669377311689*(l_T_tilde - 0.995))
identical(inverse(fT), TendonForceLengthInverseDeGroote2016)
#> [1] TRUE
doit(TendonForceLengthDeGroote2016(x, Symbol("c0"), Symbol("c1"), Symbol("c2"), Symbol("c3")))
#> c0*exp(c3*(-c1 + x)) - c2
doit(TendonForceLengthInverseDeGroote2016_with_defaults(x))
#> 0.02946663003430684*log(5.0*x + 1.25) + 0.995
doit(TendonForceLengthInverseDeGroote2016(x, 1, 2, 3, 4))
#> log(x + 3)/4 + 2
doit(FiberForceLengthPassiveDeGroote2016_with_defaults(x))
#> 2.37439874427164e-5*exp(6.66666666666667*x) - 0.0186573603637741
doit(FiberForceLengthPassiveDeGroote2016(x, 1, 2))
#> (exp(2*x - 2) - 1)/(-1 + exp(2))
doit(FiberForceLengthPassiveInverseDeGroote2016_with_defaults(x))
#> 0.15*log(53.5981500331442*x + 1) + 1
doit(FiberForceLengthPassiveInverseDeGroote2016(x, 1, 2))
#> log(x*(-1 + exp(2)) + 1)/2 + 1
FiberForceLengthActiveDeGroote2016_with_defaults(x)
#> FiberForceLengthActiveDeGroote2016(x, 0.814, 1.06, 0.162, 0.0633, 0.433, 0.717, -0.0299, 0.2, 0.1, 1.0, 0.354, 0.0)
FiberForceLengthActiveDeGroote2016(x, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
#> FiberForceLengthActiveDeGroote2016(x, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
doit(FiberForceVelocityDeGroote2016_with_defaults(x))
#> 0.886 - 0.318*log(-8.149*x + sqrt((-8.149*x - 0.374)^2 + 1) - 0.374)
doit(FiberForceVelocityDeGroote2016(x, 1, 2, 3, 4))
#> log(2*x + sqrt((2*x + 3)^2 + 1) + 3) + 4
doit(FiberForceVelocityInverseDeGroote2016_with_defaults(x))
#> 0.122714443489999*sinh(3.14465408805031*x - 2.78616352201258) - 0.0458952018652595
doit(FiberForceVelocityInverseDeGroote2016(x, 1, 2, 3, 4))
#> sinh(x - 4)/2 - 3/2
curves <- CharacteristicCurveCollection(
  TendonForceLengthDeGroote2016, TendonForceLengthInverseDeGroote2016,
  FiberForceLengthPassiveDeGroote2016, FiberForceLengthPassiveInverseDeGroote2016,
  FiberForceLengthActiveDeGroote2016, FiberForceVelocityDeGroote2016,
  FiberForceVelocityInverseDeGroote2016)
curves$tendon_force_length(Symbol("l"), 1, 2, 3, 4)
#> TendonForceLengthDeGroote2016(l, 1, 2, 3, 4)
```

A musculotendon model combines a pathway, activation dynamics and the
curves;
[`MusculotendonFormulation()`](https://adamrsophie.github.io/symcompR/reference/MusculotendonFormulation.md)
selects its state variables.

``` r

MusculotendonFormulation("RIGID_TENDON")
#> <MusculotendonFormulation.RIGID_TENDON: 0>
MusculotendonFormulation(2) == 2
#> [1] TRUE
Om <- PhysPoint("Om"); Im <- PhysPoint("Im")
Im$set_pos(Om, q * N$x)
muscle <- MusculotendonDeGroote2016("m", LinearPathway(Om, Im),
                                    FirstOrderActivationDeGroote2016("m"))
muscle$x
#> Matrix(1, 1)
#> [a_m(t)]
muscle$r
#> Matrix(1, 1)
#> [e_m(t)]
MusculotendonDeGroote2016_with_defaults("m", LinearPathway(Om, Im),
                                        FirstOrderActivationDeGroote2016("m"))$constants
#> Matrix(6, 1)
#> [l_T_slack_m]
#> [F_M_max_m  ]
#> [l_M_opt_m  ]
#> [tau_a_m    ]
#> [tau_d_m    ]
#> [b_m        ]
```

[`ActivationBase()`](https://adamrsophie.github.io/symcompR/reference/activation-dynamics.md),
[`CharacteristicCurveFunction()`](https://adamrsophie.github.io/symcompR/reference/characteristic-curves.md)
and
[`MusculotendonBase()`](https://adamrsophie.github.io/symcompR/reference/musculotendon.md)
are abstract base classes, which cannot be instantiated:

``` r

tryCatch(ActivationBase("a"), error = conditionMessage)
#> [1] "Can't instantiate abstract class ActivationBase without an implementation for abstract methods 'F', 'M', 'constants', 'input_vars', 'order', 'p', 'r', 'rhs', 'state_vars', 'with_defaults', 'x'"
tryCatch(CharacteristicCurveFunction(), error = conditionMessage)
#> [1] "Cannot directly instantiate 'CharacteristicCurveFunction', instances of characteristic curves must be of a concrete subclass."
tryCatch(MusculotendonBase("m"), error = conditionMessage)
#> [1] "Can't instantiate abstract class MusculotendonBase without an implementation for abstract method 'curves'"
```
