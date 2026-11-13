# Physical units

symcompR ports `sympy.physics.units`: units and physical constants are
symbolic quantities that combine with ordinary expressions, convert
between each other exactly, and carry dimensions that can be checked.

``` r

library(symcompR)
```

## The unit registry

[`physics_units()`](https://adamrsophie.github.io/symcompR/reference/physics_units.md)
returns a list of all predefined units, prefixes, constants, dimensions,
and unit systems. Index it with `$`, or evaluate code inside it with
[`with()`](https://rdrr.io/r/base/with.html).

``` r

u <- physics_units()
5 * u$kilometer
#> 5*kilometer
u$kilo * u$meter
#> 1000*meter
u$kilo * u$kilo
#> M
with(u, 3 * newton * 2 * meter)
#> 6*meter*newton
```

Prefixes multiply into units, and two prefixes combine into a single one
(`kilo * kilo` is mega, printed `M`).

## Converting units

`convert_to(expr, target)` rewrites an expression in the target units.
The conversion is exact: the factor between a mile and a kilometer is
the rational number 1.609344.

``` r

convert_to(u$mile, u$kilometer)
#> 25146*kilometer/15625
N(convert_to(u$mile, u$kilometer))
#> 1.609344*kilometer
convert_to(u$inch, u$centimeter)
#> 127*centimeter/50
convert_to(100 * u$kilometer / u$hour, u$meter / u$second)
#> 250*meter/(9*second)
convert_to(u$day, u$second)
#> 86400*second
convert_to(u$liter, u$meter^3)
#> meter^3/1000
convert_to(u$atmosphere, u$pascal)
#> 101325*pascal
convert_to(u$watt * u$hour, u$joule)
#> 3600*joule
convert_to(u$kilo * u$gram, u$pound)
#> 100000000*pound/45359237
convert_to(3 * u$newton * 2 * u$meter, u$joule)
#> 6*joule
convert_to(u$ohm * u$ampere, u$volt)
#> volt
```

The target may be a list of units, in which case the result is expressed
in terms of all of them. This is a convenient way to see a derived unit
in base units.

``` r

convert_to(u$newton, list(u$kilogram, u$meter, u$second))
#> kilogram*meter/second^2
convert_to(u$volt, list(u$kilogram, u$meter, u$second, u$ampere))
#> kilogram*meter^2/(ampere*second^3)
```

## Physical constants

The constants use their CODATA values; the SI defining constants are
exact.

``` r

convert_to(u$speed_of_light, list(u$meter, u$second))
#> 299792458*meter/second
convert_to(u$planck, list(u$joule, u$second))
#> 6.62607015e-34*joule*second
convert_to(u$hbar, list(u$joule, u$second))
#> 3.313035075e-34*joule*second/pi
convert_to(u$elementary_charge, u$coulomb)
#> 1.602176634e-19*coulomb
convert_to(u$boltzmann, list(u$joule, u$kelvin))
#> 1.380649e-23*joule/kelvin
convert_to(u$gravitational_constant, list(u$meter, u$kilogram, u$second))
#> 6.6743e-11*meter^3/(kilogram*second^2)
convert_to(u$joule, u$electronvolt)
#> 6.24150907446076e+18*electronvolt
```

A light year is the speed of light times a year:

``` r

convert_to(u$speed_of_light * u$year, u$meter)
#> 9.46052817742682e+15*meter
```

## Simplifying quantities

[`quantity_simplify()`](https://adamrsophie.github.io/symcompR/reference/convert_to.md)
combines prefixes and units of the same dimension.

``` r

quantity_simplify(u$kilo * u$foot * u$inch)
#> 250*foot^2/3
```

## Dimensions

Every quantity has a dimension in a unit system. Dimensions form their
own algebra, and the SI dimension system expresses each one in terms of
the base dimensions.

``` r

u$length / u$time
#> Dimension(length/time)
dimension_name(u$length)
#> length
dimension_symbol(u$length)
#> L
Dimension("length", "L") == u$length
#> [1] TRUE
get_quantity_dimension(u$SI, u$newton)
#> Dimension(force, F)
quantity_dimension(u$newton)
#> Dimension(force, F)
get_quantity_scale_factor(u$SI, u$kilometer)
#> 1000
equivalent_dims(u$dimsys_SI, u$energy, u$force * u$length)
#> [1] TRUE
SI <- get_unit_system("SI")
get_dimensional_expr(SI, u$newton / u$meter)
#> force/length
```

[`collect_factor_and_dimension()`](https://adamrsophie.github.io/symcompR/reference/UnitSystem.md)
returns the numeric scale factor of an expression together with its
dimension. Adding quantities of different dimensions is an error, which
makes this a dimensional analysis check.

``` r

collect_factor_and_dimension(SI, u$meter + u$kilometer)
#> [[1]]
#> 1001
#> 
#> [[2]]
#> Dimension(length, L)
collect_factor_and_dimension(SI, 5 * u$kilometer / (2 * u$hour))
#> [[1]]
#> 25/36
#> 
#> [[2]]
#> Dimension(length/time)
tryCatch(collect_factor_and_dimension(SI, u$meter + u$second),
         error = function(e) conditionMessage(e))
#> [1] "Dimension of \"second\" is Dimension(time, T), but it should be Dimension(length, L)"
```

## Defining new units

[`Quantity()`](https://adamrsophie.github.io/symcompR/reference/Quantity.md)
creates a new unit, and
[`set_quantity_dimension()`](https://adamrsophie.github.io/symcompR/reference/set_quantity_dimension.md)
and
[`set_quantity_scale_factor()`](https://adamrsophie.github.io/symcompR/reference/set_quantity_dimension.md)
register it in a unit system.

``` r

furlong <- Quantity("furlong")
set_quantity_dimension(u$SI, furlong, u$length)
set_quantity_scale_factor(u$SI, furlong, 220 * u$yard)
convert_to(furlong, u$meter)
#> 25146*meter/125
```

[`Unit()`](https://adamrsophie.github.io/symcompR/reference/Quantity.md)
is another name for
[`Quantity()`](https://adamrsophie.github.io/symcompR/reference/Quantity.md).
[`set_global_dimension()`](https://adamrsophie.github.io/symcompR/reference/Quantity.md)
and
[`set_global_relative_scale_factor()`](https://adamrsophie.github.io/symcompR/reference/Quantity.md)
register a quantity in every unit system at once, relative to another
quantity.

``` r

rod <- Unit("rod")
set_global_relative_scale_factor(rod, Rational(33, 2), u$foot)
convert_to(rod, u$meter)
#> 12573*meter/2500
cfm <- Quantity("cfm")
set_global_dimension(cfm, u$length^3 / u$time)
quantity_name(rod)
#> rod
quantity_abbrev(u$meter)
#> m
is_prefixed(u$kilometer)
#> [1] TRUE
is_prefixed(u$meter)
#> [1] FALSE
scale_factor(u$kilometer)
#> 1000
```

[`PhysicalConstant()`](https://adamrsophie.github.io/symcompR/reference/Quantity.md)
creates a constant, and
[`is_physical_constant()`](https://adamrsophie.github.io/symcompR/reference/Quantity.md)
tests for one.

``` r

is_physical_constant(u$speed_of_light)
#> [1] TRUE
is_physical_constant(u$meter)
#> [1] FALSE
k_e <- PhysicalConstant("coulomb_constant", abbrev = "k_e")
is_physical_constant(k_e)
#> [1] TRUE
```

## Prefixes

[`Prefix()`](https://adamrsophie.github.io/symcompR/reference/Prefix.md)
defines a new prefix and
[`prefix_unit()`](https://adamrsophie.github.io/symcompR/reference/Prefix.md)
builds the prefixed versions of a unit.

``` r

u$kilo
#> k
scale_factor(u$kilo)
#> 1000
Prefix("zz", "Z", 2)
#> Z
prefix_unit(u$gram, list(k = u$kilo, m = u$milli))
#> [[1]]
#> kilogram
#> 
#> [[2]]
#> milligram
```

## Dimension systems

A `DimensionSystem` lists base dimensions and derived dimensions with
their dependencies on the base ones.

``` r

ds <- DimensionSystem(list(u$length, u$mass, u$time),
                      list(Dimension("velocity"), Dimension("force")),
                      list(velocity = list(length = 1, time = -1),
                           force = list(mass = 1, length = 1, time = -2)))
base_dims(ds)
#> [[1]]
#> Dimension(length, L)
#> 
#> [[2]]
#> Dimension(mass, M)
#> 
#> [[3]]
#> Dimension(time, T)
derived_dims(ds)
#> [[1]]
#> Dimension(force)
#> 
#> [[2]]
#> Dimension(velocity)
list_can_dims(ds)
#> [[1]]
#> Dimension(length, L)
#> 
#> [[2]]
#> Dimension(mass, M)
#> 
#> [[3]]
#> Dimension(time, T)
can_transf_matrix(ds)
#> Matrix(3, 3)
#> [1, 0, 0]
#> [0, 1, 0]
#> [0, 0, 1]
inv_can_transf_matrix(ds)
#> Matrix(3, 3)
#> [1, 0, 0]
#> [0, 1, 0]
#> [0, 0, 1]
dim_can_vector(ds, u$length)
#> Matrix(3, 1)
#> [1]
#> [0]
#> [0]
dim_vector(ds, u$length)
#> Matrix(3, 1)
#> [1]
#> [0]
#> [0]
print_dim_base(ds, Dimension("force"))
#> L*M/T^2
dimensional_dependencies(ds)
#> {Dimension(force): {Dimension(length): 1, Dimension(mass): 1, Dimension(time): -2}, Dimension(length, L): {Dimension(length, L): 1}, Dimension(mass, M): {Dimension(mass, M): 1}, Dimension(time, T): {Dimension(time, T): 1}, Dimension(velocity): {Dimension(length): 1, Dimension(time): -1}}
is_consistent(ds)
#> [1] TRUE
system_dim(ds)
#> [1] 3
```

The SI dimension system is `u$dimsys_SI`:

``` r

get_dimensional_dependencies(u$dimsys_SI, u$force)
#> {Dimension(mass): 1, Dimension(length): 1, Dimension(time): -2}
is_dimensionless(u$dimsys_SI, u$length / u$length)
#> [1] TRUE
has_integer_powers(u$length, u$dimsys_SI)
#> [1] TRUE
```

## Unit systems

A `UnitSystem` chooses base units and a dimension system.

``` r

SI
#> SI
get_default_unit_system()
#> SI
get_dimension_system(SI)
#> DimensionSystem((Dimension(amount_of_substance), Dimension(current, I), Dimension(length, L), Dimension(luminous_intensity), Dimension(mass, M), Dimension(temperature, T), Dimension(time, T)), (Dimension(acceleration), Dimension(action), Dimension(area), Dimension(capacitance), Dimension(charge), Dimension(conductance), Dimension(energy), Dimension(force), Dimension(frequency), Dimension(impedance), Dimension(inductance), Dimension(magnetic_density), Dimension(magnetic_flux), Dimension(momentum), Dimension(power), Dimension(pressure), Dimension(velocity), Dimension(voltage), Dimension(volume)), {Dimension(acceleration): {Dimension(length): 1, Dimension(time): -2}, Dimension(action): {Dimension(length): 2, Dimension(mass): 1, Dimension(time): -1}, Dimension(amount_of_substance): {Dimension(amount_of_substance): 1}, Dimension(area): {Dimension(length): 2}, Dimension(capacitance): {Dimension(current): 2, Dimension(length): -2, Dimension(mass): -1, Dimension(time): 4}, Dimension(charge): {Dimension(current): 1, Dimension(time): 1}, Dimension(conductance): {Dimension(current): 2, Dimension(length): -2, Dimension(mass): -1, Dimension(time): 3}, Dimension(current, I): {Dimension(current, I): 1}, Dimension(energy): {Dimension(length): 2, Dimension(mass): 1, Dimension(time): -2}, Dimension(force): {Dimension(length): 1, Dimension(mass): 1, Dimension(time): -2}, Dimension(frequency): {Dimension(time): -1}, Dimension(impedance): {Dimension(current): -2, Dimension(length): 2, Dimension(mass): 1, Dimension(time): -3}, Dimension(inductance): {Dimension(current): -2, Dimension(length): 2, Dimension(mass): 1, Dimension(time): -2}, Dimension(length, L): {Dimension(length, L): 1}, Dimension(luminous_intensity): {Dimension(luminous_intensity): 1}, Dimension(magnetic_density): {Dimension(current): -1, Dimension(mass): 1, Dimension(time): -2}, Dimension(magnetic_flux): {Dimension(current): -1, Dimension(length): 2, Dimension(mass): 1, Dimension(time): -2}, Dimension(mass, M): {Dimension(mass, M): 1}, Dimension(momentum): {Dimension(length): 1, Dimension(mass): 1, Dimension(time): -1}, Dimension(power): {Dimension(length): 2, Dimension(mass): 1, Dimension(time): -3}, Dimension(pressure): {Dimension(length): -1, Dimension(mass): 1, Dimension(time): -2}, Dimension(temperature, T): {Dimension(temperature, T): 1}, Dimension(time, T): {Dimension(time, T): 1}, Dimension(velocity): {Dimension(length): 1, Dimension(time): -1}, Dimension(voltage): {Dimension(current): -1, Dimension(length): 2, Dimension(mass): 1, Dimension(time): -3}, Dimension(volume): {Dimension(length): 3}})
system_dim(SI)
#> [1] 7
is_consistent(SI)
#> [1] TRUE
length(derived_units(SI)) > 0
#> [1] TRUE
any(vapply(get_units_non_prefixed(SI), function(q) q == u$volt, TRUE))
#> [1] TRUE
UnitSystem(list(u$meter, u$second), name = "MS")
#> MS
extend(UnitSystem(list(u$meter, u$second), list(u$speed_of_light)),
       list(u$kilogram))
#> UnitSystem((meter, second, kilogram))
```

## Checking and finding units

[`check_dimensions()`](https://adamrsophie.github.io/symcompR/reference/convert_to.md)
signals an error when an expression adds quantities of different
dimensions and returns the expression otherwise.
[`find_unit()`](https://adamrsophie.github.io/symcompR/reference/convert_to.md)
lists the predefined units with a given name or dimension.

``` r

check_dimensions(u$inch + u$meter)
#> inch + meter
tryCatch(check_dimensions(u$meter + u$second),
         error = function(e) conditionMessage(e))
#> [1] "addends have incompatible dimensions: (((Dimension(length, L), 1),), ((Dimension(time, T), 1),))"
find_unit("charge")
#> [1] "C"                 "coulomb"           "coulombs"         
#> [4] "planck_charge"     "elementary_charge"
find_unit(u$inch)
#>  [1] "m"                  "au"                 "cm"                
#>  [4] "dm"                 "ft"                 "km"                
#>  [7] "ly"                 "mi"                 "mm"                
#> [10] "nm"                 "pm"                 "um"                
#> [13] "yd"                 "nmi"                "feet"              
#> [16] "foot"               "inch"               "mile"              
#> [19] "yard"               "meter"              "miles"             
#> [22] "yards"              "inches"             "meters"            
#> [25] "micron"             "microns"            "angstrom"          
#> [28] "angstroms"          "decimeter"          "kilometer"         
#> [31] "lightyear"          "nanometer"          "picometer"         
#> [34] "centimeter"         "decimeters"         "kilometers"        
#> [37] "lightyears"         "micrometer"         "millimeter"        
#> [40] "nanometers"         "picometers"         "centimeters"       
#> [43] "micrometers"        "millimeters"        "nautical_mile"     
#> [46] "planck_length"      "nautical_miles"     "astronomical_unit" 
#> [49] "astronomical_units"
```
