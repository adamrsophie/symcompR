# Categories

`Object("A")` is an object; `NamedMorphism(A, B, "f")`,
`IdentityMorphism(A)` and `CompositeMorphism(f, g)` (also `g * f`) are
morphisms with `morphism_domain()`, `codomain()`, `components()` and
`flatten_morphism()`. `Category(name, objects, commutative_diagrams)` is
a category. `Diagram(premises, conclusions)` closes its premises under
composition and adds identities, as SymPy's `Diagram`; premises and
conclusions are lists of morphisms or named lists (a list of pairs
`list(morphism, properties)`). `hom()`, `is_subdiagram()` and
`subdiagram_from_objects()` are SymPy's methods.

## Usage

``` r
Object(name)

IdentityMorphism(domain)

NamedMorphism(domain, codomain, name)

morphism_domain(m)

codomain(m)

morphism_name(m)

components(m)

CompositeMorphism(...)

flatten_morphism(m, new_name)

Category(name, objects = S_EmptySet, commutative_diagrams = list())

category_objects(x)

commutative_diagrams(x)

Diagram(premises = list(), conclusions = list())

premises(x)

conclusions(x)

diagram_objects(x)

hom(x, A, B)

is_subdiagram(x, diagram)

subdiagram_from_objects(x, objects)
```

## Arguments

- name:

  A name.

- domain, codomain, A, B:

  Objects.

- m, x:

  A morphism, category or diagram.

- ...:

  Morphisms.

- new_name:

  The name of the flattened morphism.

- objects:

  A set of objects.

- commutative_diagrams:

  Diagrams.

- premises, conclusions:

  Morphisms, or pairs of a morphism and its properties.

- diagram:

  A diagram.

## Value

A sym object.

## Examples

``` r
A <- Object("A"); B <- Object("B"); C <- Object("C")
f <- NamedMorphism(A, B, "f"); g <- NamedMorphism(B, C, "g")
g * f
#> CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g")))
Diagram(list(f, g))
#> Diagram({CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g"))): EmptySet, IdentityMorphism(Object("A")): EmptySet, IdentityMorphism(Object("B")): EmptySet, IdentityMorphism(Object("C")): EmptySet, NamedMorphism(Object("A"), Object("B"), "f"): EmptySet, NamedMorphism(Object("B"), Object("C"), "g"): EmptySet}, {}, {Object("A"), Object("B"), Object("C")})
```
