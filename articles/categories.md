# Categories

``` r

library(symcompR)
```

This vignette covers `sympy.categories`: objects, morphisms, their
composition, diagrams, categories, and the drawing of commutative
diagrams.

## Objects and morphisms

``` r

Ao <- Object("A"); Bo <- Object("B"); Co <- Object("C")
f <- NamedMorphism(Ao, Bo, "f"); g <- NamedMorphism(Bo, Co, "g")
g * f
#> CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g")))
morphism_domain(g * f)
#> Object("A")
codomain(g * f)
#> Object("C")
```

[`IdentityMorphism()`](https://adamrsophie.github.io/symcompR/reference/categories.md)
is the identity of an object and
[`CompositeMorphism()`](https://adamrsophie.github.io/symcompR/reference/categories.md)
a composition of morphisms;
[`components()`](https://adamrsophie.github.io/symcompR/reference/categories.md)
lists the factors and
[`flatten_morphism()`](https://adamrsophie.github.io/symcompR/reference/categories.md)
replaces a composite by a single named morphism.

``` r

Do <- Object("D"); h <- NamedMorphism(Co, Do, "h")
morphism_name(f)
#> [1] "f"
IdentityMorphism(Ao)
#> IdentityMorphism(Object("A"))
codomain(IdentityMorphism(Ao))
#> Object("A")
cm <- CompositeMorphism(f, g)
components(cm)
#> [[1]]
#> NamedMorphism(Object("A"), Object("B"), "f")
#> 
#> [[2]]
#> NamedMorphism(Object("B"), Object("C"), "g")
flatten_morphism(h * g * f, "u")
#> NamedMorphism(Object("A"), Object("D"), "u")
```

## Diagrams and categories

A `Diagram` has premises (the morphisms that are given) and conclusions
(the morphisms that are claimed to exist), each with a set of
properties.

``` r

d <- Diagram(list(f, g))
premises(d)
#> {CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g"))): EmptySet, IdentityMorphism(Object("A")): EmptySet, IdentityMorphism(Object("B")): EmptySet, IdentityMorphism(Object("C")): EmptySet, NamedMorphism(Object("A"), Object("B"), "f"): EmptySet, NamedMorphism(Object("B"), Object("C"), "g"): EmptySet}
conclusions(d)
#> {}
diagram_objects(d)
#> {Object("A"), Object("B"), Object("C")}
hom(d, Ao, Bo)
#> [[1]]
#> {NamedMorphism(Object("A"), Object("B"), "f")}
#> 
#> [[2]]
#> EmptySet
is_subdiagram(Diagram(list(f, g, h)), d)
#> [1] TRUE
d2 <- Diagram(list(list(f, "unique"), list(g, "unique")),
              list(list(g * f, "unique")))
conclusions(d2)
#> {CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g"))): {unique}}
subdiagram_from_objects(d2, FiniteSet(Ao, Bo))
#> Diagram({IdentityMorphism(Object("A")): EmptySet, IdentityMorphism(Object("B")): EmptySet, NamedMorphism(Object("A"), Object("B"), "f"): {unique}}, {}, {Object("A"), Object("B")})
K <- Category("K", FiniteSet(Ao, Bo), list(d))
K
#> Category("K")
category_objects(K)
#> Class({Object("A"), Object("B")})
commutative_diagrams(K)
#> FiniteSet(Diagram({CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g"))): EmptySet, IdentityMorphism(Object("A")): EmptySet, IdentityMorphism(Object("B")): EmptySet, IdentityMorphism(Object("C")): EmptySet, NamedMorphism(Object("A"), Object("B"), "f"): EmptySet, NamedMorphism(Object("B"), Object("C"), "g"): EmptySet}, {}, {Object("A"), Object("B"), Object("C")}))
```

## Drawing diagrams

[`DiagramGrid()`](https://adamrsophie.github.io/symcompR/reference/diagram_drawing.md)
places the objects of a diagram on a grid, and
[`xypic_draw_diagram()`](https://adamrsophie.github.io/symcompR/reference/diagram_drawing.md)
renders the diagram as Xy-pic code for LaTeX.
[`XypicDiagramDrawer()`](https://adamrsophie.github.io/symcompR/reference/diagram_drawing.md)
is the drawer object behind it, and
[`ArrowStringDescription()`](https://adamrsophie.github.io/symcompR/reference/diagram_drawing.md)
describes one arrow.

``` r

grid <- DiagramGrid(Diagram(list(f, g)))
grid$width
#> [1] 2
grid$height
#> [1] 2
grid_get(grid, 1, 1)
#> Object("C")
cat(xypic_draw_diagram(Diagram(list(f, g))), "\n")
#> \xymatrix{
#> A \ar[r]^{f} & B \ar[d]^{g} \\
#> & C 
#> }
#> 
class(XypicDiagramDrawer())
#> [1] "XypicDiagramDrawer"
cat(format(ArrowStringDescription("mm", "", 0, "", "", "d", "r", "_", "f")), "\n")
#> \ar[dr]_{f}
```
