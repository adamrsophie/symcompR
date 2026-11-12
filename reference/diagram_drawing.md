# Drawing diagrams

`DiagramGrid(diagram, groups, layout, transpose)` lays the objects of a
diagram out on a grid as SymPy's `DiagramGrid`; `grid_get(g, i, j)`
reads a cell (from 0) and `g$width`, `g$height`, `g$morphisms` describe
it. `xypic_draw_diagram(diagram, masked, diagram_format, groups, ...)`
returns the Xy-pic code of a diagram. `ArrowStringDescription()` and
`XypicDiagramDrawer()` are SymPy's building blocks.

## Usage

``` r
DiagramGrid(diagram, groups = NULL, layout = NULL, transpose = FALSE)

grid_get(g, i, j)

ArrowStringDescription(
  unit,
  curving,
  curving_amount,
  looping_start,
  looping_end,
  horizontal_direction,
  vertical_direction,
  label_position,
  label
)

XypicDiagramDrawer()

xypic_draw_diagram(
  diagram,
  masked = NULL,
  diagram_format = "",
  groups = NULL,
  ...
)
```

## Arguments

- diagram:

  A diagram.

- groups:

  Groups of objects laid out together.

- layout:

  `"sequential"` for the sequential layout.

- transpose:

  Transpose the grid.

- g:

  A grid.

- i, j:

  Indices from 0.

- unit, curving, curving_amount, looping_start, looping_end,
  horizontal_direction, vertical_direction, label_position, label:

  Fields of an arrow description.

- masked:

  Morphisms to leave out.

- diagram_format:

  Xy-pic format options.

- ...:

  Hints for `DiagramGrid()`.

## Value

A grid, a string or a drawer.

## Examples

``` r
A <- Object("A"); B <- Object("B"); C <- Object("C")
f <- NamedMorphism(A, B, "f"); g <- NamedMorphism(B, C, "g")
cat(xypic_draw_diagram(Diagram(list(f, g))))
#> \xymatrix{
#> A \ar[r]^{f} & B \ar[d]^{g} \\
#> & C 
#> }
```
