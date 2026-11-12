# MathML and Graphviz output

`mathml(e, printer)` writes content MathML (default) or presentation
MathML (`printer = "presentation"`); `print_mathml()` prints it.
`dotprint(e)` writes the expression tree in the Graphviz dot language,
with SymPy's node naming.

## Usage

``` r
mathml(e, printer = "content")

print_mathml(e, printer = "content")

dotprint(e)
```

## Arguments

- e:

  An expression.

- printer:

  `"content"` or `"presentation"`.

## Value

A character scalar.

## Examples

``` r
x <- Symbol("x")
mathml(x^2 + 1)
#> [1] "<apply><plus/><apply><power/><ci>x</ci><cn>2</cn></apply><cn>1</cn></apply>"
mathml(x^2 + 1, printer = "presentation")
#> [1] "<mrow><msup><mi>x</mi><mn>2</mn></msup><mo>+</mo><mn>1</mn></mrow>"
cat(dotprint(x + 2))
#> digraph{
#> 
#> # Graph style
#> "ordering"="out"
#> "rankdir"="TD"
#> 
#> #########
#> # Nodes #
#> #########
#> 
#> "Add(Integer(2), Symbol('x'))_()" ["color"="black", "label"="Add", "shape"="ellipse"];
#> "Integer(2)_(0,)" ["color"="black", "label"="2", "shape"="ellipse"];
#> "Symbol('x')_(1,)" ["color"="black", "label"="x", "shape"="ellipse"];
#> 
#> #########
#> # Edges #
#> #########
#> 
#> "Add(Integer(2), Symbol('x'))_()" -> "Integer(2)_(0,)";
#> "Add(Integer(2), Symbol('x'))_()" -> "Symbol('x')_(1,)";
#> }
```
