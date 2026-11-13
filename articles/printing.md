# Printing

symcompR expressions can be printed as plain strings, 2D ASCII or
Unicode art, LaTeX, and MathML. Code generation and parsing have their
own vignettes.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); k <- Symbol("k")
f <- Function("f")
```

## String forms

[`print()`](https://rdrr.io/r/base/print.html) uses
[`sstr()`](https://adamrsophie.github.io/symcompR/reference/sstr.md),
SymPy’s [`str()`](https://rdrr.io/r/utils/str.html).
[`srepr()`](https://adamrsophie.github.io/symcompR/reference/sstr.md)
shows the exact tree structure of the expression.

``` r

sstr(x^2 + 2 * x + 1)
#> [1] "x^2 + 2*x + 1"
srepr(x^2 + 1)
#> [1] "Add(Integer(1), Pow(Symbol('x'), Integer(2)))"
```

## Pretty printing

[`pretty()`](https://adamrsophie.github.io/symcompR/reference/pretty.md)
returns a 2D rendering as a string and
[`pprint()`](https://adamrsophie.github.io/symcompR/reference/pretty.md)
prints it. ASCII output works everywhere; Unicode output uses box
drawing characters.

``` r

cat(pretty((x + 1) / (y - 2), use_unicode = FALSE), "\n")
#> x + 1
#> -----
#> y - 2
cat(pretty(sqrt(x^2 + 1), use_unicode = FALSE), "\n")
#>    ________
#>   /  2     
#> \/  x  + 1
cat(pretty(x^2 + 2 * x / 3, use_unicode = FALSE), "\n")
#>  2   2*x
#> x  + ---
#>       3
cat(pretty((x + 1) / (y - 2), use_unicode = TRUE), "\n")
#> x + 1
#> ─────
#> y - 2
```

## LaTeX and MathML

[`latex()`](https://adamrsophie.github.io/symcompR/reference/latex.md)
produces LaTeX source, ready to paste into a document or to render in R
Markdown with `$$ ... $$`.

``` r

latex(x^2 / 2 + sin(x))
#> [1] "\\frac{x^{2}}{2} + \\sin{\\left(x \\right)}"
latex(sqrt(x + 1) / (y - 2))
#> [1] "\\frac{\\sqrt{x + 1}}{y - 2}"
latex(Rational(3, 4) * sym_pi)
#> [1] "\\frac{3 \\pi}{4}"
latex(Sum(1 / k^2, list(k, 1, oo)))
#> [1] "\\sum_{k=1}^{\\infty} \\frac{1}{k^{2}}"
latex(Derivative(f(x), x, 2))
#> [1] "\\frac{d^{2}}{d x^{2}} f{\\left(x \\right)}"
latex(Eq(x^2 + 1, 0))
#> [1] "x^{2} + 1 = 0"
latex(Matrix(list(list(1, x), list(y, 2))))
#> [1] "\\left[\\begin{matrix}1 & x\\\\y & 2\\end{matrix}\\right]"
latex(Interval(0, 1, TRUE))
#> [1] "\\left(0, 1\\right]"
latex(gamma(x))
#> [1] "\\Gamma\\left(x\\right)"
```

For example, the sum above renders as \sum\_{k=1}^{\infty}
\frac{1}{k^{2}}

[`mathml()`](https://adamrsophie.github.io/symcompR/reference/mathml.md)
produces content MathML by default, or presentation MathML:

``` r

mathml(x^2 + 1)
#> [1] "<apply><plus/><apply><power/><ci>x</ci><cn>2</cn></apply><cn>1</cn></apply>"
mathml(x^2 + 1, printer = "presentation")
#> [1] "<mrow><msup><mi>x</mi><mn>2</mn></msup><mo>+</mo><mn>1</mn></mrow>"
```

[`TableForm()`](https://adamrsophie.github.io/symcompR/reference/TableForm.md)
lays out a table of expressions, as text or as LaTeX:

``` r

TableForm(list(list(5, 7), list(4, 2), list(10, 3)), headings = "automatic")
#>   | 1  2
#> --------
#> 1 | 5  7
#> 2 | 4  2
#> 3 | 10 3
cat(latex(TableForm(list(list(x^2, 1 / x)), headings = list(NULL, c("a", "b")))), "\n")
#> \begin{tabular}{l l}
#> a & b \\
#> \hline
#> $x^{2}$ & $\frac{1}{x}$ \\
#> \end{tabular}
```

[`dotprint()`](https://adamrsophie.github.io/symcompR/reference/mathml.md)
returns a Graphviz description of the expression tree:

``` r

cat(dotprint(x + 2), "\n")
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

## More printing functions

[`print_latex()`](https://adamrsophie.github.io/symcompR/reference/printing-helpers.md),
[`print_mathml()`](https://adamrsophie.github.io/symcompR/reference/mathml.md),
and
[`print_python()`](https://adamrsophie.github.io/symcompR/reference/more-code-printers.md)
print instead of returning a string.
[`python()`](https://adamrsophie.github.io/symcompR/reference/more-code-printers.md)
writes Python code that recreates the expression,
[`multiline_latex()`](https://adamrsophie.github.io/symcompR/reference/more-code-printers.md)
splits a long equation over several lines, and
[`print_tree()`](https://adamrsophie.github.io/symcompR/reference/printing-helpers.md)
shows the expression tree.

``` r

print_latex(x^2 / 2)
#> \frac{x^{2}}{2}
print_mathml(x + 1)
#> <apply><plus/><ci>x</ci><cn>1</cn></apply>
python(Rational(-1, 3))
#> [1] "e = Rational(-1, 3)"
print_python(x^2 + y)
#> x = Symbol('x')
#> y = Symbol('y')
#> e = x**2 + y
cat(multiline_latex(x, x^2 + y + 3), "\n")
#> \begin{align*}
#> x = & x^{2} \\
#> & + y \\
#> & + 3 
#> \end{align*}
print_tree(x + 2 * y)
#> Add
#>   |-Mul
#>   | |-Integer: 2
#>   | +-Symbol: y
#>   +-Symbol: x
sstrrepr(x^2)
#> [1] "x^2"
```

[`pretty_print()`](https://adamrsophie.github.io/symcompR/reference/pretty.md)
is another name for
[`pprint()`](https://adamrsophie.github.io/symcompR/reference/pretty.md),
and
[`pager_print()`](https://adamrsophie.github.io/symcompR/reference/pretty.md)
sends the pretty form through R’s pager, set by `options(pager = )`.
Here the pager is a function that prints the file into the document. The
Unicode setting can be read and changed:

``` r

pretty_print(x^2 + 1, use_unicode = FALSE)
#>  2    
#> x  + 1
pprint_try_use_unicode()
old <- pprint_use_unicode(FALSE)
pprint_use_unicode()
#> [1] FALSE
pprint_use_unicode(old)
#> [1] FALSE
old_pager <- options(pager = function(files, header, title, delete.file) {
  cat(readLines(files), sep = "\n")
  if (delete.file) unlink(files)
})
pager_print(x^2, use_unicode = FALSE)
#>  2
#> x
options(old_pager)
```
