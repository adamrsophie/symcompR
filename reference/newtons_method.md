# Newton's method as an abstract syntax tree

`newtons_method(expr, wrt)` builds the
[`CodeBlock()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
of a Newton iteration solving `expr == 0` for `wrt`, as SymPy's
`newtons_method()`. `newtons_method_function(expr, wrt)` wraps it in a
[`FunctionDefinition()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
returning the root.

## Usage

``` r
newtons_method(
  expr,
  wrt,
  atol = 1e-12,
  delta = NULL,
  rtol = 4e-16,
  debug = FALSE,
  itermax = NULL,
  counter = NULL,
  delta_fn = function(e, x) -e/diff(e, x),
  cse = FALSE,
  handle_nan = NULL,
  bounds = NULL
)

newtons_method_function(
  expr,
  wrt,
  params = NULL,
  func_name = "newton",
  attrs = list(),
  delta = NULL,
  ...
)
```

## Arguments

- expr:

  The expression whose root is sought.

- wrt:

  The symbol to solve for.

- atol:

  The absolute tolerance.

- delta:

  `NULL` (a
  [`Dummy()`](https://adamrsophie.github.io/symcompR/reference/Dummy.md)
  inside a
  [`Scope()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)),
  or the symbol of the step.

- rtol:

  The relative tolerance.

- debug:

  Whether to print the iterates.

- itermax:

  `NULL`, or the largest number of iterations.

- counter:

  `NULL`, or the symbol counting the iterations.

- delta_fn:

  An R function of the expression and the variable giving the step; the
  default is the Newton step `-e/diff(e, x)`.

- cse:

  Whether to apply common subexpression elimination to the step.

- handle_nan:

  `NULL`, or a node executed when the step is not a number.

- bounds:

  `NULL`, or a list of the lower and upper bounds of `wrt`.

- params:

  `NULL` (only `wrt`), or a list of the parameters (symbols,
  [`Variable()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
  or
  [`Pointer()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
  nodes).

- func_name:

  The name of the function.

- attrs:

  The attributes of the function.

- ...:

  Further arguments of `newtons_method()`.

## Value

A
[`CodeBlock()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
(in a
[`Scope()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md)
when `delta` is `NULL`), or for `newtons_method_function()` a
[`FunctionDefinition()`](https://adamrsophie.github.io/symcompR/reference/codegen-ast.md).

## Examples

``` r
x <- Symbol("x", real = TRUE)
algo <- newtons_method(cos(x) - x^3, x)
cat(ccode(newtons_method_function(cos(x) - x^3, x)), "\n")
#> double newton(double x){
#>    double d_x = INFINITY;
#>    while (fabs(d_x) > 3.9999999999999999e-16*fabs(x) + 9.9999999999999998e-13) {
#>       d_x = (pow(x, 3) - cos(x))/(-3*pow(x, 2) - sin(x));
#>       x += d_x;
#>    };
#>    return x;
#> } 
```
