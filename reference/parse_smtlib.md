# Parse SMT-LIB

`parse_smtlib(s)` parses the text of an SMT-LIB script, as
`sympy.parsing.smtlib.parse_smtlib()`, and returns the declared symbols
and the assertions. Declaring a constant of sort `Int` or `Real` adds
the assumption `Q$integer()` or `Q$real()` to the assertions. Terms
built from the core, integer and real theories become expressions (`ite`
gives a `Piecewise`, `div` and `mod` are Euclidean), `define-fun` macros
are expanded and `push`, `pop`, `reset` and `reset-assertions` act on
the state. Invalid syntax signals an `SMTLibSyntaxError`, an unknown
command an `UnknownSMTLibCommandError`, an unknown operator an
`UnknownSMTLibOperatorError`, and quantifiers, datatypes, recursive or
parametric definitions and the bitvector, string and array operators a
`NotImplementedError`.

## Usage

``` r
parse_smtlib(s)
```

## Arguments

- s:

  A character string with the SMT-LIB text.

## Value

A list with `symbols`, a named list of the symbols by name, and
`assertions`, a list of expressions and predicates.

## Examples

``` r
r <- parse_smtlib("(declare-const x Int) (declare-const y Int) (assert (= (+ x y) 10))")
r$symbols$x
#> x
r$assertions
#> [[1]]
#> Q.integer(x)
#> 
#> [[2]]
#> Q.integer(y)
#> 
#> [[3]]
#> Eq(x + y, 10)
#> 
```
