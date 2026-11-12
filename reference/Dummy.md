# Dummy and wild symbols

`Dummy()` creates a symbol that is unique even when another symbol has
the same name. `Wild()` creates a pattern variable for
[`sym_match()`](https://adamrsophie.github.io/symcompR/reference/sym_match.md)
and
[`sym_replace()`](https://adamrsophie.github.io/symcompR/reference/sym_match.md);
`exclude` lists expressions the wild may not contain and `properties` is
a list of predicate functions the match must satisfy.

## Usage

``` r
Dummy(name = NULL, ...)

Wild(name, exclude = list(), properties = list())
```

## Arguments

- name:

  Display name; by default `Dummy_<n>`, numbered by a session counter as
  SymPy does.

- ...:

  Assumptions, as for
  [`Symbol()`](https://adamrsophie.github.io/symcompR/reference/Symbol.md).

- exclude:

  A list of expressions the match must not contain.

- properties:

  A list of functions returning `TRUE` for valid matches.

## Value

A `sym` object.
