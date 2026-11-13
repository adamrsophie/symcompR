# Inner products

`InnerProduct(bra, ket)` is the unevaluated `<bra|ket>`;
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
evaluates it when the states know their inner product.

## Usage

``` r
InnerProduct(bra, ket)
```

## Arguments

- bra:

  A bra.

- ket:

  A ket.

## Value

An `InnerProduct`.

## Examples

``` r
InnerProduct(Bra("a"), Ket("b"))
#> <a|b>
doit(InnerProduct(OrthogonalBra(1), OrthogonalKet(2)))
#> 0
```
