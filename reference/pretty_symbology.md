# Pretty printing symbology

The helpers of `sympy.printing.pretty.pretty_symbology` that the pretty
printer draws with. `pretty_symbol(name)` renders a symbol name with
Greek letters and Unicode sub- and superscripts; `xsym(s)` gives the
character for an operator such as `"=="` or `"*"`; `hobj(symb, width)`
and `vobj(symb, height)` draw a horizontal or vertical object such as
`"-"` or `"("`; `center_pad(wstring, wtarget)` the padding that centers
a string, `center(string, width)` the centered string;
`annotated(letter)` the stylized `"F"` or `"G"` of hypergeometric and
Meijer G functions as `list(size, top, bottom, add, picture)`.
`greek_unicode` maps Greek letter names to their characters and `sup`
characters to their superscripts.

## Usage

``` r
pretty_symbol(name, bold_name = FALSE, use_unicode = TRUE)

xsym(s, use_unicode = TRUE)

hobj(symb, width, use_unicode = TRUE)

vobj(symb, height, use_unicode = TRUE)

center_pad(wstring, wtarget, fillchar = " ")

center(string, width, fillchar = " ")

annotated(letter, use_unicode = TRUE)

greek_unicode

sup
```

## Format

`greek_unicode` and `sup` are named lists of strings.

## Arguments

- name:

  A symbol name.

- bold_name:

  Whether to use bold letters.

- use_unicode:

  Whether to use Unicode characters.

- s:

  An operator name.

- symb:

  An object name.

- width, height:

  The size.

- wstring, wtarget:

  The string width and the target width.

- fillchar:

  The padding character.

- string:

  A string.

- letter:

  `"F"` or `"G"`.

## Value

A string, a pair of strings, or a list.

## Examples

``` r
pretty_symbol("alpha_1")
#> [1] "α₁"
xsym("==")
#> [1] "="
hobj("-", 5)
#> [1] "─────"
cat(vobj("(", 3), "\n")
#> ⎛
#> ⎜
#> ⎝ 
center("ab", 6)
#> [1] "  ab  "
```
