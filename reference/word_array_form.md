# Words in free groups

Methods of SymPy's `FreeGroupElement`. A word `w` is made from the
generators returned by
[`free_group()`](https://adamrsophie.github.io/symcompR/reference/FreeGroup.md).
`word_array_form()` gives the pairs (generator, exponent),
`word_letter_form()` the letters (an inverse letter is a negated
symbol), `word_ext_rep()` the flattened array form and
`word_letter_form_elm()` the letters as words. `word_len()` is the
length, `word_inverse()` the inverse, `word_order()` the order (1 or
`oo`) and `word_commutator(w, v)` is `w^-1*v^-1*w*v`.
`is_identity_word()` tests for the identity. The other functions are
named as in SymPy; positions are zero based and ranges exclude their
end, as in Python.

## Usage

``` r
word_array_form(w)

word_letter_form(w)

word_letter_form_elm(w)

word_ext_rep(w)

word_len(w)

is_identity_word(w)

word_inverse(w)

word_order(w)

word_commutator(w, v)

word_prod(G, words)

word_dict(words, by)

eliminate_words(w, words, all = FALSE, inverse = TRUE)

eliminate_word(w, gen, by = NULL, all = FALSE, inverse = TRUE)

subword_index(w, word, start = 0L)

is_dependent(w, word)

is_independent(w, word)

subword(w, from_i, to_j, strict = TRUE)

contains_generators(w)

cyclic_subword(w, from_i, to_j)

cyclic_conjugates(w)

is_cyclic_conjugate(w, v)

number_syllables(w)

exponent_syllable(w, i)

generator_syllable(w, i)

sub_syllables(w, from_i, to_j)

substituted_word(w, from_i, to_j, by)

is_cyclically_reduced(w)

identity_cyclic_reduction(w)

cyclic_reduction(w, removed = FALSE)

power_of(w, v)

exponent_sum(w, gen)

generator_count(w, gen)

word_index(w, gen)

word_contains(w, gen)
```

## Arguments

- w, v, word:

  Words of the same free group.

- G:

  A free group.

- words:

  A list of words, or for `eliminate_words()` a two column specification
  `list(keys, values)` made with `word_dict()`.

- by:

  The replacement word (default the identity).

- all:

  Replace all occurrences (SymPy's `_all`).

- inverse:

  Also replace occurrences of the inverse of `gen`.

- gen:

  A word (for most functions a generator or its inverse).

- start:

  First position to search.

- from_i, to_j:

  Zero based positions.

- strict:

  When `FALSE`, clip the positions to the word.

- i:

  A zero based syllable position.

- removed:

  Also return the removed conjugating word.

## Value

A word, a list, an integer or a logical.

## Examples

``` r
G <- free_group("x, y, z")
x <- G[[2]]; y <- G[[3]]; z <- G[[4]]
w <- x^5 * y * x^2 * y^-4 * x
fp_str(eliminate_word(w, x, x^2))
#> [1] "x^10*y*x^4*y^-4*x^2"
fp_str(subword(w, 2, 6))
#> [1] "x^3*y"
number_syllables(w)
#> [1] 5
exponent_sum(x^2 * y^3, x)
#> [1] 2
is_cyclic_conjugate(x * y, y * x)
#> [1] TRUE
```
