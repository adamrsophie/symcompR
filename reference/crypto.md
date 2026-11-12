# Cryptography

Ports of SymPy's `sympy.crypto` functions: shift, affine, Atbash,
substitution, Vigenere, Hill, bifid (5x5, 6x6 and general squares) and
rail fence ciphers; RSA (multiprime, Euler or Carmichael totient), Kid
RSA, ElGamal, Diffie-Hellman, Goldwasser-Micali and Blum-Goldwasser
schemes; Morse code; linear feedback shift registers. Integers that do
not fit a double exactly are returned as
[`gmp::bigz`](https://rdrr.io/pkg/gmp/man/biginteger.html).

## Usage

``` r
AZ(s = NULL)

bifid5

bifid6

bifid10

padded_key(key, symbols)

check_and_join(phrase, symbols = NULL, filter = NULL)

cycle_list(k, n)

encipher_shift(msg, key, symbols = NULL)

decipher_shift(msg, key, symbols = NULL)

encipher_rot13(msg, symbols = NULL)

decipher_rot13(msg, symbols = NULL)

encipher_affine(msg, key, symbols = NULL, .inverse = FALSE)

decipher_affine(msg, key, symbols = NULL)

encipher_atbash(msg, symbols = NULL)

decipher_atbash(msg, symbols = NULL)

encipher_substitution(msg, old, new = NULL)

encipher_vigenere(msg, key, symbols = NULL)

decipher_vigenere(msg, key, symbols = NULL)

inv_mod(M, m)

encipher_hill(msg, key, symbols = NULL, pad = "Q")

decipher_hill(msg, key, symbols = NULL)

encipher_bifid(msg, key, symbols = NULL)

decipher_bifid(msg, key, symbols = NULL)

bifid_square(key)

encipher_bifid5(msg, key)

decipher_bifid5(msg, key)

encipher_bifid6(msg, key)

decipher_bifid6(msg, key)

bifid5_square(key = NULL)

bifid6_square(key = NULL)

rsa_public_key(..., totient = "Euler", index = NULL, multipower = NULL)

rsa_private_key(..., totient = "Euler", index = NULL, multipower = NULL)

encipher_rsa(i, key, factors = NULL)

decipher_rsa(i, key, factors = NULL)

kid_rsa_public_key(a, b, A, B)

kid_rsa_private_key(a, b, A, B)

encipher_kid_rsa(msg, key)

decipher_kid_rsa(msg, key)

encode_morse(msg, sep = "|", mapping = NULL)

decode_morse(msg, sep = "|", mapping = NULL)

lfsr_sequence(key, fill, n, p = attr(key, "modulus"))

lfsr_autocorrelation(L, P, k)

lfsr_connection_polynomial(s, p = attr(s, "modulus"))

elgamal_private_key(digit = 10, seed = NULL)

elgamal_public_key(key)

encipher_elgamal(i, key, seed = NULL)

decipher_elgamal(msg, key)

dh_private_key(digit = 10, seed = NULL)

dh_public_key(key)

dh_shared_key(key, b)

gm_private_key(p, q, a = NULL)

gm_public_key(p, q, a = NULL, seed = NULL)

encipher_gm(i, key, seed = NULL)

decipher_gm(message, key)

encipher_railfence(message, rails)

decipher_railfence(ciphertext, rails)

bg_private_key(p, q)

bg_public_key(p, q)

encipher_bg(i, key, seed = NULL)

decipher_bg(message, key)
```

## Arguments

- s:

  A string or vector of strings; for LFSRs a sequence.

- key:

  A key: an integer, a pair, a string or a matrix.

- symbols:

  The alphabet, a string or character vector.

- filter:

  Drop characters not in `symbols` instead of failing.

- k, n:

  Integers.

- msg, message, ciphertext, phrase:

  A string (or for public key schemes an integer or the cipher output).

- .inverse:

  Internal: decipher instead of encipher.

- old, new:

  Characters to replace and their replacements; `old` may be a named
  list of substitutions.

- M:

  A square matrix.

- m:

  A modulus.

- pad:

  The padding character.

- ...:

  Primes and exponent (`rsa_public_key(p, q, e)`).

- totient:

  `"Euler"` or `"Carmichael"`.

- index:

  Adds `index` times the totient to the exponent.

- multipower:

  Silence the warning about repeated primes.

- i:

  An integer message.

- factors:

  The prime factors of the modulus, for decryption by the Chinese
  remainder theorem.

- a, b, A, B:

  Integers.

- sep:

  The letter separator.

- mapping:

  A named character vector from characters to codes (or codes to
  characters for decoding).

- fill:

  The initial fill.

- p:

  The modulus of the field elements (SymPy reads it from the `FF`
  elements of `key`).

- L:

  A sequence.

- P:

  The period.

- digit:

  The key size: the prime is the next one after `2^digit`.

- seed:

  See Details.

- q:

  A prime.

- rails:

  The number of rails.

## Value

A string, integer, list or expression.

## Details

`seed` is `NULL` (draw from the generator), an integer (seed it as SymPy
does) or a vector of values to use in order.

## Examples

``` r
encipher_shift("ABC GO TEAM", 3)
#> [1] "DEFJRWHDP"
encipher_vigenere("meet me on monday", "encrypt")
#> [1] "QRGKKTHRZQEBPR"
encode_morse("SOS")
#> [1] "...|---|..."
rsa_public_key(3, 5, 7)
#> [[1]]
#> [1] 15
#> 
#> [[2]]
#> [1] 7
#> 
```
