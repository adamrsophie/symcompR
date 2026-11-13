# Cryptography

``` r

library(symcompR)
```

This vignette covers `sympy.crypto.crypto`: classical ciphers, public
key systems, key exchange, linear feedback shift registers, and Morse
code. The schemes are textbook versions meant for teaching and
experimentation, not for protecting data.

## Classical ciphers

The classical ciphers work on uppercase letters by default, dropping
other characters.
[`AZ()`](https://adamrsophie.github.io/symcompR/reference/crypto.md)
normalizes a string in the same way, and
[`check_and_join()`](https://adamrsophie.github.io/symcompR/reference/crypto.md)
joins and filters text against an alphabet.

``` r

AZ("Hello, world!")
#> [1] "HELLOWORLD"
AZ()
#> [1] "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
check_and_join(c("ab", "c"))
#> [1] "abc"
check_and_join("abcx", "abc", filter = TRUE)
#> [1] "abc"
padded_key("glee", "abcdefghijklmnopqrstuvwxyz")
#> [1] "gleabcdfhijkmnopqrstuvwxyz"
cycle_list(3, 5)
#> [1] 3 4 0 1 2
encipher_shift("GONAVYBEATARMY", 1)
#> [1] "HPOBWZCFBUBSNZ"
decipher_shift("HPOBWZCFBUBSNZ", 1)
#> [1] "GONAVYBEATARMY"
encipher_rot13("Hello World")
#> [1] "URYYBJBEYQ"
decipher_rot13("URYYBJBEYQ")
#> [1] "HELLOWORLD"
encipher_affine("GO NAVY!", c(3, 1))
#> [1] "TROBMV"
decipher_affine("QTW", c(3, 16))
#> [1] "ABC"
encipher_atbash("ABC")
#> [1] "ZYX"
decipher_atbash("ZYX")
#> [1] "ABC"
encipher_substitution("ABC", "BAC", "ABC")
#> [1] "BAC"
ct <- encipher_vigenere("meet me on monday", "encrypt")
ct
#> [1] "QRGKKTHRZQEBPR"
decipher_vigenere(ct, "encrypt")
#> [1] "MEETMEONMONDAY"
```

The Hill cipher multiplies blocks of letters by an invertible matrix
modulo 26:

``` r

key <- Matrix(list(list(1, 2), list(3, 5)))
h <- encipher_hill("meet me on monday", key)
h
#> [1] "UEQDUEODOCTCWQ"
decipher_hill(h, key)
#> [1] "MEETMEONMONDAY"
inv_mod(key, 26)
#> Matrix(2, 2)
#> [21, 2 ]
#> [3 , 25]
```

Textbook RSA with the primes 3 and 11 and public exponent 7:

``` r

pub <- rsa_public_key(3, 11, 7)
priv <- rsa_private_key(3, 11, 7)
unlist(pub)
#> [1] 33  7
unlist(priv)
#> [1] 33  3
decipher_rsa(encipher_rsa(12, pub), priv)
#> [1] 12
```

The rail fence cipher writes the message in a zigzag over several rails:

``` r

encipher_railfence("hello world", 3)
#> [1] "horel ollwd"
decipher_railfence("horel ollwd", 3)
#> [1] "hello world"
```

## Bifid ciphers

The bifid cipher uses a Polybius square. `bifid5` (no J), `bifid6`
(letters and digits), and `bifid10` (printable characters) are the
standard alphabets, and the `*_square()` functions build the key square.

``` r

bifid5
#> [1] "ABCDEFGHIKLMNOPQRSTUVWXYZ"
bifid6
#> [1] "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
nchar(bifid10)
#> [1] 100
bifid5_square("gold bug")
#> Matrix(5, 5)
#> [G, O, L, D, B]
#> [U, A, C, E, F]
#> [H, I, K, M, N]
#> [P, Q, R, S, T]
#> [V, W, X, Y, Z]
bifid6_square("gold bug")
#> Matrix(6, 6)
#> [G, O, L, D, B, U]
#> [A, C, E, F, H, I]
#> [J, K, M, N, P, Q]
#> [R, S, T, V, W, X]
#> [Y, Z, 0, 1, 2, 3]
#> [4, 5, 6, 7, 8, 9]
bifid_square("abcd")
#> Matrix(2, 2)
#> [a, b]
#> [c, d]
encipher_bifid("abc", "b", "abcd")
#> [1] "bdb"
decipher_bifid("bdb", "b", "abcd")
#> [1] "abc"
encipher_bifid5("meet me on friday", "gold bug")
#> [1] "IEILEHFSTSFXEE"
decipher_bifid5("IEILEHFSTSFXEE", "gold bug")
#> [1] "MEETMEONFRIDAY"
encipher_bifid6("AB", "CD")
#> [1] "CP"
decipher_bifid6("CP", "CD")
#> [1] "AB"
```

## Kid RSA

Kid RSA is a simplified RSA for teaching:

``` r

kpub <- kid_rsa_public_key(1, 2, 1, 2)
kpriv <- kid_rsa_private_key(1, 2, 1, 2)
unlist(kpub)
#> [1] 7 2
unlist(kpriv)
#> [1] 7 4
encipher_kid_rsa(3, kpub)
#> [1] 6
decipher_kid_rsa(encipher_kid_rsa(3, kpub), kpriv)
#> [1] 3
```

## ElGamal and Diffie-Hellman

The key generators take a seed so that examples are reproducible.

``` r

priv <- elgamal_private_key(5, seed = list(3))
unlist(priv)
#> [1] 37  2  3
pub <- elgamal_public_key(priv)
unlist(pub)
#> [1] 37  2  8
msg <- encipher_elgamal(36, pub, seed = list(3))
unlist(msg)
#> [1] 8 6
decipher_elgamal(msg, priv)
#> [1] 36
```

Diffie-Hellman key exchange: both parties arrive at the same shared key.

``` r

dh <- list(1031, 14, 565)
unlist(dh_public_key(dh))
#> [1] 1031   14  666
dh_shared_key(dh_public_key(dh), 1000)
#> [1] 505
unlist(dh_private_key(5, seed = list(3)))
#> [1] 37  2  3
```

## Goldwasser-Micali and Blum-Goldwasser

Two probabilistic public key systems. Goldwasser-Micali encrypts one bit
at a time and Blum-Goldwasser uses a Blum Blum Shub generator.

``` r

unlist(gm_private_key(7, 11))
#> [1]  7 11
gpub <- gm_public_key(7, 11, seed = 0)
unlist(gpub)
#> [1] 62 77
c12 <- encipher_gm(12, gpub, seed = 1)
unlist(c12)
#> [1] 54 10 64 23
decipher_gm(c12, list(7, 11))
#> [1] 12
unlist(bg_private_key(23, 31))
#> [1] 23 31
bg_public_key(23, 31)
#> [1] 713
cbg <- encipher_bg(15, 713, seed = 3)
cbg
#> [[1]]
#> [[1]][[1]]
#> [1] 0
#> 
#> [[1]][[2]]
#> [1] 0
#> 
#> [[1]][[3]]
#> [1] 0
#> 
#> [[1]][[4]]
#> [1] 1
#> 
#> 
#> [[2]]
#> [1] 164
decipher_bg(cbg, list(23, 31))
#> [1] 15
```

## Linear feedback shift registers

[`lfsr_sequence()`](https://adamrsophie.github.io/symcompR/reference/crypto.md)
generates the output of a linear feedback shift register over
\mathbb{F}\_p.
[`lfsr_connection_polynomial()`](https://adamrsophie.github.io/symcompR/reference/crypto.md)
recovers the shortest register that produces a sequence (the
Berlekamp-Massey algorithm), and
[`lfsr_autocorrelation()`](https://adamrsophie.github.io/symcompR/reference/crypto.md)
measures how random the sequence looks.

``` r

s <- lfsr_sequence(c(1, 1, 0), c(0, 1, 1), 20, p = 2)
unlist(lapply(s, as.integer))
#>  [1] 0 1 1 1 0 0 1 0 1 1 1 0 0 1 0 1 1 1 0 0
lfsr_connection_polynomial(s)
#> x^3 + x^2 + 1
lfsr_autocorrelation(s, 7, 3)
#> -1/7
```

## Morse code

``` r

m <- encode_morse("RUN SOS")
m
#> [1] ".-.|..-|-.||...|---|..."
decode_morse(m)
#> [1] "RUN SOS"
decode_morse("-.-|.|-.--")
#> [1] "KEY"
```
