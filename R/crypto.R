# Classical ciphers and textbook public key schemes, after
# sympy/crypto/crypto.py. Seeds follow sympy.core.random: an integer seeds
# Python's Mersenne Twister (R/utilities-pyrandom.R), so seeded results
# match SymPy's; a vector supplies the "random" values in order.

#' Cryptography
#'
#' Ports of SymPy's `sympy.crypto` functions: shift, affine, Atbash,
#' substitution, Vigenere, Hill, bifid (5x5, 6x6 and general squares) and
#' rail fence ciphers; RSA (multiprime, Euler or Carmichael totient), Kid
#' RSA, ElGamal, Diffie-Hellman, Goldwasser-Micali and Blum-Goldwasser
#' schemes; Morse code; linear feedback shift registers. Integers that do not
#' fit a double exactly are returned as `gmp::bigz`.
#'
#' `seed` is `NULL` (draw from the generator), an integer (seed it as SymPy
#' does) or a vector of values to use in order.
#'
#' @param msg,message,ciphertext,phrase A string (or for public key schemes
#'   an integer or the cipher output).
#' @param key A key: an integer, a pair, a string or a matrix.
#' @param symbols The alphabet, a string or character vector.
#' @param .inverse Internal: decipher instead of encipher.
#' @param i An integer message.
#' @param s A string or vector of strings; for LFSRs a sequence.
#' @param ... Primes and exponent (`rsa_public_key(p, q, e)`).
#' @param seed See Details.
#' @return A string, integer, list or expression.
#' @examples
#' encipher_shift("ABC GO TEAM", 3)
#' encipher_vigenere("meet me on monday", "encrypt")
#' encode_morse("SOS")
#' rsa_public_key(3, 5, 7)
#' @name crypto
NULL

.chars <- function(s) if (!length(s) || !nzchar(s)) character() else strsplit(s, "", fixed = TRUE)[[1L]]
.join <- function(x) paste(x, collapse = "")
.uniq_chars <- function(s) .join(unique(.chars(s)))
.int_of <- function(v) if (is_sym(v)) gmp::as.bigz(as.character(v$value)) else gmp::as.bigz(v)

# sympy.utilities.misc.translate: c deletes characters, a dict replaces
# substrings, then a -> b maps single characters
.py_translate <- function(s, a = NULL, b = NULL, c = NULL) {
  mr <- list()
  if (is.null(a)) {
    if (!is.null(c)) stop("c should be None when a=None is passed, instead got ", c, call. = FALSE)
    if (is.null(b)) return(s)
    c <- b; a <- ""; b <- ""
  } else if (is.list(a)) {
    short <- Filter(function(k) nchar(k) == 1L && nchar(a[[k]]) == 1L, names(a))
    mr <- a[setdiff(names(a), short)]
    c <- b
    if (length(short)) { a2 <- .join(short); b <- .join(unlist(a[short])); a <- a2 } else { a <- ""; b <- "" }
  } else if (nchar(a) != nchar(b)) {
    stop("oldchars and newchars have different lengths", call. = FALSE)
  }
  ch <- .chars(s)
  if (!is.null(c) && nzchar(c)) ch <- ch[!(ch %in% .chars(c))]
  s <- .join(ch)
  if (length(mr)) {
    # replace(): longest keys first
    ks <- names(mr)[order(-nchar(names(mr)), names(mr))]
    for (k in ks) s <- gsub(k, mr[[k]], s, fixed = TRUE)
  }
  if (!nzchar(a)) return(s)
  from <- .chars(a); to <- .chars(b)
  # str.maketrans: a later duplicate wins
  map <- stats::setNames(to, from)
  map <- map[!duplicated(from, fromLast = TRUE)]
  ch <- .chars(s)
  hit <- ch %in% names(map)
  ch[hit] <- map[ch[hit]]
  .join(ch)
}

#' @rdname crypto
#' @export
AZ <- function(s = NULL) {
  if (is.null(s) || (length(s) == 1L && !nzchar(s))) return(.join(LETTERS))
  rv <- vapply(s, function(i) check_and_join(strsplit(toupper(i), "[[:space:]]+")[[1L]], .join(LETTERS), filter = TRUE), "")
  unname(rv)
}

#' @rdname crypto
#' @export
bifid5 <- .join(setdiff(LETTERS, "J"))
#' @rdname crypto
#' @export
bifid6 <- .join(c(LETTERS, 0:9))
#' @rdname crypto
#' @export
bifid10 <- .join(c(0:9, letters, LETTERS, .chars("!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~"), " ", "\t", "\n", "\r", "\v", "\f"))

#' @rdname crypto
#' @export
padded_key <- function(key, symbols) {
  syms <- unique(.chars(symbols))
  if (length(syms) != nchar(symbols)) {
    sc <- .chars(symbols)
    stop("duplicate characters in symbols: ", .join(sort(unique(sc[duplicated(sc)]), method = "radix")), call. = FALSE)
  }
  extra <- setdiff(.chars(key), syms)
  if (length(extra)) stop("characters in key but not symbols: ", .join(sort(extra, method = "radix")), call. = FALSE)
  key0 <- .uniq_chars(key)
  paste0(key0, .py_translate(.join(syms), NULL, key0))
}

#' @rdname crypto
#' @param filter Drop characters not in `symbols` instead of failing.
#' @export
check_and_join <- function(phrase, symbols = NULL, filter = NULL) {
  rv <- .join(phrase)
  if (!is.null(symbols)) {
    symbols <- check_and_join(symbols)
    missing <- .join(sort(setdiff(unique(.chars(rv)), .chars(symbols)), method = "radix"))
    if (nzchar(missing)) {
      if (!isTRUE(filter)) stop("characters in phrase but not symbols: \"", missing, "\"", call. = FALSE)
      rv <- .py_translate(rv, NULL, missing)
    }
  }
  rv
}

.prep <- function(msg, key, alp, default = NULL) {
  if (is.null(alp) || !length(alp) || !nzchar(.join(alp))) {
    if (is.null(default)) {
      alp <- AZ(); msg <- AZ(msg); key <- AZ(key)
    } else alp <- default
  } else alp <- .join(alp)
  key <- check_and_join(key, alp, filter = TRUE)
  msg <- check_and_join(msg, alp, filter = TRUE)
  list(msg, key, alp)
}

#' @rdname crypto
#' @param k,n Integers.
#' @export
cycle_list <- function(k, n) {
  k <- k %% n
  c(if (k < n) k:(n - 1L), if (k > 0) 0:(k - 1L))
}

#' @rdname crypto
#' @export
encipher_shift <- function(msg, key, symbols = NULL) {
  p <- .prep(msg, "", symbols); msg <- p[[1L]]; A <- p[[3L]]
  N <- nchar(A)
  shift <- N - key %% N
  k <- paste0(substring(A, shift + 1L), substr(A, 1L, shift))
  .py_translate(msg, k, A)
}
#' @rdname crypto
#' @export
decipher_shift <- function(msg, key, symbols = NULL) encipher_shift(msg, -key, symbols)
#' @rdname crypto
#' @export
encipher_rot13 <- function(msg, symbols = NULL) encipher_shift(msg, 13, symbols)
#' @rdname crypto
#' @export
decipher_rot13 <- function(msg, symbols = NULL) decipher_shift(msg, 13, symbols)

#' @rdname crypto
#' @export
encipher_affine <- function(msg, key, symbols = NULL, .inverse = FALSE) {
  p <- .prep(msg, "", symbols); msg <- p[[1L]]; A <- p[[3L]]
  N <- nchar(A)
  a <- key[[1L]]; b <- key[[2L]]
  if (as.numeric(gmp::gcd.bigz(a, N)) != 1) stop("the key multiplier must be coprime to the alphabet size", call. = FALSE)
  if (.inverse) {
    cc <- as.numeric(gmp::inv.bigz(a, N)); d <- -b * cc
    a <- cc; b <- d
  }
  B <- .join(.chars(A)[((a * (0:(N - 1L)) + b) %% N) + 1L])
  .py_translate(msg, A, B)
}
#' @rdname crypto
#' @export
decipher_affine <- function(msg, key, symbols = NULL) encipher_affine(msg, key, symbols, .inverse = TRUE)
#' @rdname crypto
#' @export
encipher_atbash <- function(msg, symbols = NULL) encipher_affine(msg, c(25, 25), symbols)
#' @rdname crypto
#' @export
decipher_atbash <- function(msg, symbols = NULL) decipher_affine(msg, c(25, 25), symbols)

#' @rdname crypto
#' @param old,new Characters to replace and their replacements; `old` may be
#'   a named list of substitutions.
#' @export
encipher_substitution <- function(msg, old, new = NULL) .py_translate(msg, old, new)

#' @rdname crypto
#' @export
encipher_vigenere <- function(msg, key, symbols = NULL) {
  p <- .prep(msg, key, symbols); m <- .chars(p[[1L]]); A <- .chars(p[[3L]])
  K <- match(.chars(p[[2L]]), A) - 1L
  M <- match(m, A) - 1L
  .join(A[((M + K[((seq_along(M) - 1L) %% length(K)) + 1L]) %% length(A)) + 1L])
}
#' @rdname crypto
#' @export
decipher_vigenere <- function(msg, key, symbols = NULL) {
  p <- .prep(msg, key, symbols); A <- .chars(p[[3L]])
  K <- match(.chars(p[[2L]]), A) - 1L
  C <- match(.chars(p[[1L]]), A) - 1L
  .join(A[((-K[((seq_along(C) - 1L) %% length(K)) + 1L] + C) %% length(A)) + 1L])
}

#' @rdname crypto
#' @param M A square matrix.
#' @param m A modulus.
#' @export
inv_mod <- function(M, m) {
  M <- Matrix(M)
  if (M$nrow != M$ncol) stop("Only square matrices have inverses.", call. = FALSE)
  N <- M$nrow
  dk <- .int_of(det(M))
  di <- tryCatch(gmp::inv.bigz(dk, m), error = function(e) NA)
  if (is.na(di)) stop("Matrix is not invertible (mod ", m, ")", call. = FALSE)
  adj <- adjugate(M)
  Matrix(N, N, f = function(i, j) as.character((di * .int_of(.mget(adj, i, j))) %% m))
}

.hill_apply <- function(key, P, k, N, A) {
  m <- length(P) %/% k
  out <- character()
  for (j in seq_len(m)) {
    v <- Matrix(lapply(P[(k * (j - 1L) + 1L):(k * j)], as_sym))
    r <- key %*% v
    out <- c(out, vapply(r$data, function(c) A[as.numeric(.int_of(c) %% N) + 1L], ""))
  }
  .join(out)
}

#' @rdname crypto
#' @param pad The padding character.
#' @export
encipher_hill <- function(msg, key, symbols = NULL, pad = "Q") {
  key <- Matrix(key)
  if (key$nrow != key$ncol) stop("key must be a square matrix", call. = FALSE)
  if (nchar(pad) != 1L) stop("pad must be a single character", call. = FALSE)
  p <- .prep(msg, pad, symbols); A <- .chars(p[[3L]])
  P <- match(.chars(p[[1L]]), A) - 1L
  k <- key$ncol; r <- length(P) %% k
  if (r) P <- c(P, rep(match(p[[2L]], A) - 1L, k - r))
  .hill_apply(key, P, k, length(A), A)
}
#' @rdname crypto
#' @export
decipher_hill <- function(msg, key, symbols = NULL) {
  key <- Matrix(key)
  if (key$nrow != key$ncol) stop("key must be a square matrix", call. = FALSE)
  p <- .prep(msg, "", symbols); A <- .chars(p[[3L]])
  C <- match(.chars(p[[1L]]), A) - 1L
  k <- key$ncol; r <- length(C) %% k
  if (r) C <- c(C, rep(0L, k - r))
  .hill_apply(inv_mod(key, length(A)), C, k, length(A), A)
}

.bifid_long_key <- function(key, A) {
  N <- sqrt(nchar(A))
  if (N != floor(N)) stop("Length of alphabet (", nchar(A), ") is not a square number.", call. = FALSE)
  lk <- .chars(.uniq_chars(key))
  if (!length(lk)) lk <- .chars(A)
  if (length(lk) < N^2) lk <- c(lk, setdiff(.chars(A), lk))
  list(lk = lk, N = as.integer(N))
}

#' @rdname crypto
#' @export
encipher_bifid <- function(msg, key, symbols = NULL) {
  p <- .prep(msg, key, symbols, bifid10)
  b <- .bifid_long_key(p[[2L]], p[[3L]]); lk <- b$lk; N <- b$N
  pos <- match(.chars(p[[1L]]), lk) - 1L
  rc <- c(pos %/% N, pos %% N)
  rr <- rc[c(TRUE, FALSE)]; cc <- rc[c(FALSE, TRUE)]
  .join(lk[rr * N + cc + 1L])
}
#' @rdname crypto
#' @export
decipher_bifid <- function(msg, key, symbols = NULL) {
  p <- .prep(msg, "", symbols, bifid10)
  b <- .bifid_long_key(key, p[[3L]]); lk <- b$lk; N <- b$N
  pos <- match(.chars(p[[1L]]), lk) - 1L
  rc <- as.vector(rbind(pos %/% N, pos %% N))
  n <- length(pos)
  .join(lk[rc[seq_len(n)] * N + rc[n + seq_len(n)] + 1L])
}

#' @rdname crypto
#' @export
bifid_square <- function(key) {
  A <- .chars(.uniq_chars(.join(key)))
  n <- sqrt(length(A))
  if (n != floor(n)) stop("Length of alphabet (", length(A), ") is not a square number.", call. = FALSE)
  Matrix(n, n, f = function(i, j) Symbol(A[n * (i - 1L) + j]))
}

.bifid_n <- function(msg, key, alpha, fn) {
  p <- .prep(toupper(msg), toupper(key), NULL, alpha)
  fn(p[[1L]], "", padded_key(p[[2L]], alpha))
}
#' @rdname crypto
#' @export
encipher_bifid5 <- function(msg, key) .bifid_n(msg, key, bifid5, encipher_bifid)
#' @rdname crypto
#' @export
decipher_bifid5 <- function(msg, key) .bifid_n(msg, key, bifid5, decipher_bifid)
#' @rdname crypto
#' @export
encipher_bifid6 <- function(msg, key) .bifid_n(msg, key, bifid6, encipher_bifid)
#' @rdname crypto
#' @export
decipher_bifid6 <- function(msg, key) .bifid_n(msg, key, bifid6, decipher_bifid)
.bifid_sq <- function(key, alpha) {
  if (is.null(key) || !nzchar(key)) return(bifid_square(alpha))
  bifid_square(padded_key(.prep("", toupper(key), NULL, alpha)[[2L]], alpha))
}
#' @rdname crypto
#' @export
bifid5_square <- function(key = NULL) .bifid_sq(key, bifid5)
#' @rdname crypto
#' @export
bifid6_square <- function(key = NULL) .bifid_sq(key, bifid6)

# ---------------------------------------------------------------------------
# RSA

.rsa_key <- function(args, public = TRUE, private = TRUE, totient = "Euler", index = NULL, multipower = NULL) {
  if (length(args) < 2L) return(FALSE)
  if (!(totient %in% c("Euler", "Carmichael")))
    stop("The argument totient=", totient, " should either be 'Euler', 'Carmichalel'.", call. = FALSE)
  if (!is.null(index) && totient != "Carmichael")
    stop("Setting the 'index' keyword argument requires totientnotation to be specified as 'Carmichael'.", call. = FALSE)
  primes <- lapply(args[-length(args)], gmp::as.bigz); e <- gmp::as.bigz(args[[length(args)]])
  if (!all(vapply(primes, function(p) isTRUE(isprime(p)), NA))) {
    primes <- unlist(lapply(primes, function(p) as.character(factorint(p, multiple = TRUE))))
    primes <- lapply(primes, gmp::as.bigz)
  }
  n <- Reduce(`*`, primes)
  tally <- table(vapply(primes, as.character, ""))
  if (!all(tally == 1L) && !isTRUE(multipower))
    warning("Non-distinctive primes found in the factors ", paste(vapply(primes, as.character, ""), collapse = ", "),
            ". The cipher may not be decryptable for some numbers in the complete residue system Z[", as.character(n),
            "], but the cipher can still be valid if you restrict the domain to be the reduced residue system Z*[",
            as.character(n), "]. You can pass the flag multipower=True if you want to suppress this warning.", call. = FALSE)
  # totient(tally) and reduced_totient(tally) of the factor multiset
  pk <- names(tally)
  phis <- lapply(pk, function(k) {
    pr <- gmp::as.bigz(k); e <- as.integer(tally[[k]])
    if (totient == "Carmichael" && pr == 2 && e >= 3) return(gmp::as.bigz(2)^(e - 2))
    if (totient == "Carmichael" && pr == 2) return(gmp::as.bigz(2)^(e - 1))
    pr^(e - 1) * (pr - 1)
  })
  phi <- if (totient == "Euler") Reduce(`*`, phis) else Reduce(function(a, b) a * b / gmp::gcd.bigz(a, b), phis)
  phi <- gmp::as.bigz(phi)
  if (gmp::gcd.bigz(e, phi) == 1) {
    if (public && !private) {
      if (!is.null(index)) e <- e %% phi + gmp::as.bigz(index) * phi
      return(list(.as_out(n), .as_out(e)))
    }
    if (private && !public) {
      d <- gmp::inv.bigz(e, phi)
      if (!is.null(index)) d <- d + gmp::as.bigz(index) * phi
      return(list(.as_out(n), .as_out(d)))
    }
  }
  FALSE
}

#' @rdname crypto
#' @param totient `"Euler"` or `"Carmichael"`.
#' @param index Adds `index` times the totient to the exponent.
#' @param multipower Silence the warning about repeated primes.
#' @export
rsa_public_key <- function(..., totient = "Euler", index = NULL, multipower = NULL)
  .rsa_key(list(...), public = TRUE, private = FALSE, totient = totient, index = index, multipower = multipower)
#' @rdname crypto
#' @export
rsa_private_key <- function(..., totient = "Euler", index = NULL, multipower = NULL)
  .rsa_key(list(...), public = FALSE, private = TRUE, totient = totient, index = index, multipower = multipower)

.rsa_apply <- function(i, key, factors = NULL) {
  n <- gmp::as.bigz(key[[1L]]); d <- gmp::as.bigz(key[[2L]])
  if (!length(factors)) return(.as_out(.powm(i, d, n)))
  fs <- lapply(factors, gmp::as.bigz)
  coprime <- TRUE
  for (a in seq_along(fs)) for (b in seq_along(fs)) if (a < b && gmp::gcd.bigz(fs[[a]], fs[[b]]) != 1) coprime <- FALSE
  if (Reduce(`*`, fs) == n && coprime) {
    mods <- vapply(fs, function(p) as.character(.powm(i, d, p)), "")
    r <- crt(vapply(fs, as.character, ""), mods)
    if (is.null(r) || identical(r, FALSE)) stop("CRT failed", call. = FALSE)
    return(.as_out(gmp::as.bigz(r[[1L]])))
  }
  .rsa_apply(i, key, NULL)
}
#' @rdname crypto
#' @param factors The prime factors of the modulus, for decryption by the
#'   Chinese remainder theorem.
#' @export
encipher_rsa <- function(i, key, factors = NULL) .rsa_apply(i, key, factors)
#' @rdname crypto
#' @export
decipher_rsa <- function(i, key, factors = NULL) .rsa_apply(i, key, factors)

.kid <- function(a, b, A, B) {
  M <- a * b - 1; e <- A * M + a; d <- B * M + b
  list(n = (e * d - 1) %/% M, e = e, d = d)
}
#' @rdname crypto
#' @param a,b,A,B Integers.
#' @export
kid_rsa_public_key <- function(a, b, A, B) { k <- .kid(a, b, A, B); list(k$n, k$e) }
#' @rdname crypto
#' @export
kid_rsa_private_key <- function(a, b, A, B) { k <- .kid(a, b, A, B); list(k$n, k$d) }
#' @rdname crypto
#' @export
encipher_kid_rsa <- function(msg, key) (msg * key[[2L]]) %% key[[1L]]
#' @rdname crypto
#' @export
decipher_kid_rsa <- function(msg, key) (msg * key[[2L]]) %% key[[1L]]

# ---------------------------------------------------------------------------
# Morse code

.morse_char <- c(".-" = "A", "-..." = "B", "-.-." = "C", "-.." = "D", "." = "E", "..-." = "F",
  "--." = "G", "...." = "H", ".." = "I", ".---" = "J", "-.-" = "K", ".-.." = "L", "--" = "M",
  "-." = "N", "---" = "O", ".--." = "P", "--.-" = "Q", ".-." = "R", "..." = "S", "-" = "T",
  "..-" = "U", "...-" = "V", ".--" = "W", "-..-" = "X", "-.--" = "Y", "--.." = "Z",
  "-----" = "0", ".----" = "1", "..---" = "2", "...--" = "3", "....-" = "4", "....." = "5",
  "-...." = "6", "--..." = "7", "---.." = "8", "----." = "9", ".-.-.-" = ".", "--..--" = ",",
  "---..." = ":", "-.-.-." = ";", "..--.." = "?", "-....-" = "-", "..--.-" = "_", "-.--." = "(",
  "-.--.-" = ")", ".----." = "'", "-...-" = "=", ".-.-." = "+", "-..-." = "/", ".--.-." = "@",
  "...-..-" = "$", "-.-.--" = "!")
.char_morse <- stats::setNames(names(.morse_char), unname(.morse_char))

#' @rdname crypto
#' @param sep The letter separator.
#' @param mapping A named character vector from characters to codes (or
#'   codes to characters for decoding).
#' @export
encode_morse <- function(msg, sep = "|", mapping = NULL) {
  mapping <- if (is.null(mapping)) .char_morse else mapping
  if (sep %in% names(mapping)) stop("Separator '", sep, "' is already used in the mapping.", call. = FALSE)
  word_sep <- strrep(sep, 2L)
  mapping[" "] <- word_sep
  ch <- .chars(msg)
  suffix <- length(ch) && grepl("[[:space:]]", ch[length(ch)])
  words <- strsplit(trimws(msg), "[[:space:]]+")[[1L]]
  msg <- paste(words, collapse = " ")
  chars <- unique(.chars(gsub("[[:space:]]", "", msg)))
  msg <- .py_translate(msg, NULL, .join(setdiff(chars, names(mapping))))
  words <- Filter(nzchar, strsplit(msg, "[[:space:]]+")[[1L]])
  out <- vapply(words, function(w) paste(unname(mapping[.chars(w)]), collapse = sep), "")
  paste0(paste(out, collapse = word_sep), if (suffix) word_sep else "")
}

#' @rdname crypto
#' @export
decode_morse <- function(msg, sep = "|", mapping = NULL) {
  mapping <- if (is.null(mapping)) .morse_char else mapping
  word_sep <- strrep(sep, 2L)
  # str.strip(word_sep) strips any of its characters from both ends
  esc <- gsub("([][{}()+*^$|\\\\?.-])", "\\\\\\1", sep)
  msg <- gsub(paste0("^(", esc, ")+|(", esc, ")+$"), "", msg)
  words <- strsplit(msg, word_sep, fixed = TRUE)[[1L]]
  paste(vapply(words, function(w) .join(unname(mapping[strsplit(w, sep, fixed = TRUE)[[1L]]])), ""), collapse = " ")
}

# ---------------------------------------------------------------------------
# Linear feedback shift registers

#' @rdname crypto
#' @param fill The initial fill.
#' @param p The modulus of the field elements (SymPy reads it from the `FF`
#'   elements of `key`).
#' @export
lfsr_sequence <- function(key, fill, n, p = attr(key, "modulus")) {
  if (is.null(p)) p <- 2
  key <- as.numeric(key) %% p; s <- as.numeric(fill) %% p
  k <- length(s)
  L <- numeric(n)
  for (i in seq_len(n)) {
    s0 <- s
    L[i] <- s[1L]
    s <- s[-1L]
    s <- c(s, sum(key * s0) %% p)
  }
  structure(L, modulus = p)
}

#' @rdname crypto
#' @param L A sequence.
#' @param P The period.
#' @export
lfsr_autocorrelation <- function(L, P, k) {
  P <- as.integer(P); k <- as.integer(k)
  L0 <- as.integer(L)[seq_len(P)]
  L1 <- c(L0, L0[seq_len(k)])
  tot <- sum((-1)^(L1[seq_len(P)] + L1[seq_len(P) + k]))
  Rational(tot, P)
}

#' @rdname crypto
#' @export
lfsr_connection_polynomial <- function(s, p = attr(s, "modulus")) {
  if (is.null(p)) p <- 2
  s <- as.numeric(s)
  x <- Symbol("x")
  # coefficient vectors, index 1 for x^0
  C <- 1; B <- 1; m <- 1; b <- 1; L <- 0; N <- 0
  shift_sub <- function(C, d, b, m, B) {
    f <- (d * as.numeric(.powm(b, p - 2, p))) %% p
    out <- numeric(max(length(C), length(B) + m)); out[seq_along(C)] <- C
    idx <- seq_along(B) + m
    out[idx] <- out[idx] - f * B
    out
  }
  while (N < length(s)) {
    if (L > 0) {
      dC <- max(which(C != 0)) - 1
      r <- min(L + 1, dC + 1)
      # s[N - i] wraps around for a negative index, as Python's indexing
      d <- (s[N + 1] + sum(vapply(seq_len(r - 1), function(i) C[i + 1] * s[((N - i) %% length(s)) + 1], 0))) %% p
    }
    if (L == 0) d <- s[N + 1]
    if (d == 0) { m <- m + 1; N <- N + 1 }
    if (d > 0) {
      if (2 * L > N) {
        C <- shift_sub(C, d, b, m, B); m <- m + 1; N <- N + 1
      } else {
        Tm <- C
        C <- shift_sub(C, d, b, m, B)
        L <- N + 1 - L; m <- 1; b <- d; B <- Tm; N <- N + 1
      }
    }
  }
  dC <- max(which(C != 0)) - 1
  .add(lapply(0:dC, function(i) .mul(list(.int(C[i + 1] %% p), .pow(x, .int(i))))))
}

# ---------------------------------------------------------------------------
# ElGamal, Diffie-Hellman, Goldwasser-Micali, Blum-Goldwasser, rail fence

#' @rdname crypto
#' @param digit The key size: the prime is the next one after `2^digit`.
#' @export
elgamal_private_key <- function(digit = 10, seed = NULL) {
  rr <- .sympy_randrange(seed)
  p <- .int_of(nextprime(gmp::as.bigz(2)^digit))
  list(.as_out(p), .as_out(.int_of(primitive_root(p))), rr(2, p))
}
#' @rdname crypto
#' @export
elgamal_public_key <- function(key) list(key[[1L]], key[[2L]], .as_out(.powm(key[[2L]], key[[3L]], key[[1L]])))
#' @rdname crypto
#' @export
encipher_elgamal <- function(i, key, seed = NULL) {
  p <- key[[1L]]; r <- key[[2L]]; e <- key[[3L]]
  if (i < 0 || i >= p) stop("Message (", i, ") should be in range(", p, ")", call. = FALSE)
  a <- .sympy_randrange(seed)(2, p)
  list(.as_out(.powm(r, a, p)), .as_out((gmp::as.bigz(i) * .powm(e, a, p)) %% gmp::as.bigz(p)))
}
#' @rdname crypto
#' @export
decipher_elgamal <- function(msg, key) {
  p <- key[[1L]]; d <- key[[3L]]
  u <- .powm(msg[[1L]], -gmp::as.bigz(d), p)
  .as_out((u * gmp::as.bigz(msg[[2L]])) %% gmp::as.bigz(p))
}

#' @rdname crypto
#' @export
dh_private_key <- function(digit = 10, seed = NULL) {
  p <- .int_of(nextprime(gmp::as.bigz(2)^digit))
  g <- .int_of(primitive_root(p))
  a <- .sympy_randrange(seed)(2, p)
  list(.as_out(p), .as_out(g), a)
}
#' @rdname crypto
#' @export
dh_public_key <- function(key) list(key[[1L]], key[[2L]], .as_out(.powm(key[[2L]], key[[3L]], key[[1L]])))
#' @rdname crypto
#' @export
dh_shared_key <- function(key, b) {
  p <- key[[1L]]
  if (1 >= b || b >= p) stop("Value of b should be greater 1 and less than prime ", p, ".", call. = FALSE)
  .as_out(.powm(key[[3L]], b, p))
}

.legendre2 <- function(a, p) {
  sig <- .powm(a, (gmp::as.bigz(p) - 1) %/% 2, p)
  if (sig == 1) 1L else if (sig == 0) 0L else -1L
}

#' @rdname crypto
#' @param q A prime.
#' @export
gm_private_key <- function(p, q, a = NULL) {
  if (p == q) stop("expected distinct primes, got two copies of ", p, call. = FALSE)
  if (!isTRUE(isprime(p)) || !isTRUE(isprime(q))) stop("first two arguments must be prime, got ", p, " of ", q, call. = FALSE)
  if (p == 2 || q == 2) stop("first two arguments must not be even, got ", p, " of ", q, call. = FALSE)
  list(p, q)
}
#' @rdname crypto
#' @export
gm_public_key <- function(p, q, a = NULL, seed = NULL) {
  gm_private_key(p, q)
  N <- p * q
  if (is.null(a)) {
    rr <- .sympy_randrange(seed)
    repeat {
      a <- rr(N)
      if (.legendre2(a, p) == -1L && .legendre2(a, q) == -1L) break
    }
  } else if (.legendre2(a, p) != -1L || .legendre2(a, q) != -1L) return(FALSE)
  list(a, N)
}
#' @rdname crypto
#' @export
encipher_gm <- function(i, key, seed = NULL) {
  if (i < 0) stop("message must be a non-negative integer: got ", i, " instead", call. = FALSE)
  a <- key[[1L]]; N <- key[[2L]]
  bits <- integer()
  i <- gmp::as.bigz(i)
  while (i > 0) { bits <- c(bits, as.integer(i %% 2)); i <- i %/% 2 }
  rr <- NULL
  nxt <- function() {
    if (is.null(rr)) rr <<- .sympy_randrange(seed)
    repeat { y <- rr(N); if (gmp::gcd.bigz(y, N) == 1) return(y) }
  }
  lapply(rev(bits), function(b) .as_out((gmp::as.bigz(nxt())^2 * gmp::as.bigz(a)^b) %% gmp::as.bigz(N)))
}
#' @rdname crypto
#' @export
decipher_gm <- function(message, key) {
  p <- key[[1L]]; q <- key[[2L]]
  res <- function(m, p) .legendre2(m, p) > 0
  m <- gmp::as.bigz(0)
  for (c in message) m <- m * 2 + as.integer(!(res(c, p) && res(c, q)))
  .as_out(m)
}

.rail_order <- function(n, rails) {
  r <- 0:(rails - 1L)
  pat <- c(r, if (rails > 2L) rev(r)[-c(1L, rails)])
  if (rails == 1L) pat <- 0L
  order(rep_len(pat, n), method = "radix")
}
#' @rdname crypto
#' @param rails The number of rails.
#' @export
encipher_railfence <- function(message, rails) .join(.chars(message)[.rail_order(nchar(message), rails)])
#' @rdname crypto
#' @export
decipher_railfence <- function(ciphertext, rails) {
  ch <- .chars(ciphertext)
  res <- character(length(ch))
  res[.rail_order(length(ch), rails)] <- ch
  .join(res)
}

#' @rdname crypto
#' @export
bg_private_key <- function(p, q) {
  if (!isTRUE(isprime(p)) || !isTRUE(isprime(q))) stop("the two arguments must be prime, got ", p, " and ", q, call. = FALSE)
  if (p == q) stop("the two arguments must be distinct, got two copies of ", p, ". ", call. = FALSE)
  if ((p - 3) %% 4 != 0 || (q - 3) %% 4 != 0)
    stop("the two arguments must be congruent to 3 mod 4, got ", p, " and ", q, call. = FALSE)
  list(p, q)
}
#' @rdname crypto
#' @export
bg_public_key <- function(p, q) { bg_private_key(p, q); p * q }
#' @rdname crypto
#' @export
encipher_bg <- function(i, key, seed = NULL) {
  if (i < 0) stop("message must be a non-negative integer: got ", i, " instead", call. = FALSE)
  bits <- integer(); i <- gmp::as.bigz(i)
  while (i > 0) { bits <- c(bits, as.integer(i %% 2)); i <- i %/% 2 }
  bits <- rev(bits); L <- length(bits)
  key <- gmp::as.bigz(key)
  r <- .sympy_randint(seed)(2, as.numeric(key - 1))
  x <- (gmp::as.bigz(r)^2) %% key
  xL <- .powm(x, gmp::as.bigz(2)^L, key)
  rb <- integer(L)
  for (k in seq_len(L)) { rb[k] <- as.integer(x %% 2); x <- (x^2) %% key }
  list(as.list(bitwXor(bits, rb)), .as_out(xL))
}
#' @rdname crypto
#' @export
decipher_bg <- function(message, key) {
  p <- key[[1L]]; q <- key[[2L]]
  em <- unlist(message[[1L]]); y <- gmp::as.bigz(message[[2L]])
  N <- gmp::as.bigz(p) * gmp::as.bigz(q)
  L <- length(em)
  # SymPy computes these exponents in floating point
  pt <- gmp::as.bigz(trunc(((p + 1) / 4)^L)); qt <- gmp::as.bigz(trunc(((q + 1) / 4)^L))
  rp <- .powm(y, pt, p); rq <- .powm(y, qt, q)
  x <- (q * gmp::inv.bigz(q, p) * rp + p * gmp::inv.bigz(p, q) * rq) %% N
  ob <- integer(L)
  for (k in seq_len(L)) { ob[k] <- as.integer(x %% 2); x <- (x^2) %% N }
  m <- gmp::as.bigz(0)
  for (k in seq_len(L)) m <- m * 2 + bitwXor(em[k], ob[k])
  .as_out(m)
}
