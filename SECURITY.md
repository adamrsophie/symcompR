# Security Policy

## Supported Versions

Security fixes are considered for the active development branch and the
most recent release.

| Version   | Supported |
|-----------|-----------|
| `0.1.x`   | Yes       |
| `< 0.1.0` | No        |

## Reporting a Vulnerability

Please do not report security vulnerabilities in public GitHub issues.

Report suspected vulnerabilities by email to `adamrsophie@gmail.com`, or
use GitHub private vulnerability reporting if it is enabled for the
repository.

Include as much detail as you can:

- affected version or commit,
- operating system and R version,
- minimal reproduction steps,
- any proof-of-concept code.

## Response Process

The maintainer will acknowledge reports as soon as practical,
investigate the issue, and coordinate a fix if the report is confirmed.
Public disclosure should wait until a fix or mitigation is available,
unless there is an overriding user safety reason to disclose earlier.

## Scope

Security-sensitive areas include:

- parsing of strings:
  [`sympify()`](https://adamrsophie.github.io/symcompR/reference/sympify.md)
  and
  [`parse_expr()`](https://adamrsophie.github.io/symcompR/reference/parse_expr.md)
  parse text and map names through a fixed table instead of evaluating
  it, but they do call any R function supplied through `locals`, and
  [`lambdify()`](https://adamrsophie.github.io/symcompR/reference/lambdify.md)
  returns an R function built from an expression; a way to run code an
  input did not ask for would be a vulnerability,
- arbitrary file writes or path traversal, for example in code
  generation or
  [`pager_print()`](https://adamrsophie.github.io/symcompR/reference/pretty.md),
- unbounded memory or time on small inputs, for example in factorization
  or exact arithmetic on huge exponents,
- exposure of private data in examples, tests, or logs.

A mathematically wrong answer that does not create a software-security
vulnerability should be opened as a normal GitHub issue instead.
