# QASM circuits

`Qasm(...)` reads lines of QASM (`qubit`, `h`, `x`, `z`, `s`, `t`,
`cnot`, `swap`, `cphase`, `toffoli`, `c-x`, `c-z`, `measure`, `def`) as
in `sympy.physics.quantum.qasm`. The result is an environment whose
methods are called with `$`, for example `q$get_circuit()` and
`q$get_labels()`. `read_qasm()` reads a string and `read_qasm_file()` a
file.

## Usage

``` r
Qasm(...)

read_qasm(lines)

read_qasm_file(filename)
```

## Arguments

- ...:

  Lines of QASM.

- lines:

  A string of QASM lines.

- filename:

  A file of QASM lines.

## Value

A `Qasm` object.

## Examples

``` r
q <- Qasm("qubit q0", "qubit q1", "h q0", "cnot q0,q1")
q$get_circuit()
#> CNOT(1,0)*H(1)
```
