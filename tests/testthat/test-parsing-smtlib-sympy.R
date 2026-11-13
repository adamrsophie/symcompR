skip_unless_full_tests()

# Expected values are SymPy's outputs for the same calls (sympy/parsing/tests/test_smtlib.py);
# the rows come from the battery that compares this port with SymPy.

test_that("smtlib: test assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("\n    (set-option :print-success false)\n    (set-logic QF_UF)\n    (declare-sort x5 0)\n    (declare-fun x4 () x5)\n    (declare-fun x6 () x5)\n    (declare-fun x7 () x5)\n    (declare-fun x1 () x5)\n    (declare-fun x3 () x5)\n    (assert (= x6 x7))\n    (assert (not (= x4 x3)))\n    (assert (not (or (= x4 x1) (= x1 x3))))\n    (check-sat)\n    (exit)\n    ")$assertions })), "[Eq(x6, x7), Ne(x4, x3), ~(Eq(x1, x3) | Eq(x4, x1))]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (set-option :print-success false)\n    (set-logic QF_UF)\n    (declare-sort x0 0)\n    (declare-fun x2 () x0)\n    (declare-fun x4 () x0)\n    (declare-fun x7 () x0)\n    (declare-fun x8 () x0)\n    (declare-fun x1 () x0)\n    (declare-fun x3 () x0)\n    (declare-fun x5 () x0)\n    (declare-fun x6 () x0)\n    (assert (not (not (= x1 x3))))\n    (assert (=> (= x1 x3) (not (= x5 x6))))\n    (assert (not (= x2 x4)))\n    (assert (= x7 x8))\n    (check-sat)\n    (exit)\n    ")$assertions })), "[Eq(x1, x3), Implies(Eq(x1, x3), Ne(x5, x6)), Ne(x2, x4), Eq(x7, x8)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (set-logic LRA)\n    (declare-fun x () Real)\n    (declare-fun y () Real)\n    (assert (let ((a (+ x 1)) (b (- y 1))) (and (> a 0) (< b 0))))\n    ")$assertions })), "[Q.real(x), Q.real(y), (x + 1 > 0) & (y - 1 < 0)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (set-logic LRA)\n    (define-fun max ((a Real) (b Real)) Real (ite (> a b) a b))\n    (declare-fun x () Real)\n    (declare-fun y () Real)\n    (assert (= (max x y) 10))\n    ")$assertions })), "[Q.real(x), Q.real(y), Eq(Piecewise((x, x > y), (y, True)), 10)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-fun |a b| () String)\n    (assert (= |a b| \"hello world\"))\n    ")$assertions })), "[Eq(a b, hello world)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (assert (= #x0A 10))\n    (assert (= #b1010 10))\n    ")$assertions })), "[True, True]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-sort A 0)\n    (declare-fun y () A)\n    ")$assertions })), "[]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (define-fun-rec f ((x Int)) Int (ite (> x 0) 1 0))\n    ")$assertions })), "ERR NotImplementedError: Recursive macros are not supported in SymPy because AST generation evaluates bottom-up.")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-fun x () Real)\n    (declare-fun y () Real)\n    (declare-fun z () Real)\n    (assert (= (/ x y z) 1.0))\n    (assert false)\n    ")$assertions })), "[Q.real(x), Q.real(y), Q.real(z), Eq(x/(y*z), 1.0), False]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (assert (and true))\n    ")$assertions })), "[True]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-fun x () Int)\n    (declare-fun y () Int)\n    (declare-fun z () Int)\n    (assert (= x y z))\n    (assert (distinct x y z))\n    (assert (xor x y z))\n    (assert (= (- x) (abs y)))\n    ")$assertions })), "[Q.integer(x), Q.integer(y), Q.integer(z), Eq(x, y) & Eq(y, z), Ne(x, y) & Ne(x, z) & Ne(y, z), x ^ y ^ z, Eq(-x, Abs(y))]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (set-info :smt-lib-version 2.6)\n    (set-logic QF_LIA)\n    (declare-fun x () Int)\n    (push 1)\n    (assert (> x 0))\n    (pop 1)\n    (get-info :name)\n    ")$assertions })), "[Q.integer(x)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-fun x () Int)\n    (assert (= x 10))\n    (reset-assertions)\n    ")$assertions })), "[]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-fun y () Int)\n    (reset)\n    ")$assertions })), "[]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-sort S 0)\n    (declare-fun a () S)\n    (declare-fun b () S)\n    (declare-fun c () S)\n    (assert (=> a b c))\n    ")$assertions })), "[Implies(a, Implies(b, c))]")
  expect_identical(.pyrun(local({ parse_smtlib("(assert (< 1 2 0))")$assertions })), "[False]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-fun a () Int)\n    (declare-fun b () Int)\n    (define-fun f ((a Int)(b Int)) Int (- a b))\n    (assert (= (f b a) 0))\n    ")$assertions })), "[Q.integer(a), Q.integer(b), Eq(-a + b, 0)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-fun a () Int)\n    (declare-fun b () Int)\n    (assert (let ((b a) (a 1)) (= b 1)))\n    ")$assertions })), "[Q.integer(a), Q.integer(b), Eq(a, 1)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (set-option :print-success false)\n    (set-logic QF_LRA)\n    (declare-fun x1 () Real)\n    (declare-fun x3 () Real)\n    (declare-fun x2 () Real)\n    (assert (not (=> (and (<= x3 (+ 10 x2)) (and (and (>= x3 1) (and (= 10 x3) (= 0 x1))) (<= x3 12))) (or (>= 10 x3) (= x1 3)))))\n    (check-sat)\n    (exit)\n    ")$assertions })), "[Q.real(x1), Q.real(x3), Q.real(x2), ~(Implies(Eq(0, x1) & Eq(10, x3) & (x3 >= 1) & (x3 <= 12) & (x3 <= x2 + 10), Eq(x1, 3) | (10 >= x3)))]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (set-option :print-success false)\n    (set-logic QF_LRA)\n    (declare-fun x1 () Real)\n    (declare-fun x3 () Real)\n    (declare-fun x2 () Real)\n    (assert (not (=> (and (and (and (>= x3 1) (and (= 5 x3) (= 2 x2))) (<= x3 12)) (>= x3 (- 5 (* 2 x1)))) (or (>= x3 5) (= x2 1)))))\n    (check-sat)\n    (exit)\n    ")$assertions })), "[Q.real(x1), Q.real(x3), Q.real(x2), ~(Implies(Eq(2, x2) & Eq(5, x3) & (x3 >= 1) & (x3 <= 12) & (x3 >= 5 - 2*x1), Eq(x2, 1) | (x3 >= 5)))]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (set-logic QF_LRA)\n    (declare-fun x () Real)\n    (declare-fun y () Real)\n    (assert (distinct x y))\n    (assert (= (/ x y) 10.0))\n    ")$assertions })), "[Q.real(x), Q.real(y), Ne(x, y), Eq(x/y, 10.0)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (set-logic QF_LRA)\n    (declare-fun x () Real)\n    (declare-fun y () Real)\n    (declare-fun z () Real)\n    (assert (<= x y z))\n    (assert (> x y z))\n    (assert (>= x y z))\n    ")$assertions })), "[Q.real(x), Q.real(y), Q.real(z), (x <= y) & (y <= z), (x > y) & (y > z), (x >= y) & (y >= z)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (define-const n Int 5)\n    (declare-const x Int)\n    (assert (= x n))\n    ")$assertions })), "[Q.integer(x), Eq(x, 5)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (define-sort MyReal () Real)\n    (define-sort AlsoMyReal () MyReal)\n    (declare-const x AlsoMyReal)\n    ")$assertions })), "[Q.real(x)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-const x Int)\n    (assert (= x (_ bv5 32)))\n    ")$assertions })), "[Q.integer(x), Eq(x, 5)]")
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= |x| \"say \"\"hi\"\"\"))")$assertions })), "[Eq(x, say \"hi\")]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-sort S 0)\n    (declare-fun f (S) S)\n    (declare-fun x () S)\n    (assert (= (f (as x S)) x))\n    ")$assertions })), "[Eq(f(x), x)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (assert (= (div 7 2) 0))\n    (assert (= (div (- 7) 2) 0))\n    (assert (= (div 7 (- 2)) 0))\n    (assert (= (mod 7 2) 0))\n    (assert (= (mod (- 7) 2) 0))\n    (assert (= (mod 7 (- 2)) 0))\n    ")$assertions })), "[False, False, False, False, False, False]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-const a Int)\n    (declare-const b Int)\n    (assert (= (mod a b) 0))\n    ")$assertions })), "[Q.integer(a), Q.integer(b), Eq(Mod(a, Abs(b)), 0)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-const x Real)\n    (assert (= (to_int x) 0))\n    (assert (= (to_real 5) 5.0))\n    (assert (is_int x))\n    ")$assertions })), "[Q.real(x), Eq(floor(x), 0), True, Eq(x, floor(x))]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (assert true)\n    (assert (let ((a 1)) true))\n    ")$assertions })), "[True, True]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-const True Int)\n    (assert (= True 1))\n    ")$assertions })), "[Q.integer(True), Eq(True, 1)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-const p Int)\n    (check-sat-assuming (p (not p)))\n    (echo \"hello\")\n    ")$assertions })), "[Q.integer(p)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-const x Real)\n    (assert (= x -1))\n    (assert (< x -2.5))\n    (assert (= (- 1) -1))\n    ")$assertions })), "[Q.real(x), Eq(x, -1), x < -2.5, True]")
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const a-1 Int)")$assertions })), "[Q.integer(a-1)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n    (declare-const x Real)\n    (assert (= (pow x 2) 4))\n    (assert (= (^ x 3) 8))\n    ")$assertions })), "[Q.real(x), Eq(x^2, 4), Eq(x^3, 8)]")
  expect_identical(.pyrun(local({ parse_smtlib("\n        (set-logic LRA)\n        (declare-fun x () Real)\n        (assert (forall ((a Real) (b Real)) (=> (> a b) (> (+ a x) (+ b x)))))\n        ")$assertions })), "ERR NotImplementedError: Quantifiers (forall) are not supported since SymPy has no representation for them")
  expect_identical(.pyrun(local({ parse_smtlib("\n        (set-logic LRA)\n        (declare-fun x () Real)\n        (assert (exists ((c Real)) (= c (* x 2))))\n        ")$assertions })), "ERR NotImplementedError: Quantifiers (exists) are not supported since SymPy has no representation for them")
  expect_identical(.pyrun(local({ parse_smtlib("\n        (set-logic QF_BV)\n        (declare-fun a () (_ BitVec 32))\n        (declare-fun b () (_ BitVec 32))\n        (assert (= (bvadd a b) #x00000000))\n        ")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvadd' has no SymPy equivalent and is not supported")
  expect_identical(.pyrun(local({ parse_smtlib("(assert (unknown-op x y))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: unknown-op")
  expect_identical(.pyrun(local({ parse_smtlib("(asssert (= x 1))")$assertions })), "ERR UnknownSMTLibCommandError: Unknown SMT-LIB command: asssert")
  expect_identical(.pyrun(local({ parse_smtlib("(some-unknown-command 1 2)")$assertions })), "ERR UnknownSMTLibCommandError: Unknown SMT-LIB command: some-unknown-command")
  expect_identical(.pyrun(local({ parse_smtlib("(assert (match 10 ((x 20))))")$assertions })), "ERR NotImplementedError: Match expressions are not fully supported yet.")
  expect_identical(.pyrun(local({ parse_smtlib("\n        (declare-sort S 0)\n        (declare-fun x () S)\n        (assert (= ((_ extract 7 0) x) 0))\n        ")$assertions })), "ERR NotImplementedError: Indexed identifier (_ extract ...) is not supported")
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= ((as const Array) 1) 1))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'const' has no SymPy equivalent and is not supported")
  expect_identical(.pyrun(local({ parse_smtlib("(define-sort MyArray (T) (Array Int T))")$assertions })), "ERR NotImplementedError: Parametric sort definitions (define-sort with parameters) are not supported")
  expect_identical(.pyrun(local({ parse_smtlib("(assert (AND x y))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: AND")
  expect_identical(.pyrun(local({ parse_smtlib("(declare-datatype Color ((red) (green) (blue)))")$assertions })), "ERR NotImplementedError: Algebraic datatypes (declare-datatype) are not supported")
  expect_identical(.pyrun(local({ parse_smtlib("\n        (declare-datatypes ((Color 0)) (((red) (green) (blue))))\n        ")$assertions })), "ERR NotImplementedError: Algebraic datatypes (declare-datatypes) are not supported")
})

test_that("smtlib: test symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (set-option :print-success false)\n    (set-logic QF_UF)\n    (declare-sort x5 0)\n    (declare-fun x4 () x5)\n    (declare-fun x6 () x5)\n    (declare-fun x7 () x5)\n    (declare-fun x1 () x5)\n    (declare-fun x3 () x5)\n    (assert (= x6 x7))\n    (assert (not (= x4 x3)))\n    (assert (not (or (= x4 x1) (= x1 x3))))\n    (check-sat)\n    (exit)\n    ")$symbols) })), "{'x4': x4, 'x6': x6, 'x7': x7, 'x1': x1, 'x3': x3}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (set-option :print-success false)\n    (set-logic QF_UF)\n    (declare-sort x0 0)\n    (declare-fun x2 () x0)\n    (declare-fun x4 () x0)\n    (declare-fun x7 () x0)\n    (declare-fun x8 () x0)\n    (declare-fun x1 () x0)\n    (declare-fun x3 () x0)\n    (declare-fun x5 () x0)\n    (declare-fun x6 () x0)\n    (assert (not (not (= x1 x3))))\n    (assert (=> (= x1 x3) (not (= x5 x6))))\n    (assert (not (= x2 x4)))\n    (assert (= x7 x8))\n    (check-sat)\n    (exit)\n    ")$symbols) })), "{'x2': x2, 'x4': x4, 'x7': x7, 'x8': x8, 'x1': x1, 'x3': x3, 'x5': x5, 'x6': x6}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (set-logic LRA)\n    (declare-fun x () Real)\n    (declare-fun y () Real)\n    (assert (let ((a (+ x 1)) (b (- y 1))) (and (> a 0) (< b 0))))\n    ")$symbols) })), "{'x': x, 'y': y}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (set-logic LRA)\n    (define-fun max ((a Real) (b Real)) Real (ite (> a b) a b))\n    (declare-fun x () Real)\n    (declare-fun y () Real)\n    (assert (= (max x y) 10))\n    ")$symbols) })), "{'x': x, 'y': y}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-fun |a b| () String)\n    (assert (= |a b| \"hello world\"))\n    ")$symbols) })), "{'a b': a b}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (assert (= #x0A 10))\n    (assert (= #b1010 10))\n    ")$symbols) })), "{}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-sort A 0)\n    (declare-fun y () A)\n    ")$symbols) })), "{'y': y}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-fun x () Real)\n    (declare-fun y () Real)\n    (declare-fun z () Real)\n    (assert (= (/ x y z) 1.0))\n    (assert false)\n    ")$symbols) })), "{'x': x, 'y': y, 'z': z}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (assert (and true))\n    ")$symbols) })), "{}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-fun x () Int)\n    (declare-fun y () Int)\n    (declare-fun z () Int)\n    (assert (= x y z))\n    (assert (distinct x y z))\n    (assert (xor x y z))\n    (assert (= (- x) (abs y)))\n    ")$symbols) })), "{'x': x, 'y': y, 'z': z}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (set-info :smt-lib-version 2.6)\n    (set-logic QF_LIA)\n    (declare-fun x () Int)\n    (push 1)\n    (assert (> x 0))\n    (pop 1)\n    (get-info :name)\n    ")$symbols) })), "{'x': x}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-fun x () Int)\n    (assert (= x 10))\n    (reset-assertions)\n    ")$symbols) })), "{'x': x}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-fun y () Int)\n    (reset)\n    ")$symbols) })), "{}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-sort S 0)\n    (declare-fun a () S)\n    (declare-fun b () S)\n    (declare-fun c () S)\n    (assert (=> a b c))\n    ")$symbols) })), "{'a': a, 'b': b, 'c': c}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (< 1 2 0))")$symbols) })), "{}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-fun a () Int)\n    (declare-fun b () Int)\n    (define-fun f ((a Int)(b Int)) Int (- a b))\n    (assert (= (f b a) 0))\n    ")$symbols) })), "{'a': a, 'b': b}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-fun a () Int)\n    (declare-fun b () Int)\n    (assert (let ((b a) (a 1)) (= b 1)))\n    ")$symbols) })), "{'a': a, 'b': b}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (set-option :print-success false)\n    (set-logic QF_LRA)\n    (declare-fun x1 () Real)\n    (declare-fun x3 () Real)\n    (declare-fun x2 () Real)\n    (assert (not (=> (and (<= x3 (+ 10 x2)) (and (and (>= x3 1) (and (= 10 x3) (= 0 x1))) (<= x3 12))) (or (>= 10 x3) (= x1 3)))))\n    (check-sat)\n    (exit)\n    ")$symbols) })), "{'x1': x1, 'x3': x3, 'x2': x2}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (set-option :print-success false)\n    (set-logic QF_LRA)\n    (declare-fun x1 () Real)\n    (declare-fun x3 () Real)\n    (declare-fun x2 () Real)\n    (assert (not (=> (and (and (and (>= x3 1) (and (= 5 x3) (= 2 x2))) (<= x3 12)) (>= x3 (- 5 (* 2 x1)))) (or (>= x3 5) (= x2 1)))))\n    (check-sat)\n    (exit)\n    ")$symbols) })), "{'x1': x1, 'x3': x3, 'x2': x2}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (set-logic QF_LRA)\n    (declare-fun x () Real)\n    (declare-fun y () Real)\n    (assert (distinct x y))\n    (assert (= (/ x y) 10.0))\n    ")$symbols) })), "{'x': x, 'y': y}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (set-logic QF_LRA)\n    (declare-fun x () Real)\n    (declare-fun y () Real)\n    (declare-fun z () Real)\n    (assert (<= x y z))\n    (assert (> x y z))\n    (assert (>= x y z))\n    ")$symbols) })), "{'x': x, 'y': y, 'z': z}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (define-const n Int 5)\n    (declare-const x Int)\n    (assert (= x n))\n    ")$symbols) })), "{'n': 5, 'x': x}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (define-sort MyReal () Real)\n    (define-sort AlsoMyReal () MyReal)\n    (declare-const x AlsoMyReal)\n    ")$symbols) })), "{'x': x}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-const x Int)\n    (assert (= x (_ bv5 32)))\n    ")$symbols) })), "{'x': x}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= |x| \"say \"\"hi\"\"\"))")$symbols) })), "{}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-sort S 0)\n    (declare-fun f (S) S)\n    (declare-fun x () S)\n    (assert (= (f (as x S)) x))\n    ")$symbols) })), "{'x': x}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (assert (= (div 7 2) 0))\n    (assert (= (div (- 7) 2) 0))\n    (assert (= (div 7 (- 2)) 0))\n    (assert (= (mod 7 2) 0))\n    (assert (= (mod (- 7) 2) 0))\n    (assert (= (mod 7 (- 2)) 0))\n    ")$symbols) })), "{}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-const a Int)\n    (declare-const b Int)\n    (assert (= (mod a b) 0))\n    ")$symbols) })), "{'a': a, 'b': b}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-const x Real)\n    (assert (= (to_int x) 0))\n    (assert (= (to_real 5) 5.0))\n    (assert (is_int x))\n    ")$symbols) })), "{'x': x}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (assert true)\n    (assert (let ((a 1)) true))\n    ")$symbols) })), "{}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-const True Int)\n    (assert (= True 1))\n    ")$symbols) })), "{'True': True}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-const p Int)\n    (check-sat-assuming (p (not p)))\n    (echo \"hello\")\n    ")$symbols) })), "{'p': p}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-const x Real)\n    (assert (= x -1))\n    (assert (< x -2.5))\n    (assert (= (- 1) -1))\n    ")$symbols) })), "{'x': x}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const a-1 Int)")$symbols) })), "{'a-1': a-1}")
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("\n    (declare-const x Real)\n    (assert (= (pow x 2) 4))\n    (assert (= (^ x 3) 8))\n    ")$symbols) })), "{'x': x}")
})

test_that("smtlib: test syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (= x 1)")) })), "ERR SMTLibSyntaxError")
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (= x 1)))")) })), "ERR SMTLibSyntaxError")
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (= x \"unclosed))")) })), "ERR SMTLibSyntaxError")
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-fun |unclosed () S)")) })), "ERR SMTLibSyntaxError")
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-sort B 1 2)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: roundtrip 0 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (= (+ x (* 2 y)) 2))")$assertions })), "[Q.real(x), Q.real(y), Eq(x + 2*y, 2)]")
})

test_that("smtlib: roundtrip 0 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (= (+ x (* 2 y)) 2))")$symbols) })), "{'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 1 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (< (* x (pow y -1)) 5))")$assertions })), "[Q.real(x), Q.real(y), x/y < 5]")
})

test_that("smtlib: roundtrip 1 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (< (* x (pow y -1)) 5))")$symbols) })), "{'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 2 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(assert (= (pow x 2) 4))")$assertions })), "[Q.real(x), Eq(x^2, 4)]")
})

test_that("smtlib: roundtrip 2 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(assert (= (pow x 2) 4))")$symbols) })), "{'x': x}")
})

test_that("smtlib: roundtrip 3 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const n Int)\n(declare-const x Real)\n(assert (>= n 5))\n(assert (= x (/ 1 3)))\n(assert (<= x 2.5))")$assertions })), "[Q.integer(n), Q.real(x), n >= 5, Eq(x, 1/3), x <= 2.5]")
})

test_that("smtlib: roundtrip 3 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const n Int)\n(declare-const x Real)\n(assert (>= n 5))\n(assert (= x (/ 1 3)))\n(assert (<= x 2.5))")$symbols) })), "{'n': n, 'x': x}")
})

test_that("smtlib: roundtrip 4 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(assert (>= x -5))\n(assert (= x -2.5))")$assertions })), "[Q.real(x), x >= -5, Eq(x, -2.5)]")
})

test_that("smtlib: roundtrip 4 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(assert (>= x -5))\n(assert (= x -2.5))")$symbols) })), "{'x': x}")
})

test_that("smtlib: roundtrip 5 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const b1 Bool)\n(declare-const b2 Bool)\n(declare-const x Real)\n(assert (=> (= x 1) (xor b1 b2)))")$assertions })), "[Q.real(x), Implies(Eq(x, 1), b1 ^ b2)]")
})

test_that("smtlib: roundtrip 5 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const b1 Bool)\n(declare-const b2 Bool)\n(declare-const x Real)\n(assert (=> (= x 1) (xor b1 b2)))")$symbols) })), "{'b1': b1, 'b2': b2, 'x': x}")
})

test_that("smtlib: roundtrip 6 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const b1 Bool)\n(declare-const b2 Bool)\n(declare-const x Real)\n(assert (not (and b1 (or b2 (> x 0)))))")$assertions })), "[Q.real(x), ~(b1 & (b2 | (x > 0)))]")
})

test_that("smtlib: roundtrip 6 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const b1 Bool)\n(declare-const b2 Bool)\n(declare-const x Real)\n(assert (not (and b1 (or b2 (> x 0)))))")$symbols) })), "{'b1': b1, 'b2': b2, 'x': x}")
})

test_that("smtlib: roundtrip 7 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (not (= x y)))")$assertions })), "[Q.real(x), Q.real(y), Ne(x, y)]")
})

test_that("smtlib: roundtrip 7 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (not (= x y)))")$symbols) })), "{'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 8 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(assert (< (ite (> x 0) x (* -1 x)) 2))")$assertions })), "[Q.real(x), Piecewise((x, x > 0), (-x, True)) < 2]")
})

test_that("smtlib: roundtrip 8 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(assert (< (ite (> x 0) x (* -1 x)) 2))")$symbols) })), "{'x': x}")
})

test_that("smtlib: roundtrip 9 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(assert (= (abs x) 1))")$assertions })), "[Q.real(x), Eq(Abs(x), 1)]")
})

test_that("smtlib: roundtrip 9 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(assert (= (abs x) 1))")$symbols) })), "{'x': x}")
})

test_that("smtlib: roundtrip 10 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const n Int)\n(declare-const x Real)\n(declare-const y Real)\n(assert (= (+ (* -3 x) (* x y)) n))")$assertions })), "[Q.integer(n), Q.real(x), Q.real(y), Eq(x*y - 3*x, n)]")
})

test_that("smtlib: roundtrip 10 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const n Int)\n(declare-const x Real)\n(declare-const y Real)\n(assert (= (+ (* -3 x) (* x y)) n))")$symbols) })), "{'n': n, 'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 11 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (> (+ (pow x 3) (pow y 2)) 1))")$assertions })), "[Q.real(x), Q.real(y), x^3 + y^2 > 1]")
})

test_that("smtlib: roundtrip 11 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (> (+ (pow x 3) (pow y 2)) 1))")$symbols) })), "{'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 12 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (= (+ (* (/ 1 2) y) (* (/ 1 3) x)) 1))")$assertions })), "[Q.real(x), Q.real(y), Eq(x/3 + y/2, 1)]")
})

test_that("smtlib: roundtrip 12 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (= (+ (* (/ 1 2) y) (* (/ 1 3) x)) 1))")$symbols) })), "{'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 13 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const n Int)\n(declare-const x Real)\n(declare-const y Real)\n(assert (or (> x 1) (< y -1)))\n(assert (= n 2))")$assertions })), "[Q.integer(n), Q.real(x), Q.real(y), (x > 1) | (y < -1), Eq(n, 2)]")
})

test_that("smtlib: roundtrip 13 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const n Int)\n(declare-const x Real)\n(declare-const y Real)\n(assert (or (> x 1) (< y -1)))\n(assert (= n 2))")$symbols) })), "{'n': n, 'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 14 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (<= (+ x (* -1 y)) 0))")$assertions })), "[Q.real(x), Q.real(y), x - y <= 0]")
})

test_that("smtlib: roundtrip 14 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (<= (+ x (* -1 y)) 0))")$symbols) })), "{'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 15 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (= (+ (/ 3 4) (* -7 y) (* 2 x)) 0))")$assertions })), "[Q.real(x), Q.real(y), Eq(2*x - 7*y + 3/4, 0)]")
})

test_that("smtlib: roundtrip 15 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (= (+ (/ 3 4) (* -7 y) (* 2 x)) 0))")$symbols) })), "{'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 16 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (and (> x 0) (> y 0) (< (+ x y) 1)))")$assertions })), "[Q.real(x), Q.real(y), (x > 0) & (y > 0) & (x + y < 1)]")
})

test_that("smtlib: roundtrip 16 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (and (> x 0) (> y 0) (< (+ x y) 1)))")$symbols) })), "{'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 17 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const b1 Bool)\n(declare-const b2 Bool)\n(assert (=> b1 b2))\n(assert (xor b1 b2))")$assertions })), "[Implies(b1, b2), b1 ^ b2]")
})

test_that("smtlib: roundtrip 17 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const b1 Bool)\n(declare-const b2 Bool)\n(assert (=> b1 b2))\n(assert (xor b1 b2))")$symbols) })), "{'b1': b1, 'b2': b2}")
})

test_that("smtlib: roundtrip 18 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (= (abs (+ x (* -1 y))) 2))")$assertions })), "[Q.real(x), Q.real(y), Eq(Abs(x - y), 2)]")
})

test_that("smtlib: roundtrip 18 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)\n(declare-const y Real)\n(assert (= (abs (+ x (* -1 y))) 2))")$symbols) })), "{'x': x, 'y': y}")
})

test_that("smtlib: roundtrip 19 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const b1 Real)\n(declare-const x Real)\n(declare-const y Real)\n(assert (= (ite b1 x y) 1))")$assertions })), "[Q.real(b1), Q.real(x), Q.real(y), Eq(Piecewise((x, b1), (y, True)), 1)]")
})

test_that("smtlib: roundtrip 19 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const b1 Real)\n(declare-const x Real)\n(declare-const y Real)\n(assert (= (ite b1 x y) 1))")$symbols) })), "{'b1': b1, 'x': x, 'y': y}")
})

test_that("smtlib: roundtrip f assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const n Int)\n(declare-fun f (Int) Int)\n(assert (= (f n) 0))")$assertions })), "[Q.integer(n), Eq(f(n), 0)]")
})

test_that("smtlib: roundtrip f symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const n Int)\n(declare-fun f (Int) Int)\n(assert (= (f n) 0))")$symbols) })), "{'n': n}")
})

test_that("smtlib: extra 0 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("")$assertions })), "[]")
})

test_that("smtlib: extra 0 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("")$symbols) })), "{}")
})

test_that("smtlib: extra 1 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("; only a comment")$assertions })), "[]")
})

test_that("smtlib: extra 1 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("; only a comment")$symbols) })), "{}")
})

test_that("smtlib: extra 2 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(check-sat)")$assertions })), "[]")
})

test_that("smtlib: extra 2 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(check-sat)")$symbols) })), "{}")
})

test_that("smtlib: extra 3 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(exit)")$assertions })), "[]")
})

test_that("smtlib: extra 3 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(exit)")$symbols) })), "{}")
})

test_that("smtlib: extra 4 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(get-model)")$assertions })), "[]")
})

test_that("smtlib: extra 4 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(get-model)")$symbols) })), "{}")
})

test_that("smtlib: extra 5 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Int)")$assertions })), "[Q.integer(x)]")
})

test_that("smtlib: extra 5 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Int)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 6 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)")$assertions })), "[Q.real(x)]")
})

test_that("smtlib: extra 6 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 7 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Bool)")$assertions })), "[]")
})

test_that("smtlib: extra 7 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Bool)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 8 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x (Array Int Int))")$assertions })), "[]")
})

test_that("smtlib: extra 8 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x (Array Int Int))")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 9 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x (_ BitVec 8))")$assertions })), "[]")
})

test_that("smtlib: extra 9 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x (_ BitVec 8))")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 10 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Int)(declare-const x Real)")$assertions })), "[Q.integer(x), Q.real(x)]")
})

test_that("smtlib: extra 10 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Int)(declare-const x Real)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 11 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const |a b| Int)(assert (> |a b| 0))")$assertions })), "[Q.integer(a b), a b > 0]")
})

test_that("smtlib: extra 11 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const |a b| Int)(assert (> |a b| 0))")$symbols) })), "{'a b': a b}")
})

test_that("smtlib: extra 12 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x 01))")$assertions })), "[False]")
})

test_that("smtlib: extra 12 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x 01))")$symbols) })), "{}")
})

test_that("smtlib: extra 13 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x 1x))")$symbols) })), "{}")
})

test_that("smtlib: extra 14 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const -5 Int)")$assertions })), "[Q.integer(-5)]")
})

test_that("smtlib: extra 14 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const -5 Int)")$symbols) })), "{'-5': -5}")
})

test_that("smtlib: extra 15 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const -5 Int)(assert (= -5 0))")$assertions })), "[Q.integer(-5), False]")
})

test_that("smtlib: extra 15 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const -5 Int)(assert (= -5 0))")$symbols) })), "{'-5': -5}")
})

test_that("smtlib: extra 16 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x -0))")$assertions })), "[Eq(x, 0)]")
})

test_that("smtlib: extra 16 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x -0))")$symbols) })), "{}")
})

test_that("smtlib: extra 17 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x 1.5))")$assertions })), "[Eq(x, 1.5)]")
})

test_that("smtlib: extra 17 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x 1.5))")$symbols) })), "{}")
})

test_that("smtlib: extra 18 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x -0.25))")$assertions })), "[Eq(x, -0.25)]")
})

test_that("smtlib: extra 18 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x -0.25))")$symbols) })), "{}")
})

test_that("smtlib: extra 19 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x 0.1))")$assertions })), "[Eq(x, 0.1)]")
})

test_that("smtlib: extra 19 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x 0.1))")$symbols) })), "{}")
})

test_that("smtlib: extra 20 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x 12345678901234567890))")$assertions })), "[Eq(x, 12345678901234567890)]")
})

test_that("smtlib: extra 20 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x 12345678901234567890))")$symbols) })), "{}")
})

test_that("smtlib: extra 21 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x 3.14159265358979323846))")$assertions })), "[Eq(x, 3.14159265358979323846)]")
})

test_that("smtlib: extra 21 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x 3.14159265358979323846))")$symbols) })), "{}")
})

test_that("smtlib: extra 22 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x #xFF))")$assertions })), "[Eq(x, 255)]")
})

test_that("smtlib: extra 22 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x #xFF))")$symbols) })), "{}")
})

test_that("smtlib: extra 23 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x #b0))")$assertions })), "[Eq(x, 0)]")
})

test_that("smtlib: extra 23 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x #b0))")$symbols) })), "{}")
})

test_that("smtlib: extra 24 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x #b111))")$assertions })), "[Eq(x, 7)]")
})

test_that("smtlib: extra 24 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x #b111))")$symbols) })), "{}")
})

test_that("smtlib: extra 25 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x \"a;b\"))")$assertions })), "[Eq(x, a;b)]")
})

test_that("smtlib: extra 25 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x \"a;b\"))")$symbols) })), "{}")
})

test_that("smtlib: extra 26 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x \"\"))")$assertions })), "[Eq(x, )]")
})

test_that("smtlib: extra 26 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x \"\"))")$symbols) })), "{}")
})

test_that("smtlib: extra 27 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x \"\"\"\"))")$assertions })), "[Eq(x, \")]")
})

test_that("smtlib: extra 27 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x \"\"\"\"))")$symbols) })), "{}")
})

test_that("smtlib: extra 28 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x \"x\"\"y\"))")$assertions })), "[Eq(x, x\"y)]")
})

test_that("smtlib: extra 28 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x \"x\"\"y\"))")$symbols) })), "{}")
})

test_that("smtlib: extra 29 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (+ 1 2 3))")$assertions })), "[6]")
})

test_that("smtlib: extra 29 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (+ 1 2 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 30 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (- 10 1 2 3))")$assertions })), "[4]")
})

test_that("smtlib: extra 30 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (- 10 1 2 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 31 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (- 5))")$assertions })), "[-5]")
})

test_that("smtlib: extra 31 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (- 5))")$symbols) })), "{}")
})

test_that("smtlib: extra 32 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (* 2 3 4))")$assertions })), "[24]")
})

test_that("smtlib: extra 32 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (* 2 3 4))")$symbols) })), "{}")
})

test_that("smtlib: extra 33 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (/ 1 2))")$assertions })), "[1/2]")
})

test_that("smtlib: extra 33 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (/ 1 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 34 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (/ 12 2 3))")$assertions })), "[2]")
})

test_that("smtlib: extra 34 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (/ 12 2 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 35 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (/ x 2))")$assertions })), "[x/2]")
})

test_that("smtlib: extra 35 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (/ x 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 36 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (+ x x))")$assertions })), "[2*x]")
})

test_that("smtlib: extra 36 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (+ x x))")$symbols) })), "{}")
})

test_that("smtlib: extra 37 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (* x x))")$assertions })), "[x^2]")
})

test_that("smtlib: extra 37 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (* x x))")$symbols) })), "{}")
})

test_that("smtlib: extra 38 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (- x x))")$assertions })), "[0]")
})

test_that("smtlib: extra 38 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (- x x))")$symbols) })), "{}")
})

test_that("smtlib: extra 39 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (+ x (* 2 y)))")$assertions })), "[x + 2*y]")
})

test_that("smtlib: extra 39 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (+ x (* 2 y)))")$symbols) })), "{}")
})

test_that("smtlib: extra 40 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (div 7 2 3))")$assertions })), "[1]")
})

test_that("smtlib: extra 40 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (div 7 2 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 41 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (div x 2))")$assertions })), "[x/2 - Mod(x, 2)/2]")
})

test_that("smtlib: extra 41 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (div x 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 42 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (div x y))")$assertions })), "[(x - Mod(x, Abs(y)))/y]")
})

test_that("smtlib: extra 42 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (div x y))")$symbols) })), "{}")
})

test_that("smtlib: extra 43 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (div (- 9) (- 4)))")$assertions })), "[3]")
})

test_that("smtlib: extra 43 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (div (- 9) (- 4)))")$symbols) })), "{}")
})

test_that("smtlib: extra 44 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (div 0 5))")$assertions })), "[0]")
})

test_that("smtlib: extra 44 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (div 0 5))")$symbols) })), "{}")
})

test_that("smtlib: extra 45 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (mod 10 3))")$assertions })), "[1]")
})

test_that("smtlib: extra 45 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (mod 10 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 46 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (mod (- 10) 3))")$assertions })), "[2]")
})

test_that("smtlib: extra 46 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (mod (- 10) 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 47 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (mod 10 (- 3)))")$assertions })), "[1]")
})

test_that("smtlib: extra 47 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (mod 10 (- 3)))")$symbols) })), "{}")
})

test_that("smtlib: extra 48 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (mod (- 10) (- 3)))")$assertions })), "[2]")
})

test_that("smtlib: extra 48 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (mod (- 10) (- 3)))")$symbols) })), "{}")
})

test_that("smtlib: extra 49 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (mod x 3))")$assertions })), "[Mod(x, 3)]")
})

test_that("smtlib: extra 49 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (mod x 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 50 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (abs (- 4)))")$assertions })), "[4]")
})

test_that("smtlib: extra 50 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (abs (- 4)))")$symbols) })), "{}")
})

test_that("smtlib: extra 51 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (abs x))")$assertions })), "[Abs(x)]")
})

test_that("smtlib: extra 51 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (abs x))")$symbols) })), "{}")
})

test_that("smtlib: extra 52 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (to_real x))")$assertions })), "[x]")
})

test_that("smtlib: extra 52 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (to_real x))")$symbols) })), "{}")
})

test_that("smtlib: extra 53 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (to_int 2.5))")$assertions })), "[2]")
})

test_that("smtlib: extra 53 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (to_int 2.5))")$symbols) })), "{}")
})

test_that("smtlib: extra 54 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (to_int (- 2.5)))")$assertions })), "[-3]")
})

test_that("smtlib: extra 54 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (to_int (- 2.5)))")$symbols) })), "{}")
})

test_that("smtlib: extra 55 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (is_int 2))")$assertions })), "[True]")
})

test_that("smtlib: extra 55 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (is_int 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 56 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (is_int 2.5))")$assertions })), "[False]")
})

test_that("smtlib: extra 56 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (is_int 2.5))")$symbols) })), "{}")
})

test_that("smtlib: extra 57 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (to_int x))")$assertions })), "[floor(x)]")
})

test_that("smtlib: extra 57 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (to_int x))")$symbols) })), "{}")
})

test_that("smtlib: extra 58 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (pow 2 10))")$assertions })), "[1024]")
})

test_that("smtlib: extra 58 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (pow 2 10))")$symbols) })), "{}")
})

test_that("smtlib: extra 59 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (^ x y))")$assertions })), "[x^y]")
})

test_that("smtlib: extra 59 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (^ x y))")$symbols) })), "{}")
})

test_that("smtlib: extra 60 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (pow x 0.5))")$assertions })), "[x^0.5]")
})

test_that("smtlib: extra 60 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (pow x 0.5))")$symbols) })), "{}")
})

test_that("smtlib: extra 61 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (ite true 1 2))")$assertions })), "[1]")
})

test_that("smtlib: extra 61 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (ite true 1 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 62 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (ite false 1 2))")$assertions })), "[2]")
})

test_that("smtlib: extra 62 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (ite false 1 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 63 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (ite (> x 0) x (- x)))")$assertions })), "[Piecewise((x, x > 0), (-x, True))]")
})

test_that("smtlib: extra 63 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (ite (> x 0) x (- x)))")$symbols) })), "{}")
})

test_that("smtlib: extra 64 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (ite b (ite c 1 2) 3))")$assertions })), "[Piecewise((Piecewise((1, c), (2, True)), b), (3, True))]")
})

test_that("smtlib: extra 64 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (ite b (ite c 1 2) 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 65 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (not true))")$assertions })), "[False]")
})

test_that("smtlib: extra 65 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (not true))")$symbols) })), "{}")
})

test_that("smtlib: extra 66 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (not false))")$assertions })), "[True]")
})

test_that("smtlib: extra 66 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (not false))")$symbols) })), "{}")
})

test_that("smtlib: extra 67 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (not (not b)))")$assertions })), "[b]")
})

test_that("smtlib: extra 67 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (not (not b)))")$symbols) })), "{}")
})

test_that("smtlib: extra 68 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (and))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 69 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (and a b))")$assertions })), "[a & b]")
})

test_that("smtlib: extra 69 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (and a b))")$symbols) })), "{}")
})

test_that("smtlib: extra 70 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (and a b c d))")$assertions })), "[a & b & c & d]")
})

test_that("smtlib: extra 70 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (and a b c d))")$symbols) })), "{}")
})

test_that("smtlib: extra 71 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (or a))")$assertions })), "[a]")
})

test_that("smtlib: extra 71 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (or a))")$symbols) })), "{}")
})

test_that("smtlib: extra 72 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (or a b false))")$assertions })), "[a | b]")
})

test_that("smtlib: extra 72 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (or a b false))")$symbols) })), "{}")
})

test_that("smtlib: extra 73 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (and a b true))")$assertions })), "[a & b]")
})

test_that("smtlib: extra 73 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (and a b true))")$symbols) })), "{}")
})

test_that("smtlib: extra 74 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (or a true))")$assertions })), "[True]")
})

test_that("smtlib: extra 74 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (or a true))")$symbols) })), "{}")
})

test_that("smtlib: extra 75 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (xor a b))")$assertions })), "[a ^ b]")
})

test_that("smtlib: extra 75 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (xor a b))")$symbols) })), "{}")
})

test_that("smtlib: extra 76 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (xor true false))")$assertions })), "[True]")
})

test_that("smtlib: extra 76 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (xor true false))")$symbols) })), "{}")
})

test_that("smtlib: extra 77 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (xor a true))")$assertions })), "[~a]")
})

test_that("smtlib: extra 77 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (xor a true))")$symbols) })), "{}")
})

test_that("smtlib: extra 78 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (=> a b))")$assertions })), "[Implies(a, b)]")
})

test_that("smtlib: extra 78 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (=> a b))")$symbols) })), "{}")
})

test_that("smtlib: extra 79 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (=> a b c d))")$assertions })), "[Implies(a, Implies(b, Implies(c, d)))]")
})

test_that("smtlib: extra 79 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (=> a b c d))")$symbols) })), "{}")
})

test_that("smtlib: extra 80 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (=> true b))")$assertions })), "[b]")
})

test_that("smtlib: extra 80 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (=> true b))")$symbols) })), "{}")
})

test_that("smtlib: extra 81 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (=> a))")$assertions })), "[a]")
})

test_that("smtlib: extra 81 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (=> a))")$symbols) })), "{}")
})

test_that("smtlib: extra 82 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= a b))")$assertions })), "[Eq(a, b)]")
})

test_that("smtlib: extra 82 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= a b))")$symbols) })), "{}")
})

test_that("smtlib: extra 83 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= 1 1))")$assertions })), "[True]")
})

test_that("smtlib: extra 83 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= 1 1))")$symbols) })), "{}")
})

test_that("smtlib: extra 84 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= 1 2))")$assertions })), "[False]")
})

test_that("smtlib: extra 84 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= 1 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 85 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= a b c d))")$assertions })), "[Eq(a, b) & Eq(b, c) & Eq(c, d)]")
})

test_that("smtlib: extra 85 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= a b c d))")$symbols) })), "{}")
})

test_that("smtlib: extra 86 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= 1 1 1))")$assertions })), "[True]")
})

test_that("smtlib: extra 86 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= 1 1 1))")$symbols) })), "{}")
})

test_that("smtlib: extra 87 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (distinct a b))")$assertions })), "[Ne(a, b)]")
})

test_that("smtlib: extra 87 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (distinct a b))")$symbols) })), "{}")
})

test_that("smtlib: extra 88 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (distinct a b c d))")$assertions })), "[Ne(a, b) & Ne(a, c) & Ne(a, d) & Ne(b, c) & Ne(b, d) & Ne(c, d)]")
})

test_that("smtlib: extra 88 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (distinct a b c d))")$symbols) })), "{}")
})

test_that("smtlib: extra 89 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (distinct 1 2))")$assertions })), "[True]")
})

test_that("smtlib: extra 89 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (distinct 1 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 90 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (distinct 1 1))")$assertions })), "[False]")
})

test_that("smtlib: extra 90 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (distinct 1 1))")$symbols) })), "{}")
})

test_that("smtlib: extra 91 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (distinct 1 2 3))")$assertions })), "[True]")
})

test_that("smtlib: extra 91 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (distinct 1 2 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 92 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (< x y))")$assertions })), "[x < y]")
})

test_that("smtlib: extra 92 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (< x y))")$symbols) })), "{}")
})

test_that("smtlib: extra 93 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (<= x y))")$assertions })), "[x <= y]")
})

test_that("smtlib: extra 93 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (<= x y))")$symbols) })), "{}")
})

test_that("smtlib: extra 94 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (> x y))")$assertions })), "[x > y]")
})

test_that("smtlib: extra 94 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (> x y))")$symbols) })), "{}")
})

test_that("smtlib: extra 95 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (>= x y))")$assertions })), "[x >= y]")
})

test_that("smtlib: extra 95 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (>= x y))")$symbols) })), "{}")
})

test_that("smtlib: extra 96 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (< 1 2))")$assertions })), "[True]")
})

test_that("smtlib: extra 96 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (< 1 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 97 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (<= 2 2 3))")$assertions })), "[True]")
})

test_that("smtlib: extra 97 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (<= 2 2 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 98 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (> 3 2 1))")$assertions })), "[True]")
})

test_that("smtlib: extra 98 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (> 3 2 1))")$symbols) })), "{}")
})

test_that("smtlib: extra 99 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (>= 3 3 4))")$assertions })), "[False]")
})

test_that("smtlib: extra 99 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (>= 3 3 4))")$symbols) })), "{}")
})

test_that("smtlib: extra 100 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (< x y z w))")$assertions })), "[(x < y) & (y < z) & (z < w)]")
})

test_that("smtlib: extra 100 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (< x y z w))")$symbols) })), "{}")
})

test_that("smtlib: extra 101 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (> x 1 y))")$assertions })), "[(1 > y) & (x > 1)]")
})

test_that("smtlib: extra 101 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (> x 1 y))")$symbols) })), "{}")
})

test_that("smtlib: extra 102 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)(assert (= (* 2 x) (+ x x)))")$assertions })), "[Q.real(x), True]")
})

test_that("smtlib: extra 102 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)(assert (= (* 2 x) (+ x x)))")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 103 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Real)(assert (< (/ x 3) 1.5))")$assertions })), "[Q.real(x), x/3 < 1.5]")
})

test_that("smtlib: extra 103 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Real)(assert (< (/ x 3) 1.5))")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 104 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun f ((a Int)) Int (+ a 1))(assert (= (f 2) 3))")$assertions })), "[True]")
})

test_that("smtlib: extra 104 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun f ((a Int)) Int (+ a 1))(assert (= (f 2) 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 105 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun f ((a Int)) Int (+ a 1))(assert (= (f x) y))")$assertions })), "[Eq(x + 1, y)]")
})

test_that("smtlib: extra 105 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun f ((a Int)) Int (+ a 1))(assert (= (f x) y))")$symbols) })), "{}")
})

test_that("smtlib: extra 106 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun g ((a Int) (b Int)) Bool (> a b))(assert (g 1 2))")$assertions })), "[False]")
})

test_that("smtlib: extra 106 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun g ((a Int) (b Int)) Bool (> a b))(assert (g 1 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 107 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun g ((a Int) (b Int)) Bool (> a b))(assert (g x 2 3))")$assertions })), "[x > 2]")
})

test_that("smtlib: extra 107 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun g ((a Int) (b Int)) Bool (> a b))(assert (g x 2 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 108 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun g ((a Int) (b Int)) Bool (> a b))(assert (g x))")$assertions })), "[x > b]")
})

test_that("smtlib: extra 108 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun g ((a Int) (b Int)) Bool (> a b))(assert (g x))")$symbols) })), "{}")
})

test_that("smtlib: extra 109 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun h () Int 42)(assert (= x h))")$assertions })), "[Eq(x, 42)]")
})

test_that("smtlib: extra 109 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun h () Int 42)(assert (= x h))")$symbols) })), "{'h': 42}")
})

test_that("smtlib: extra 110 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun h () Int (+ y 1))(assert (= x h))")$assertions })), "[Eq(x, y + 1)]")
})

test_that("smtlib: extra 110 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun h () Int (+ y 1))(assert (= x h))")$symbols) })), "{'h': y + 1}")
})

test_that("smtlib: extra 111 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun sq ((a Real)) Real (* a a))(assert (< (sq (sq x)) 16))")$assertions })), "[x^4 < 16]")
})

test_that("smtlib: extra 111 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun sq ((a Real)) Real (* a a))(assert (< (sq (sq x)) 16))")$symbols) })), "{}")
})

test_that("smtlib: extra 112 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun f ((a Int) (a Int)) Int (- a 1))(assert (= (f 5 7) 0))")$assertions })), "[False]")
})

test_that("smtlib: extra 112 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun f ((a Int) (a Int)) Int (- a 1))(assert (= (f 5 7) 0))")$symbols) })), "{}")
})

test_that("smtlib: extra 113 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-fun f (Int) Int)(assert (= (f 1) 2))")$assertions })), "[Eq(f(1), 2)]")
})

test_that("smtlib: extra 113 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-fun f (Int) Int)(assert (= (f 1) 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 114 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-fun f (Int Int) Int)(assert (> (f x y) (f y x)))")$assertions })), "[f(x, y) > f(y, x)]")
})

test_that("smtlib: extra 114 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-fun f (Int Int) Int)(assert (> (f x y) (f y x)))")$symbols) })), "{}")
})

test_that("smtlib: extra 115 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-fun f (Int) Int)(assert (= f 2))")$assertions })), "[Eq(f, 2)]")
})

test_that("smtlib: extra 115 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-fun f (Int) Int)(assert (= f 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 116 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-fun p (Int) Bool)(assert (p 3))")$assertions })), "[p(3)]")
})

test_that("smtlib: extra 116 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-fun p (Int) Bool)(assert (p 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 117 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-fun p (Int) Bool)(assert (not (p x)))")$assertions })), "[~(p(x))]")
})

test_that("smtlib: extra 117 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-fun p (Int) Bool)(assert (not (p x)))")$symbols) })), "{}")
})

test_that("smtlib: extra 118 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-const c Real 2.5)(assert (> x c))")$assertions })), "[x > 2.5]")
})

test_that("smtlib: extra 118 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-const c Real 2.5)(assert (> x c))")$symbols) })), "{'c': 2.50000000000000}")
})

test_that("smtlib: extra 119 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-const c Int (+ 1 2))(declare-const c Int)(assert (= c 3))")$assertions })), "[Q.integer(c), Eq(c, 3)]")
})

test_that("smtlib: extra 119 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-const c Int (+ 1 2))(declare-const c Int)(assert (= c 3))")$symbols) })), "{'c': c}")
})

test_that("smtlib: extra 120 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-sort R () Real)(declare-const x R)")$assertions })), "[Q.real(x)]")
})

test_that("smtlib: extra 120 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-sort R () Real)(declare-const x R)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 121 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-sort I () Int)(declare-fun x () I)")$assertions })), "[Q.integer(x)]")
})

test_that("smtlib: extra 121 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-sort I () Int)(declare-fun x () I)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 122 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-sort A () B)(define-sort B () A)(declare-const x A)")$assertions })), "[]")
})

test_that("smtlib: extra 122 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-sort A () B)(define-sort B () A)(declare-const x A)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 123 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-sort S () (Array Int Int))(declare-const x S)")$assertions })), "[]")
})

test_that("smtlib: extra 123 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-sort S () (Array Int Int))(declare-const x S)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 124 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-sort P (T U) (Pair T U))")$assertions })), "ERR NotImplementedError: Parametric sort definitions (define-sort with parameters) are not supported")
})

test_that("smtlib: extra 125 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (let ((a 1)) (+ a x)))")$assertions })), "[x + 1]")
})

test_that("smtlib: extra 125 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (let ((a 1)) (+ a x)))")$symbols) })), "{}")
})

test_that("smtlib: extra 126 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (let ((a 1) (b 2)) (+ a b x)))")$assertions })), "[x + 3]")
})

test_that("smtlib: extra 126 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (let ((a 1) (b 2)) (+ a b x)))")$symbols) })), "{}")
})

test_that("smtlib: extra 127 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (let ((a 1)) (let ((a 2)) a)))")$assertions })), "[2]")
})

test_that("smtlib: extra 127 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (let ((a 1)) (let ((a 2)) a)))")$symbols) })), "{}")
})

test_that("smtlib: extra 128 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (let ((a x) (b a)) (+ a b)))")$assertions })), "[a + x]")
})

test_that("smtlib: extra 128 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (let ((a x) (b a)) (+ a b)))")$symbols) })), "{}")
})

test_that("smtlib: extra 129 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (let ((a 1) (a 2)) a))")$assertions })), "[2]")
})

test_that("smtlib: extra 129 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (let ((a 1) (a 2)) a))")$symbols) })), "{}")
})

test_that("smtlib: extra 130 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-const n Int 5)(assert (let ((n 1)) n))")$assertions })), "[5]")
})

test_that("smtlib: extra 130 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-const n Int 5)(assert (let ((n 1)) n))")$symbols) })), "{'n': 5}")
})

test_that("smtlib: extra 131 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const n Int)(assert (let ((n 1)) (> n 0)))")$assertions })), "[Q.integer(n), True]")
})

test_that("smtlib: extra 131 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const n Int)(assert (let ((n 1)) (> n 0)))")$symbols) })), "{'n': n}")
})

test_that("smtlib: extra 132 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (let ((p (> x 0))) (and p (< x 1))))")$assertions })), "[(x > 0) & (x < 1)]")
})

test_that("smtlib: extra 132 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (let ((p (> x 0))) (and p (< x 1))))")$symbols) })), "{}")
})

test_that("smtlib: extra 133 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (! (> x 0) :named a1))")$assertions })), "[x > 0]")
})

test_that("smtlib: extra 133 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (! (> x 0) :named a1))")$symbols) })), "{}")
})

test_that("smtlib: extra 134 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (! (> x 0) :named a1 :weight 3))")$assertions })), "[x > 0]")
})

test_that("smtlib: extra 134 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (! (> x 0) :named a1 :weight 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 135 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (! x :pattern ((f x))))")$assertions })), "[x]")
})

test_that("smtlib: extra 135 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (! x :pattern ((f x))))")$symbols) })), "{}")
})

test_that("smtlib: extra 136 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x (_ bv0 8)))")$assertions })), "[Eq(x, 0)]")
})

test_that("smtlib: extra 136 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x (_ bv0 8)))")$symbols) })), "{}")
})

test_that("smtlib: extra 137 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x (_ bv255 8)))")$assertions })), "[Eq(x, 255)]")
})

test_that("smtlib: extra 137 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x (_ bv255 8)))")$symbols) })), "{}")
})

test_that("smtlib: extra 138 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert ((_ bv5 32) x))")$assertions })), "ERR UnknownSMTLibOperatorError: Cannot apply arguments to expression 5")
})

test_that("smtlib: extra 139 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x (_ bv 8)))")$assertions })), "ERR NotImplementedError: Indexed identifier (_ bv ...) is not supported")
})

test_that("smtlib: extra 140 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x (_ bv1a 8)))")$assertions })), "ERR NotImplementedError: Indexed identifier (_ bv1a ...) is not supported")
})

test_that("smtlib: extra 141 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= ((_ zero_extend 4) x) 0))")$assertions })), "ERR NotImplementedError: Indexed identifier (_ zero_extend ...) is not supported")
})

test_that("smtlib: extra 142 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (as x Int) 1))")$assertions })), "[Eq(x, 1)]")
})

test_that("smtlib: extra 142 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= (as x Int) 1))")$symbols) })), "{}")
})

test_that("smtlib: extra 143 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= ((as f Int) 1) 2))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: f")
})

test_that("smtlib: extra 144 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-fun f (Int) Int)(assert (= ((as f Int) 1) 2))")$assertions })), "[Eq(f(1), 2)]")
})

test_that("smtlib: extra 144 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-fun f (Int) Int)(assert (= ((as f Int) 1) 2))")$symbols) })), "{}")
})

test_that("smtlib: extra 145 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (as (_ bv3 8) (_ BitVec 8)) x))")$assertions })), "[Eq(3, x)]")
})

test_that("smtlib: extra 145 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= (as (_ bv3 8) (_ BitVec 8)) x))")$symbols) })), "{}")
})

test_that("smtlib: extra 146 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= ((as (_ bv3 8) (_ BitVec 8)) 1) x))")$assertions })), "ERR UnknownSMTLibOperatorError: Cannot apply arguments to expression 3")
})

test_that("smtlib: extra 147 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(push 1)(declare-const x Int)(pop 1)")$assertions })), "[]")
})

test_that("smtlib: extra 147 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(push 1)(declare-const x Int)(pop 1)")$symbols) })), "{}")
})

test_that("smtlib: extra 148 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(push 2)(declare-const x Int)(pop 1)")$assertions })), "[]")
})

test_that("smtlib: extra 148 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(push 2)(declare-const x Int)(pop 1)")$symbols) })), "{}")
})

test_that("smtlib: extra 149 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(push 2)(declare-const x Int)(pop 2)")$assertions })), "[]")
})

test_that("smtlib: extra 149 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(push 2)(declare-const x Int)(pop 2)")$symbols) })), "{}")
})

test_that("smtlib: extra 150 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(push 1)(declare-const x Int)(pop 5)")$assertions })), "[]")
})

test_that("smtlib: extra 150 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(push 1)(declare-const x Int)(pop 5)")$symbols) })), "{}")
})

test_that("smtlib: extra 151 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(pop 1)")$assertions })), "[]")
})

test_that("smtlib: extra 151 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(pop 1)")$symbols) })), "{}")
})

test_that("smtlib: extra 152 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(push 0)(declare-const x Int)(pop 1)")$assertions })), "[Q.integer(x)]")
})

test_that("smtlib: extra 152 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(push 0)(declare-const x Int)(pop 1)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 153 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(push 1)(push 1)(assert a)(pop 1)(assert b)")$assertions })), "[b]")
})

test_that("smtlib: extra 153 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(push 1)(push 1)(assert a)(pop 1)(assert b)")$symbols) })), "{}")
})

test_that("smtlib: extra 154 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(push 1)(define-fun f ((a Int)) Int a)(pop 1)(assert (= (f 1) 1))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: f")
})

test_that("smtlib: extra 155 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(push 1)(declare-fun f (Int) Int)(pop 1)(assert (= (f 1) 1))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: f")
})

test_that("smtlib: extra 156 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(push 1)(define-sort R () Real)(pop 1)(declare-const x R)")$assertions })), "[]")
})

test_that("smtlib: extra 156 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(push 1)(define-sort R () Real)(pop 1)(declare-const x R)")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 157 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(push 1)(reset-assertions)(assert a)(pop 1)")$assertions })), "[]")
})

test_that("smtlib: extra 157 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(push 1)(reset-assertions)(assert a)(pop 1)")$symbols) })), "{}")
})

test_that("smtlib: extra 158 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert a)(push 1)(reset)(assert b)(pop 1)")$assertions })), "[b]")
})

test_that("smtlib: extra 158 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert a)(push 1)(reset)(assert b)(pop 1)")$symbols) })), "{}")
})

test_that("smtlib: extra 159 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Int)(reset)(declare-const y Real)")$assertions })), "[Q.real(y)]")
})

test_that("smtlib: extra 159 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const x Int)(reset)(declare-const y Real)")$symbols) })), "{'y': y}")
})

test_that("smtlib: extra 160 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert a)(reset-assertions)(assert b)")$assertions })), "[b]")
})

test_that("smtlib: extra 160 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert a)(reset-assertions)(assert b)")$symbols) })), "{}")
})

test_that("smtlib: extra 161 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(set-logic QF_LIA)(set-option :produce-models true)(set-info :status sat)")$assertions })), "[]")
})

test_that("smtlib: extra 161 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(set-logic QF_LIA)(set-option :produce-models true)(set-info :status sat)")$symbols) })), "{}")
})

test_that("smtlib: extra 162 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(set-info :source |a source|)")$assertions })), "[]")
})

test_that("smtlib: extra 162 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(set-info :source |a source|)")$symbols) })), "{}")
})

test_that("smtlib: extra 163 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(set-option :seed 3)")$assertions })), "[]")
})

test_that("smtlib: extra 163 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(set-option :seed 3)")$symbols) })), "{}")
})

test_that("smtlib: extra 164 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(set-info :x (a b c))")$assertions })), "[]")
})

test_that("smtlib: extra 164 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(set-info :x (a b c))")$symbols) })), "{}")
})

test_that("smtlib: extra 165 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(get-info :name)")$assertions })), "[]")
})

test_that("smtlib: extra 165 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(get-info :name)")$symbols) })), "{}")
})

test_that("smtlib: extra 166 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(get-option :seed)")$assertions })), "[]")
})

test_that("smtlib: extra 166 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(get-option :seed)")$symbols) })), "{}")
})

test_that("smtlib: extra 167 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(get-value (x y))")$assertions })), "[]")
})

test_that("smtlib: extra 167 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(get-value (x y))")$symbols) })), "{}")
})

test_that("smtlib: extra 168 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(get-value ((+ x 1)))")$assertions })), "[]")
})

test_that("smtlib: extra 168 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(get-value ((+ x 1)))")$symbols) })), "{}")
})

test_that("smtlib: extra 169 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(echo \"x\")")$assertions })), "[]")
})

test_that("smtlib: extra 169 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(echo \"x\")")$symbols) })), "{}")
})

test_that("smtlib: extra 170 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(check-sat-assuming ())")$assertions })), "[]")
})

test_that("smtlib: extra 170 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(check-sat-assuming ())")$symbols) })), "{}")
})

test_that("smtlib: extra 171 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(check-sat-assuming (a (not b) c))")$assertions })), "[]")
})

test_that("smtlib: extra 171 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(check-sat-assuming (a (not b) c))")$symbols) })), "{}")
})

test_that("smtlib: extra 172 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(get-assertions)(get-assignment)(get-proof)(get-unsat-assumptions)(get-unsat-core)")$assertions })), "[]")
})

test_that("smtlib: extra 172 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(get-assertions)(get-assignment)(get-proof)(get-unsat-assumptions)(get-unsat-core)")$symbols) })), "{}")
})

test_that("smtlib: extra 173 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-datatype P ((mk (fst Int) (snd Int))))")$assertions })), "ERR NotImplementedError: Algebraic datatypes (declare-datatype) are not supported")
})

test_that("smtlib: extra 174 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-datatype L (par (T) ((nil) (cons (hd T) (tl (L T))))))")$assertions })), "ERR NotImplementedError: Algebraic datatypes (declare-datatype) are not supported")
})

test_that("smtlib: extra 175 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-datatypes ((A 0) (B 0)) (((a1)) ((b1))))")$assertions })), "ERR NotImplementedError: Algebraic datatypes (declare-datatypes) are not supported")
})

test_that("smtlib: extra 176 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun-rec f ((x Int)) Int x)")$assertions })), "ERR NotImplementedError: Recursive macros are not supported in SymPy because AST generation evaluates bottom-up.")
})

test_that("smtlib: extra 177 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun-rec f () Int (bvadd 1 2))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvadd' has no SymPy equivalent and is not supported")
})

test_that("smtlib: extra 178 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (forall ((a Int)) (> a 0)))")$assertions })), "ERR NotImplementedError: Quantifiers (forall) are not supported since SymPy has no representation for them")
})

test_that("smtlib: extra 179 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (exists ((a Int) (b Int)) (> a b)))")$assertions })), "ERR NotImplementedError: Quantifiers (exists) are not supported since SymPy has no representation for them")
})

test_that("smtlib: extra 180 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (forall ((a Int)) (bvadd a 0)))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvadd' has no SymPy equivalent and is not supported")
})

test_that("smtlib: extra 181 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (match x ((a 1) ((cons h t) 2))))")$assertions })), "ERR NotImplementedError: Match expressions are not fully supported yet.")
})

test_that("smtlib: extra 182 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (match (unknown 1) ((a 1))))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: unknown")
})

test_that("smtlib: extra 183 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (unknown 1))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: unknown")
})

test_that("smtlib: extra 184 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (Not a))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: Not")
})

test_that("smtlib: extra 185 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (+ (unknown 1) (bvadd 1 2)))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: unknown")
})

test_that("smtlib: extra 186 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (bvadd (unknown 1) 2))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: unknown")
})

test_that("smtlib: extra 187 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (|and| a b))")$assertions })), "[a & b]")
})

test_that("smtlib: extra 187 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (|and| a b))")$symbols) })), "{}")
})

test_that("smtlib: extra 188 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (|x y| 1))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: x y")
})

test_that("smtlib: extra 189 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (-5 1))")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: -5")
})

test_that("smtlib: extra 190 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(foo)")$assertions })), "ERR UnknownSMTLibCommandError: Unknown SMT-LIB command: foo")
})

test_that("smtlib: extra 191 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(foo (bar 1) :k \"s\" #x1)")$assertions })), "ERR UnknownSMTLibCommandError: Unknown SMT-LIB command: foo")
})

test_that("smtlib: extra 192 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(|assert| (= x 1))")$assertions })), "ERR UnknownSMTLibCommandError: Unknown SMT-LIB command: assert")
})

test_that("smtlib: extra 193 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x 1))(foo)(assert (unknown 1))")$assertions })), "ERR UnknownSMTLibCommandError: Unknown SMT-LIB command: foo")
})

test_that("smtlib: extra 194 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (unknown 1))(foo)")$assertions })), "ERR UnknownSMTLibOperatorError: Unknown operator or function: unknown")
})

test_that("smtlib: extra 195 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const x Int)(asert x)")$assertions })), "ERR UnknownSMTLibCommandError: Unknown SMT-LIB command: asert")
})

test_that("smtlib: extra 196 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= |x| x))")$assertions })), "[True]")
})

test_that("smtlib: extra 196 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= |x| x))")$symbols) })), "{}")
})

test_that("smtlib: extra 197 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const |x| Real)(assert (> x 0))")$assertions })), "[Q.real(x), x > 0]")
})

test_that("smtlib: extra 197 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const |x| Real)(assert (> x 0))")$symbols) })), "{'x': x}")
})

test_that("smtlib: extra 198 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= True true))")$assertions })), "[Eq(True, True)]")
})

test_that("smtlib: extra 198 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= True true))")$symbols) })), "{}")
})

test_that("smtlib: extra 199 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= False false))")$assertions })), "[Eq(False, False)]")
})

test_that("smtlib: extra 199 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= False false))")$symbols) })), "{}")
})

test_that("smtlib: extra 200 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert TRUE)")$assertions })), "[TRUE]")
})

test_that("smtlib: extra 200 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert TRUE)")$symbols) })), "{}")
})

test_that("smtlib: extra 201 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-const true Int)(assert (= true 1))")$assertions })), "[Q.integer(true), Eq(true, 1)]")
})

test_that("smtlib: extra 201 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-const true Int)(assert (= true 1))")$symbols) })), "{'true': true}")
})

test_that("smtlib: extra 202 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun true () Bool false)(assert true)")$assertions })), "[False]")
})

test_that("smtlib: extra 202 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun true () Bool false)(assert true)")$symbols) })), "{'true': False}")
})

test_that("smtlib: extra 203 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(define-fun and ((a Bool) (b Bool)) Bool (or a b))(assert (and x y))")$assertions })), "[x | y]")
})

test_that("smtlib: extra 203 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(define-fun and ((a Bool) (b Bool)) Bool (or a b))(assert (and x y))")$symbols) })), "{}")
})

test_that("smtlib: extra 204 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(declare-fun + (Int Int) Int)(assert (= (+ 1 2) 3))")$assertions })), "[Eq(+(1, 2), 3)]")
})

test_that("smtlib: extra 204 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(declare-fun + (Int Int) Int)(assert (= (+ 1 2) 3))")$symbols) })), "{}")
})

test_that("smtlib: extra 205 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x 1)) ; trailing comment")$assertions })), "[Eq(x, 1)]")
})

test_that("smtlib: extra 205 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x 1)) ; trailing comment")$symbols) })), "{}")
})

test_that("smtlib: extra 206 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert ; inner comment\n (= x 1))")$assertions })), "[Eq(x, 1)]")
})

test_that("smtlib: extra 206 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert ; inner comment\n (= x 1))")$symbols) })), "{}")
})

test_that("smtlib: extra 207 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert\t(=\nx\r\n1))")$assertions })), "[Eq(x, 1)]")
})

test_that("smtlib: extra 207 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert\t(=\nx\r\n1))")$symbols) })), "{}")
})

test_that("smtlib: extra 208 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 209 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib(")")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 210 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 211 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (= x 1)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 212 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (= x 1)))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 213 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert ())")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 214 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (5 x))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 215 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(push)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 216 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(pop)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 217 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(push x)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 218 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-sort B 1 2)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 219 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-sort B)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 220 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-const x)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 221 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-const 5 Int)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 222 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-fun f Int Int)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 223 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (= x \"unclosed))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 224 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-fun |unclosed () S)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 225 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (= x #xG))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 226 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (= x #b2))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 227 assertions matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= x 1.))")$assertions })), "[Eq(1, .) & Eq(x, 1)]")
})

test_that("smtlib: extra 227 symbols matches SymPy", {
  expect_identical(.pyrun(local({ .pyD(parse_smtlib("(assert (= x 1.))")$symbols) })), "{}")
})

test_that("smtlib: extra 228 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(check-sat 1)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 229 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(exit now)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 230 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(set-logic)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 231 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (let () x))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 232 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (let (x 1) x))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 233 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (forall () x))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 234 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (! x))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 235 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (_ bv5))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 236 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (= x ,))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 237 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("assert")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 238 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert x) junk")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 239 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (= x 1)) )")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 240 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert \"a\"\")")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 241 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(declare-datatype C ())")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 242 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(get-info x)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 243 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(echo x)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 244 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(check-sat-assuming ((and a)))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 245 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(assert (as x))")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: extra 246 syntax matches SymPy", {
  expect_identical(.pyrun(local({ .pyEC(parse_smtlib("(:k)")) })), "ERR SMTLibSyntaxError")
})

test_that("smtlib: unsupported bv2int matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bv2int x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bv2int' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bv2nat matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bv2nat x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bv2nat' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvadd matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvadd x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvadd' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvand matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvand x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvand' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvashr matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvashr x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvashr' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvcomp matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvcomp x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvcomp' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvlshr matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvlshr x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvlshr' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvmul matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvmul x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvmul' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvnand matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvnand x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvnand' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvneg matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvneg x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvneg' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvnor matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvnor x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvnor' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvnot matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvnot x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvnot' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvor matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvor x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvor' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvsdiv matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvsdiv x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvsdiv' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvsge matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvsge x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvsge' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvsgt matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvsgt x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvsgt' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvshl matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvshl x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvshl' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvsle matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvsle x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvsle' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvslt matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvslt x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvslt' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvsmod matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvsmod x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvsmod' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvsrem matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvsrem x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvsrem' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvsub matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvsub x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvsub' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvudiv matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvudiv x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvudiv' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvuge matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvuge x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvuge' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvugt matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvugt x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvugt' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvule matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvule x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvule' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvult matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvult x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvult' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvurem matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvurem x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvurem' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvxnor matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvxnor x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvxnor' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported bvxor matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (bvxor x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'bvxor' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported concat matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (concat x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'concat' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported const matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (const x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'const' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported ext_rotate_left matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (ext_rotate_left x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'ext_rotate_left' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported ext_rotate_right matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (ext_rotate_right x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'ext_rotate_right' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported int.to.str matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (int.to.str x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'int.to.str' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported select matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (select x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'select' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported store matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (store x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'store' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported str.++ matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (str.++ x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'str.++' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported str.at matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (str.at x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'str.at' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported str.contains matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (str.contains x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'str.contains' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported str.indexof matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (str.indexof x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'str.indexof' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported str.len matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (str.len x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'str.len' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported str.prefixof matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (str.prefixof x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'str.prefixof' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported str.replace matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (str.replace x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'str.replace' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported str.substr matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (str.substr x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'str.substr' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported str.suffixof matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (str.suffixof x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'str.suffixof' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported str.to.int matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (str.to.int x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'str.to.int' has no SymPy equivalent and is not supported")
})

test_that("smtlib: unsupported ubv_to_int matches SymPy", {
  expect_identical(.pyrun(local({ parse_smtlib("(assert (= (ubv_to_int x y) 0))")$assertions })), "ERR NotImplementedError: The SMT-LIB operator 'ubv_to_int' has no SymPy equivalent and is not supported")
})

