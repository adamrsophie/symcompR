# Building a noncommutative object must switch on the noncommutative paths of
# Mul even when no constructor sets the state itself, so factor order is kept.

with_nc_off <- function(code) {
  nc <- symcompR:::.nc_state; qm <- symcompR:::.qm_state
  old <- c(nc$any, qm$on)
  on.exit({ nc$any <- old[[1L]] || nc$any; qm$on <- old[[2L]] || qm$on })
  nc$any <- FALSE; qm$on <- FALSE
  code
}

test_that("second quantization operators keep their order", {
  with_nc_off({
    a <- B(0); ad <- Bd(0)
    expect_true(symcompR:::.nc_state$any)
    expect_true(symcompR:::.qm_state$on)
    expect_false(identical(sstr(a * ad), sstr(ad * a)))
  })
})

test_that("quantum operators keep their order", {
  with_nc_off({
    A <- Operator("A"); Bo <- Operator("B")
    expect_equal(sstr(Bo * A), "B*A")
    expect_equal(sstr(A * Bo), "A*B")
  })
})

test_that("UnevaluatedExpr and vector fields switch the state on", {
  with_nc_off({
    UnevaluatedExpr(Symbol("x"))
    expect_true(symcompR:::.nc_state$any)
  })
  with_nc_off({
    rn <- diffgeom_rn()
    e_x <- base_vector(rn$R2_r, 0); e_y <- base_vector(rn$R2_r, 1)
    expect_true(symcompR:::.nc_state$any)
    expect_false(identical(sstr(e_x * e_y), sstr(e_y * e_x)))
  })
})

test_that("noncommutative symbols keep their order", {
  with_nc_off({
    A <- Symbol("A", commutative = FALSE); Bs <- Symbol("B", commutative = FALSE)
    expect_equal(sstr(Bs * A), "B*A")
  })
})
