source(file.path("..", "..", "R", "fonctions.R"))

test_that("CKD-EPI 2021 reproduit une valeur calculée à la main", {
  # Femme, 60 ans, 0,7 mg/dL : 142 x 1 x 1 x 0,9938^60 x 1,012
  expect_equal(ckdepi2021(0.7, 60, "F"), 142 * 0.9938^60 * 1.012)
  # Homme, 60 ans, 1,8 mg/dL : 142 x 2^-1,2 x 0,9938^60
  expect_equal(ckdepi2021(1.8, 60, "M"), 142 * 2^-1.2 * 0.9938^60)
  expect_lt(ckdepi2021(2, 70, "M"), ckdepi2021(1, 70, "M"))
})

test_that("les formules de Schoenfeld sont réciproques", {
  d <- schoenfeld_events(1.3, sigma = 1)
  expect_equal(schoenfeld_hr_min(d, sigma = 1), 1.3)
  expect_equal(d, (1.959964 + 1.281552)^2 / log(1.3)^2, tolerance = 1e-5)
  expect_gt(schoenfeld_events(1.3, 1, r2 = .3), d)
})

test_that("Rubin se réduit au cas simple sans variance inter-imputations", {
  r <- suppressWarnings(rubin(rep(.5, 5) + c(1e-9, -1e-9, 0, 0, 0), rep(.04, 5)))
  expect_equal(r$est, .5); expect_equal(r$se, .2, tolerance = 1e-6)
})

test_that("Rubin : variance totale w + (1 + 1/m) b", {
  r <- rubin(c(1, 2, 3), c(1, 1, 1))
  expect_equal(r$se^2, 1 + (4 / 3) * 1)
})

test_that("E-value : croît avec le HR et vaut 1 pour HR = 1", {
  expect_equal(evalue_hr(1), 1)
  expect_gt(evalue_hr(2), evalue_hr(1.5))
  # HR = 2 -> RR approché (1 - 0,5^sqrt(2)) / (1 - 0,5^sqrt(1/2)) = 1,612...
  rr <- (1 - 0.5^sqrt(2)) / (1 - 0.5^sqrt(0.5))
  expect_equal(evalue_hr(2), rr + sqrt(rr * (rr - 1)))
  # Un HR protecteur donne la même E-value que son inverse
  expect_equal(evalue_hr(0.5), evalue_hr(2))
})

test_that("intervalle de Rubin cohérent avec l'estimation", {
  h <- hr_ic(rubin(log(c(1.5, 1.6, 1.7)), c(.01, .01, .01)))
  expect_lt(h$bas, h$hr); expect_gt(h$haut, h$hr)
  expect_equal(h$hr, exp(mean(log(c(1.5, 1.6, 1.7)))))
})

test_that("formatage français", {
  expect_equal(fmt_hr(1.234, 1.1, 1.4), "1,23 (1,10 à 1,40)")
  expect_equal(fmt_p(0.0004), "< 0,001"); expect_equal(fmt_p(0.0512), "0,051")
})
