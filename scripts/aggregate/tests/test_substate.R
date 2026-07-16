# Validates the substate composition is wired into the master and the
# myofibroblast_expansion hypothesis is no longer dormant.
library(testthat); library(data.table)
m <- fread("results/aggregate/results_master.tsv")

test_that("substate_composition reaches the master with activated-fibroblast rows", {
  expect_true("substate_composition" %in% m$readout_class)
  expect_true("Fibroblast_activated" %in% m[readout_class == "substate_composition", unit])
})

test_that("myofibroblast_expansion is claimed and confirmatory", {
  h <- m[hypothesis == "myofibroblast_expansion"]
  expect_gt(nrow(h), 0)
  expect_true(all(h$readout_class == "substate_composition"))
  expect_true(all(h$tier == "confirmatory"))
})
