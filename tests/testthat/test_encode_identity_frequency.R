context("encode_identity_frequency")

test_that("encode_cell_identity_frequency_matrix", {
  toy_data = readRDS("toy_data.rds")
  result = encode_cell_identity_frequency_matrix(
    list(one = toy_data, two = toy_data),
    attributes = "visit",
    group_by = "library"
  )
  expect_equal(result["sample_12080431_one", "Baseline"], 12)
  expect_equal(result["sample_12080431_one", "Baseline"], result["sample_12080431_two", "Baseline"])
})

