# Siler function parameters fit to Coale & Demeny West model life tables for females.
# Each regime is a one-row data frame with columns: a1, b1, a2, a3, b3, name

# Fast mortality
CoaleDemenyWestF3 <- structure(data.frame(
  a1 = 0.558,
  b1 = 1.05,
  a2 = 0.01225,
  a3 = 0.000520,
  b3 = 0.0727,
  name = "CoaleDemenyWestF3"),
  class = c("persephone_siler_regime", "data.frame")
)

CoaleDemenyWestF5 <- structure(data.frame(
  a1 = 0.457,
  b1 = 1.07,
  a2 = 0.01037,
  a3 = 0.000359,
  b3 = 0.0763,
  name = "CoaleDemenyWestF5"),
  class = c("persephone_siler_regime", "data.frame")
)

CoaleDemenyWestF11 <- structure(data.frame(
  a1 = 0.256,
  b1 = 1.17,
  a2 = 0.00596,
  a3 = 0.000133,
  b3 = 0.086,
  name = "CoaleDemenyWestF11"),
  class = c("persephone_siler_regime", "data.frame")
)

CoaleDemenyWestF15 <- structure(data.frame(
  a1 = 0.175,
  b1 = 1.40,
  a2 = 0.00368,
  a3 = 0.000075,
  b3 = 0.0917,
  name = "CoaleDemenyWestF15"),
  class = c("persephone_siler_regime", "data.frame")
)

CoaleDemenyWestF17 <- structure(data.frame(
  a1 = 0.14,
  b1 = 1.57,
  a2 = 0.00265,
  a3 = 0.000056,
  b3 = 0.0949,
  name = "CoaleDemenyWestF17"),
  class = c("persephone_siler_regime", "data.frame")
)

# Slow mortality
CoaleDemenyWestF21 <- structure(data.frame(
  a1 = 0.091,
  b1 = 2.78,
  a2 = 0.00092,
  a3 = 0.000025,
  b3 = 0.1033,
  name = "CoaleDemenyWestF21"),
  class = c("persephone_siler_regime", "data.frame")
)

usethis::use_data(
  CoaleDemenyWestF3,
  CoaleDemenyWestF5,
  CoaleDemenyWestF11,
  CoaleDemenyWestF15,
  CoaleDemenyWestF17,
  CoaleDemenyWestF21,
  overwrite = TRUE
)
