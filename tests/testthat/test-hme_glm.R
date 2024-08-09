test_that("Iris Two-Expert MoE: expect S3 object",
  {
    data(iris)
    iris$color <- with(iris, ifelse(Species == "setosa", "blue",
                                            ifelse(Species == "versicolor", "orange", "green")))
    iris$points <- with(iris, ifelse(Species == "setosa", 1,
                                             ifelse(Species == "versicolor", 2, 3)))
     
    tree <- c("0", "0.1", "0.2")
    mod <- suppressWarnings(hme(tree, "Sepal.Width ~ Petal.Width | Petal.Width",
                                data=iris, maxiter=200, tolerance = 1e-6, trace=1))
    expect_s3_class(mod, "hme")
    expect_length(mod[["vcv"]][["OPG"]], 36)
    expect_length(mod[["gate.pars"]], 1)
    expect_length(mod[["expert.pars"]], 2)
  }
)