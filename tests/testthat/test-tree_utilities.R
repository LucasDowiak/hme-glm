sym_tree <- c(
  "0",
  "0.1", "0.2",
  "0.1.1", "0.1.2", "0.2.1", "0.2.2"
)

sym2_tree <- c(
  "0",
  "0.1", "0.2", "0.3",
  "0.2.1", "0.2.2"
)

asy_tree <- c(
  "0",
  "0.1", "0.2",
  "0.1.1", "0.1.2",
  "0.1.1.1", "0.1.1.2"
)

test_that("children() works", {
  expect_equal(children(as.list(sym2_tree), sym2_tree),
               list('0'=c("0.1", "0.2", "0.3"), '0.1'=NULL, '0.2'=c("0.2.1", "0.2.2"),
                    '0.3'=NULL, '0.2.1'=NULL, '0.2.2'=NULL))
})

test_that("parent() works", {
  expect_equal(parent(sym_tree, sym_tree),
               list('0'=NULL, '0.1'="0", '0.2'="0", '0.1.1'="0.1", '0.1.2'="0.1",
                    '0.2.1'="0.2", '0.2.2'="0.2"))
})

test_that("ancestors() works", {
  expect_equal(ancestors(asy_tree),
               list('0'="0", '0.1'=c("0", "0.1"), '0.2'=c("0", "0.2"),
                    '0.1.1'=c("0", "0.1", "0.1.1"), '0.1.2'=c("0", "0.1", "0.1.2"),
                    '0.1.1.1'=c("0", "0.1", "0.1.1", "0.1.1.1"),
                    '0.1.1.2'=c("0", "0.1", "0.1.1", "0.1.1.2")))
})

test_that("progeny() works", {
  expect_equal(progeny(sym_tree, sym_tree),
               list('0'=c("0.1", "0.1.1", "0.1.2", "0.2", "0.2.1", "0.2.2"),
                    '0.1'=c("0.1.1", "0.1.2"), '0.2'=c("0.2.1", "0.2.2"),
                    '0.1.1'=NULL, '0.1.2'=NULL, '0.2.1'=NULL, '0.2.2'=NULL))
  
  expect_equal(progeny(sym2_tree, sym2_tree), 
               list('0'=c("0.1", "0.2", "0.2.1", "0.2.2", "0.3"),
                    '0.1'=NULL, '0.2'=c("0.2.1", "0.2.2"), '0.3'=NULL,
                    '0.2.1'=NULL, '0.2.2'=NULL))
})

test_that("is_terminal() works", {
  expect_equal(is_terminal(sym_tree, sym_tree),
               list('0'=FALSE, '0.1'=FALSE, '0.2'=FALSE,
                    '0.1.1'=TRUE, '0.1.2'=TRUE, '0.2.1'=TRUE, '0.2.2'=TRUE))
  
  expect_equal(is_terminal(asy_tree, asy_tree),
               list('0'=FALSE, '0.1'=FALSE, '0.2'=TRUE,
                    '0.1.1'=FALSE, '0.1.2'=TRUE, '0.1.1.1'=TRUE, '0.1.1.2'=TRUE))
})

test_that("last_split() works", {
  expect_equal(last_split(asy_tree), list('0'=NULL, '0.1'="1", '0.2'="2", '0.1.1'="1",
                                       '0.1.2'="2", '0.1.1.1'="1", '0.1.1.2'="2"))
})

test_that("all_but_last_split() works", {
  expect_equal(all_but_last_split(asy_tree),
               list('0'=NULL, '0.1'="0", '0.2'="0", '0.1.1'="0.1",
                    '0.1.2'="0.1", '0.1.1.1'="0.1.1", '0.1.1.2'="0.1.1"))
})


test_that("expert_index() works", {
  expect_equal(expert_index(sym2_tree), c('0'=0, '0.1'=1, '0.2'=0, '0.3'=2, '0.2.1'=3, '0.2.2'=4))
})


