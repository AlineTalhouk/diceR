
context("Graphical displays")

set.seed(911)
x <- matrix(rnorm(1000), nrow = 10)
CC1 <- ConClust(x, k = 4, reps = 10, method = c("hcAEucl", "apEucl", "gmmBIC"),
                save = FALSE)
p1 <- graph_cdf(CC1)

test_that("graph_cdf object can have added/modified ggplot layers", {
  p2 <- p1 +
    labs(y = "Probability") +
    stat_ecdf(aes(colour = Method)) +
    scale_color_brewer(palette = "Set2") +
    theme(legend.position = "none")
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
  expect_false(isTRUE(all.equal(p1, p2)))
})

test_that("graph_heatmap can have same plot but different titles", {
  phm1 <- graph_heatmap(CC1)
  phm2 <- graph_heatmap(CC1, main = c("A", "B", "C"))
  expect_identical(phm1, phm2)
})

test_that("error in graph_heatmap if too few titles", {
  expect_error(graph_heatmap(CC1, main = "A"))
})
