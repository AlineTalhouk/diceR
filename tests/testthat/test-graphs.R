
context("Graphical displays")

set.seed(911)
x <- matrix(rnorm(1000), nrow = 10)
CC1 <- ConClust(x, k = 4, reps = 10, method = c("hcAEucl", "apEucl", "gmmBIC"),
                save = FALSE)
CC1.summ <- consensus_summary(CC1, k = 4)
y <- consensus_combine(CC1.summ, element = "matrix")
p1 <- graph_cdf(y)

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
