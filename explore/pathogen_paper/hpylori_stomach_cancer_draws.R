options(stringsAsFactors = F)

# per https://stats.stackexchange.com/questions/526416/how-to-get-beta-prior-distribution-based-on-given-ci
# using optimization, simulate a candidate beta distribution the matches the CI and mean and then sample 1k draws

entity<-'hpylori_stomachcancer'

# fill in Attributable Fractions (AFs) from here: https://ars.els-cdn.com/content/image/1-s2.0-S2214109X19304887-mmc1.pdf

mu <- 0.89
ci <- c(0.79,0.94)

get.beta <- function(alpha) {
  alpha * (1-mu)/mu
}

quantiles <- c(0.025, 0.975)

cost <- function(alpha) {
  qs <- qbeta(quantiles, alpha, get.beta(alpha))
  sum(abs(qs-ci))
}

opt <- optim(1, cost, method="BFGS")
estimated.alpha <- opt$par
estimated.beta  <- get.beta(estimated.alpha)
estimated.ci    <- qbeta(quantiles, estimated.alpha, estimated.beta)

draws <- as.data.frame(matrix(nrow = 1000, ncol = 2))
colnames(draws) <- c('pathogen_neoplasm', 'draw')
draws$pathogen_neoplasm <- entity
draws$draw <- rbeta(1000, estimated.alpha, estimated.beta)

write.csv(draws, file = paste0('FILEPATH', entity, '_draws.csv'), row.names = F)

hist(rbeta(1000, estimated.alpha, estimated.beta))