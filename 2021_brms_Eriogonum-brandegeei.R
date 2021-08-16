library(brms)


fit1 <- brm(formula = surv ~ RosNew + TempFall + (1 + TransectNew.num|Year.num),
                                                           data = dats, family = lognormal(),
                                                           prior = c(set_prior("normal(0,5)", class = "b"),
                                                                     set_prior("cauchy(0,2)", class = "sd"),
                                                                     set_prior("lkj(2)", class = "cor")),
                                                           warmup = 1000, iter = 2000, chains = 4,
                                                           control = list(adapt_delta = 0.95))
