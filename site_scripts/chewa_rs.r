
rm(list = ls())

source("./project_support.r")

measure <- "chewa_rs"

d0 <- read.csv("./data/chewa_rs.csv", stringsAsFactors = FALSE)

d0$male[which(d0$male == 2)] <- 0
drop <- which(d0$age < 18)
d0 <- d0[-drop, ]

# gonna use mothers only, since we have A LOT more data for them.
# Fathers USUALLY have the same RS, so imputation ought to be possible.

#just the moms
dp <- data.frame(d0$rsm, d0$agem/10, d0$idm)
names(dp) <- c("rs", "age", "id")
dp$id <- as.character(dp$id)
drop <- which(
  duplicated(dp$id) |
  is.na(dp$age * dp$rs) |
  dp$age < 0 |
  dp$rs == 0
)
# looks like some parents have negative age and 0 rs...?
dp <- dp[-drop, ]
#298 moms
#no predictors here

#just the kids
drop <- which(
  duplicated(d0$id) |
  is.na(d0$age * d0$rs) |
  d0$age < 0 |
  d0$id %in% d0$idf |
  d0$id %in% d0$idm
)
do <- d0[-drop, ]
do$age <- do$age/10
#1711 kids
D_o <- cbind(
  1 - do$male,
  (1 - do$male) * do$age,
  do$male,
  do$male * do$age
)

#folks with parents
drop <- which(
  duplicated(d0$id) |
  is.na(d0$age * d0$rs * d0$agem * d0$rsm) |
  d0$age < 0
)
d1 <- d0[-drop, ]
d1$age <- d1$age/10
d1$agem <- d1$agem/10
D_o1 <- cbind(
  1 - d1$male,
  (1 - d1$male) * d1$age,
  d1$male,
  d1$male * d1$age
)

(mean(d1$rs[which(d1$male == 1)])-mean(d1$rs[which(d1$male == 0)]))/sd(d1$rs)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  np = nrow(dp), 
  no = nrow(do), 
  n1 = nrow(d1), 
  rsp = dp$rs, 
  rso = do$rs, 
  rso1 = d1$rs, 
  rsm = d1$rsm, 
  D_o = D_o, 
  D_o1 = D_o1, 
  male_o = do$male, 
  male_o1 = d1$male
)

expect_true(data_list$np == 296)
expect_true(data_list$no == 1711)
expect_true(data_list$n1 == 167)



tic(paste("fit", measure, "stan model"))

if (!file.exists(file.path("temp", paste0(measure, ".stanfit")))){

  m1 <- stan(file = file.path("site_models", paste0(measure, ".stan")), data = data_list, 
    iter = n_iter, warmup = n_warmup, chains = 1, init = init)

  save(m1, file = file.path("temp", paste0(measure, ".stanfit")))

}
toc(log = TRUE)



load(file.path("temp", paste0(measure, ".stanfit")))

post <- extract(m1, permute=TRUE)

n_sim <- min(10000, length(post$lp__))

beta_s <- rep(NA, n_sim)
beta_d <- rep(NA, n_sim)

for (i in 1:n_sim) {

  #do the regression
  par_temp <- post$Dev_m[i, ]
  dtemp <- data.frame(post$Dev_o[i, ], par_temp)
  names(dtemp) <- c("Dev_o", "Dev_m")
  
  #sons
  dtemps <- dtemp[which(d1$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * Dev_m, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(d1$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * Dev_m, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]

}

# reported values here
beta_s_est <- (-0.41)
beta_s_se  <- 0.16
beta_d_est <- 0.082
beta_d_se  <- 0.091

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.492)
beta_diff_se <- 0.186

beta_diff <- beta_s - beta_d
expect_true(abs(mean(beta_diff) - beta_diff_est) < beta_tol)
expect_true(abs(sd(beta_diff) - beta_diff_se) < beta_tol)


report <- list(
  measure = measure,
  n_sons = sum(d1$male == 1),
  n_daughters = sum(d1$male == 0),
  beta_son_est = sprintf("%1.3f", mean(beta_s)),
  beta_son_se  = sprintf("%1.3f", sd(beta_s)),
  beta_dau_est = sprintf("%1.3f", mean(beta_d)),
  beta_dau_se  = sprintf("%1.3f", sd(beta_d)),
  beta_diff_est = sprintf("%1.3f", mean(beta_diff)),
  beta_diff_se = sprintf("%1.3f", sd(beta_diff))
)

write_json(report, path = paste0("./temp/", measure, ".json"), pretty = TRUE)
