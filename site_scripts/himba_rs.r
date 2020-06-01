
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "himba_rs"

d0 <- read.csv("./data/himba_rs.csv", stringsAsFactors = FALSE)

#simply not enough known father RS, so I'm going to use mothers only.

#just the moms
dm <- data.frame(d0$rsm, d0$agem/10, d0$mid)
names(dm) <- c("rs", "age", "id")
drop <- which(duplicated(dm$id) | is.na(dm$age * dm$rs)) 
dm <- dm[-drop, ]
#85 moms
D_m <- cbind(rep(1, nrow(dm)), dm$age)

#now all kids. none are in the parents column: good!
drop <- which(d0$age<18 | is.na(d0$age * d0$rs) | duplicated(d0$id))
do <- d0[-drop, ]
do$age <- do$age/10
D_o <- cbind(do$male, (1-do$male), do$male * do$age,
  (1-do$male) * do$age, do$male * do$age^2, (1-do$male) * do$age^2)

#now just folks with parents
drop <- which(d0$age < 18 |
  is.na(d0$age * d0$rs * d0$agem * d0$rsm) |
  duplicated(d0$id))
d1 <- d0[-drop, ]
d1$age <- d1$age/10
d1$agem <- d1$agem/10
D_m1 <- cbind(rep(1, nrow(d1)), d1$agem)
D_o1 <- cbind(d1$male, (1 - d1$male), d1$male * d1$age,
  (1 - d1$male) * d1$age, d1$male * d1$age^2, (1 - d1$male) * d1$age^2)

(mean(d1$rs[which(d1$male == 1)])-mean(d1$rs[which(d1$male == 0)]))/sd(d1$rs)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  nm = nrow(dm), 
  no = nrow(do), 
  n1 = nrow(d1), 
  rs_m = dm$rs, 
  rs_o = do$rs, 
  rs_m1 = d1$rsm, 
  rs_o1 = d1$rs, 
  D_m = D_m, 
  D_o = D_o, 
  D_m1 = D_m1, 
  D_o1 = D_o1
)

expect_true(data_list$nm == 85)
expect_true(data_list$no == 157)
expect_true(data_list$n1 == 146)



tic(paste("fit", measure, "stan model"))

init <- list(
  list(
    P_mu_p = c(7.58, -0.45), 
    P_sigma_p = 2.28,
    P_mu_o = c(-9.53, -1.24, 4.69, 1.11, -0.49, -0.12),
    P_sigma_o = c(-3.41, -2.52, 1.62, 1.45, -0.12, -0.15)
  )
)

if (!file.exists(file.path("temp", paste0(measure, ".stanfit")))){

  m1 <- stan(file = file.path("site_models", paste0(measure, ".stan")), data = data_list, 
    iter = n_iter, warmup = n_warmup, chains = 1, init = init)

  save(m1, file = file.path("temp", paste0(measure, ".stanfit")))

}
toc(log = TRUE)



tic("calculate himba_rs beta coefficients")

load(file.path("temp", paste0(measure, ".stanfit")))

post <- extract(m1, permute=TRUE)

n_sim <- min(10000, length(post$lp__))

beta_s <- rep(NA, n_sim)
beta_d <- rep(NA, n_sim)

for (i in 1:n_sim)
{
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
 
  if (i %% 100 == 0) print(i)

}

toc()

# reported values here
beta_s_est <- 0.145
beta_s_se  <- 0.132
beta_d_est <- 0.187
beta_d_se  <- 0.111

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.042)
beta_diff_se <- 0.171

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
