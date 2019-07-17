
rm(list = ls())

source("./project_support.r")

measure <- "maya_rs"

d0 <- read.csv("./data/maya_rs.csv", stringsAsFactors = FALSE)

# n_iter <- 10000

names(d0)[which(names(d0) == "sex")] <- "male"

d0$male[which(d0$male == 2)] <- 0
d0$age <- as.numeric(paste(d0$age))

#separately analyze the offspring dataset and the parent dataset

#just the kids
rsf <- rep(NA, nrow(d0))
rsm <- rep(NA, nrow(d0))
agef <- rep(NA, nrow(d0))
agem <- rep(NA, nrow(d0))
for (i in 1:nrow(d0)) {
  if (length(which(d0$id == d0$fid[i])) > 0) {
    rsf[i] <- d0$rs[which(d0$id == d0$fid[i])]
    agef[i] <- d0$age[which(d0$id == d0$fid[i])]
  }
  if (length(which(d0$id == d0$mid[i])) > 0) {
    rsm[i] <- d0$rs[which(d0$id == d0$mid[i])]
    agem[i] <- d0$age[which(d0$id == d0$mid[i])]
  }
}
do <- data.frame(d0, rsf, rsm, agef/10, agem/10)
do$age <- do$age/10
names(do) <- c(names(d0), "rsf", "rsm", "agef", "agem")
#perfect monogamy, so imputation is easy.
drop <- which(is.na(do$rsf) & is.na(do$rsm))
do <- do[-drop, ]
for (i in 1:nrow(do)) {
  if (is.na(do$rsf[i])) do$rsf[i] <- do$rsm[i]
  else if(is.na(do$rsm[i])) do$rsm[i] <- do$rsf[i]
}
drop <- which(is.na(do$rs * do$age * do$agem * do$agef))
do <- do[-drop, ]
D_o <- cbind(rep(1, nrow(do)), do$age)


#now the parents, generally, for their own model
dp <- data.frame(
  c(do$fid, do$mid),
  c(do$rsf, do$rsm),
  c(do$agef, do$agem),
  c(rep(1, nrow(do)),
    rep(0, nrow(do))))
names(dp) <- c("id", "rs", "age", "male")
drop <- which(duplicated(dp$id))
dp <- dp[-drop, ]


#make list for stan
data_list <- list(
  no = nrow(do), 
  np = nrow(dp), 
  RS_o = do$rs, 
  RS_f = do$rsf, 
  RS_m = do$rsm, 
  RS_p = dp$rs, 
  D_o = D_o
)

expect_true(data_list$no == 187)
expect_true(data_list$np == 92)



tic(paste("fit", measure, "stan model"))

init <- list(
  list(
    P_mu = c(-10, 3, -10, 3),
    P_sig = c(-9, 3.5),
    max_RS_o = 5,
    max_sig_o = 3,
    mu_p = 7,
    sigma_p = 2
  )
)

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

for(i in 1:n_sim)
{
  #parents are perfectly correlated because they have the same values AND age model!
  
  #now do the regressions...
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])
  dtemp <- data.frame(post$Dev_o[i, ], midpar_temp)
  names(dtemp) <- c("Dev_o", "midpar_temp")
  
  #sons
  dtemps <- dtemp[which(do$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(do$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]
}

# reported values here
beta_s_est <- 0.067
beta_s_se  <- 0.12
beta_d_est <- (-0.005)
beta_d_se  <- 0.109

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.072
beta_diff_se <- 0.16

beta_diff <- beta_s - beta_d
expect_true(abs(mean(beta_diff) - beta_diff_est) < beta_tol)
expect_true(abs(sd(beta_diff) - beta_diff_se) < beta_tol)


report <- list(
  measure = measure,
  n_sons = sum(do$male == 1),
  n_daughters = sum(do$male == 0),
  beta_son_est = sprintf("%1.3f", mean(beta_s)),
  beta_son_se  = sprintf("%1.3f", sd(beta_s)),
  beta_dau_est = sprintf("%1.3f", mean(beta_d)),
  beta_dau_se  = sprintf("%1.3f", sd(beta_d)),
  beta_diff_est = sprintf("%1.3f", mean(beta_diff)),
  beta_diff_se = sprintf("%1.3f", sd(beta_diff))
)

write_json(report, path = paste0("./temp/", measure, ".json"), pretty = TRUE)
