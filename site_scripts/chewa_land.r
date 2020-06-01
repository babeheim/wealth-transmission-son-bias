
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "chewa_land"

d0 <- read.csv("./data/chewa_land.csv", stringsAsFactors = FALSE)

d0$male[which(d0$male == 2)] <- 0

#all folks
#don't need to grab parents manually; already in offspring column
#remove folks who aren't involved in ag, and who aren't heads of household
drop <- which(is.na(d0$hhgardsize) | d0$mainsub != 1 | d0$relathh > 2 | d0$age<18)
d <- d0[-drop, ]
d$hhgardsize <- log(d$hhgardsize)
d$age <- d$age/10
D_mu <- cbind(rep(1, nrow(d)), d$age, d$age^2)
D_sig <- cbind(rep(1, nrow(d)), d$age)

#Folks with mothers (not enough fathers with age data, although this could be 
# imputed accurately. doesn't matter though, since they share garden with mother)
drop <- which(is.na(d0$hhgardsize * d0$agem * d0$gardsizem) | d0$mainsub != 1 | d0$relathh > 2 | d0$age<18 | d0$agem <0)
d1 <- d0[-drop, ]
d1$hhgardsize <- log(d1$hhgardsize)
d1$gardsizem <- log(d1$gardsizem)
d1$age <- d1$age/10
d1$agem <- d1$agem/10
D_mu_o <- cbind(rep(1, nrow(d1)), d1$age, d1$age^2)
D_sig_o <- cbind(rep(1, nrow(d1)), d1$age)
D_mu_m <- cbind(rep(1, nrow(d1)), d1$agem, d1$agem^2)
D_sig_m <- cbind(rep(1, nrow(d1)), d1$agem)

(mean(d1$hhgardsize[which(d1$male == 1)])-mean(d1$hhgardsize[which(d1$male == 0)]))/sd(d1$hhgardsize)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  n = nrow(d), 
  n1 = nrow(d1), 
  hhgardsize = d$hhgardsize, 
  hhgardsize_o = d1$hhgardsize, 
  hhgardsize_m = d1$gardsizem, 
  D_mu = D_mu, 
  D_sig = D_sig, 
  D_mu_o = D_mu_o, 
  D_sig_o = D_sig_o, 
  D_mu_m = D_mu_m, 
  D_sig_m = D_sig_m
)

expect_true(data_list$n == 1906)
expect_true(data_list$n1 == 151)



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
beta_s_est <- (-0.137)
beta_s_se  <- 0.251
beta_d_est <- 0.058
beta_d_se  <- 0.097

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.195)
beta_diff_se <- 0.268

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
