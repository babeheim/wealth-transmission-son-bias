
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "gambia_height"

d0 <- read.csv("./data/gambia_height.csv", stringsAsFactors = FALSE)

names(d0)[which(names(d0) == "sex")] <- "male"
d0$male[which(d0$male == 2)] <- 0

#d: anyone who has weight and age data. the purpose of this is to get a good age model.
id <- c(d0$id, d0$fid, d0$mid)
height <- c(
  sqrt(d0$weight/d0$bmi),
  sqrt(d0$weightf/d0$bmif),
  sqrt(d0$weightm/d0$bmim)) * 100
age <- c(d0$age, d0$agef, d0$agem)/10
male <- c(d0$male, rep(1, nrow(d0)), rep(0, nrow(d0)))
d <- data.frame(id, height, age, male)
drop <- which(
  is.na(d$age) |
  is.na(d$height) |
  d$age < 1.8 |
  duplicated(d$id)
)
d <- d[-drop, ]
#use logs
d$height <- log(d$height)
D_mu <- cbind(1 - d$male, d$male, d$age, d$age^2)
D_sig <- cbind(1 - d$male, d$male)

#now folks with parents
height <- log(sqrt(d0$weight/d0$bmi) * 100)
heightf <- log(sqrt(d0$weightf/d0$bmif) * 100)
heightm <- log(sqrt(d0$weightm/d0$bmim) * 100)
do <- data.frame(d0, height, heightf, heightm)
drop <- which(
  is.na(do$age * do$height * do$male * do$agef *
    do$heightm * do$heightf * do$agem/1000) |
  do$age < 18 |
  duplicated(do$id)
)
do <- do[-drop, ]
do$age <- do$age/10
do$agef <- do$agef/10
do$agem <- do$agem/10
D_mu_o <- cbind(1 - do$male, do$male, do$age, do$age^2)
D_sig_o <- cbind(1 - do$male, do$male)
D_mu_f <- cbind(rep(0, nrow(do)), rep(1, nrow(do)), do$agef, do$agef^2)
D_sig_f <- cbind(rep(0, nrow(do)), rep(1, nrow(do)))
D_mu_m <- cbind(rep(1, nrow(do)), rep(0, nrow(do)), do$agem, do$agem^2)
D_sig_m <- cbind(rep(1, nrow(do)), rep(0, nrow(do)))


#make list for stan
data_list <- list(
  n = nrow(d), 
  no = nrow(do), 
  height = d$height, 
  height_o = do$height, 
  height_f = do$heightf, 
  height_m = do$heightm, 
  D_mu = D_mu, 
  D_sig = D_sig, 
  D_mu_o = D_mu_o, 
  D_sig_o = D_sig_o, 
  D_mu_f = D_mu_f, 
  D_sig_f = D_sig_f, 
  D_mu_m = D_mu_m, 
  D_sig_m = D_sig_m
)

expect_true(data_list$n == 2355)
expect_true(data_list$no == 817)



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

for(i in 1:n_sim)
{
  #find parent correlation
  dparents <- data.frame(post$Dev_f[i, ], post$Dev_m[i, ])
  dparents <- dparents[-which(duplicated(data.frame(do$fid, do$mid))), ]
  names(dparents) <- c("Dev_f", "Dev_m")
  r <- sample.naive.posterior(mle2(Dev_f ~ dnorm(r * Dev_m, sqrt(1-r^2)), start = list(r=0), data = dparents), 1)[[1]]
  
  #now do the regression
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])/sqrt(1/2 * (1+r))
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
beta_s_est <- 0.435
beta_s_se  <- 0.039
beta_d_est <- 0.559
beta_d_se  <- 0.033

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.124)
beta_diff_se <- (0.049)

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
