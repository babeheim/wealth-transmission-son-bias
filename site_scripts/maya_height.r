
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "maya_height"

d0 <- read.csv("./data/maya_height.csv", stringsAsFactors = FALSE)

names(d0)[which(names(d0) == "sex")] <- "male"
d0$male[which(d0$male == 2)] <- 0

# gonna combine all folks under assumption that selection on height is very weak.

# there's a big outlier in this dataset - about 4.5 SD's out. gonna remove it and see what the effect is.
d0 <- d0[-88, ]
# RESULT: no big effect, but probably better removed

# using log height! all folks:
d <- data.frame(d0$male, as.numeric(paste(d0$age))/10, log(d0$height))  #notice log
names(d) <- c("male", "age", "height")
drop <- which(is.na(d$age * d$male * d$height))
d <- d[-drop, ]
D_mu <- cbind(1 - d$male, d$male, d$age)

#now folks with both parents (could save a lot by imputing...)
#construct parent heights to get sample sizes
heightf <- rep(NA, nrow(d0))
heightm <- rep(NA, nrow(d0))
agef <- rep(NA, nrow(d0))
agem <- rep(NA, nrow(d0))
for (i in 1:nrow(d0)) {
  if (length(which(d0$id == d0$fid[i])) > 0) {
    heightf[i] <- d0$height[which(d0$id == d0$fid[i])]
    agef[i] <- as.numeric(paste(d0$age[which(d0$id == d0$fid[i])]))
  }
  if (length(which(d0$id == d0$mid[i])) > 0) {
    heightm[i] <- d0$height[which(d0$id == d0$mid[i])]
    agem[i] <- as.numeric(paste(d0$age[which(d0$id == d0$mid[i])]))
  }
}
do <- data.frame(d0, log(heightf), log(heightm), agef/10, agem/10)
do$age <- as.numeric(paste(do$age))/10
do$height <- log(do$height)
names(do) <- c(names(d0), "heightf", "heightm", "agef", "agem")
drop <- which(is.na(do$male * do$height * do$heightf * do$heightm * do$age * do$agef * do$agem))
do <- do[-drop, ]
D_mu_o <- cbind(1-do$male, do$male, do$age)
D_mu_f <- cbind(rep(0, nrow(do)), rep(1, nrow(do)), do$agef)
D_mu_m <- cbind(rep(1, nrow(do)), rep(0, nrow(do)), do$agem)


#make list for stan
data_list <- list(
  n = nrow(d), 
  no = nrow(do), 
  height = d$height, 
  height_o = do$height, 
  height_f = do$heightf, 
  height_m = do$heightm, 
  D_mu = D_mu, 
  D_mu_o = D_mu_o, 
  D_mu_f = D_mu_f, 
  D_mu_m = D_mu_m
)

expect_true(data_list$n == 254)
expect_true(data_list$no == 142)



tic(paste("fit", measure, "stan model"))

if (!file.exists(file.path("temp", paste0(measure, ".stanfit")))){

  m1 <- stan(file = file.path("site_models", paste0(measure, ".stan")), data = data_list, 
    iter = n_iter, warmup = n_warmup, chains = 1, init = init)

  save(m1, file = file.path("temp", paste0(measure, ".stanfit")))

}
toc(log = TRUE)



load(file.path("temp", paste0(measure, ".stanfit")))

post <- extract(m1, permute = TRUE)

n_sim <- min(10000, length(post$lp__))

beta_s <- rep(NA, n_sim)
beta_d <- rep(NA, n_sim)

for (i in 1:n_sim) {

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
beta_s_est <- 0.4
beta_s_se  <- 0.097
beta_d_est <- 0.629
beta_d_se  <- 0.071

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.23)
beta_diff_se <- 0.109

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
