
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "datoga_rs"

d0 <- read.csv("./data/datoga_rs.csv", stringsAsFactors = FALSE)

d0$male <- d0$sex
d0$male[which(d0$male == "s")] <-1
d0$male[which(d0$male == "d")] <-0
d0$male <- as.numeric(d0$male)

#Sample looks good. Doing fathers and offspring separately, because removing 0 RS individuals takes out too much from the sample size.
#df: fathers only
df <- data.frame(d0$fid, d0$agef, d0$rsf)
names(df) <- c("id", "age", "rs")
drop <- which(duplicated(df$id))
df <- df[-drop, ]
df$age <- df$age/10
D_f <- cbind(rep(1, nrow(df)), df$age)

#do: the offspring
drop <- which(is.na(d0$rs) | duplicated(d0$id)) #going to include the 17.5'ers, for now. There are 19 of them! Not much variance, though...
do <- d0[-drop, ]
do$age <- do$age/10
D_o <- cbind(do$male, 1-do$male, do$male * do$age, (1-do$male) * do$age)

#d1: folks with parents
drop <- which(is.na(d0$rs * d0$rsf) | duplicated(d0$id)) 
d1 <- d0[-drop, ]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
D_f1 <- cbind(rep(1, nrow(d1)), d1$agef)
D_o1 <- cbind(d1$male, 1-d1$male, d1$male * d1$age, (1-d1$male) * d1$age)
(mean(d1$rs[which(d1$male == 1)])-mean(d1$rs[which(d1$male == 0)]))/sd(d1$rs)
mean(d1$age * 10)


#make list for stan
data_list <- list(
  nf = nrow(df), 
  no = nrow(do), 
  n1 = nrow(d1), 
  rs_f = df$rs, 
  rs_o = do$rs, 
  rs_f1 = d1$rsf, 
  rs_o1 = d1$rs, 
  D_f = D_f, 
  D_o = D_o, 
  D_f1 = D_f1, 
  D_o1 = D_o1, 
  male_o = do$male, 
  male_o1 = d1$male
)

expect_true(data_list$nf == 54)
expect_true(data_list$no == 135)
expect_true(data_list$n1 == 135)



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
  #do the regression
  par_temp <- post$Dev_f[i, ]
  dtemp <- data.frame(post$Dev_o[i, ], par_temp)
  names(dtemp) <- c("Dev_o", "Dev_f")
  
  #sons
  dtemps <- dtemp[which(d1$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * Dev_f, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(d1$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * Dev_f, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]
}

# reported values here
beta_s_est <- 0.176
beta_s_se  <- 0.126
beta_d_est <- 0.197
beta_d_se  <- 0.167

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.021)
beta_diff_se <- 0.208

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
