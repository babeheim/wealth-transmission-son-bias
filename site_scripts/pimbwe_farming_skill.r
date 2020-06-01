
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "pimbwe_farming_skill"

d0 <- read.csv("./data/pimbwe_farming_skill.csv", stringsAsFactors = FALSE)

# there are a few cases here where the parent and offspring ages
# are unusually close (5<x<10). But age doesn't appear to be super
# important here, so I'm not going to worry about it.
# imputation could be a good idea here. parent skills are strongly correlated
# doesn't appear to be much difference between parent and
# offspring skill-to-age relationship,
# so I'm going to fit one model for everyone.

d0$sex[which(d0$sex == "m")] <- 1
d0$sex[which(d0$sex == "f")] <- 0
d0$sex <- as.numeric(d0$sex)
d0$male <- d0$sex


# All folks, without repeats
d <- data.frame(
  c(d0$fid, d0$mid, d0$id),
  c(d0$skillf, d0$skillm, d0$skill),
  c(d0$agef, d0$agem, d0$age),
  c(rep(1, nrow(d0)), rep(0, nrow(d0)), d0$male)
)
names(d) <- c("id", "skill", "age", "male")
drop <- which(is.na(d$age * d$skill))
d <- d[-drop, ]
d <- d[-which(duplicated(d$id)), ]
d$age <- d$age/10
D_all_sig <- cbind(rep(1, nrow(d)), d$age)

#Now folks with their parents
drop <- which(is.na(d0$age * d0$skill * d0$agef * d0$skillf * d0$agem * d0$skillm))
d1 <- d0[-drop, ]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
d1$agem <- d1$agem/10
D_f_sig <- cbind(rep(1, nrow(d1)), d1$agef)
D_m_sig <- cbind(rep(1, nrow(d1)), d1$agem)
D_o_sig <- cbind(rep(1, nrow(d1)), d1$age)

(mean(d1$skill[which(d1$male == 1)])-mean(d1$skill[which(d1$male == 0)]))/sd(d1$skill)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  n = nrow(d), 
  n1 = nrow(d1), 
  skill = d$skill, 
  D_all_sig = D_all_sig, 
  skill_f = d1$skillf, 
  skill_m = d1$skillm, 
  skill_o = d1$skill, 
  D_f_sig = D_f_sig, 
  D_m_sig = D_m_sig, 
  D_o_sig = D_o_sig, 
  male = d1$male
)

expect_true(data_list$n == 508)
expect_true(data_list$n1 == 138)



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
  dparents <- dparents[-which(duplicated(data.frame(d1$fid, d1$mid))), ]
  names(dparents) <- c("Dev_f", "Dev_m")
  r <- sample.naive.posterior(mle2(Dev_f ~ dnorm(r * Dev_m, sqrt(1-r^2)), start = list(r=0), data = dparents), 1)[[1]]
  
  #now do the regression
  midpar_temp <- 0.5 * (post$Dev_f[i, ]+post$Dev_m[i, ])/sqrt(1/2 * (1+r))
  dtemp <- data.frame(post$Dev_o[i, ], midpar_temp)
  names(dtemp) <- c("Dev_o", "midpar_temp")
  
  #sons
  dtemps <- dtemp[which(d1$male == 1), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtemps)
  beta_s[i] <- sample.naive.posterior(m_temp, 1)[[1]]
  
  #daughters
  dtempd <- dtemp[which(d1$male == 0), ]
  m_temp <- mle2(Dev_o ~ dnorm(beta * midpar_temp, sqrt(1-beta^2)), start = list(beta=0), data = dtempd)
  beta_d[i] <- sample.naive.posterior(m_temp, 1)[[1]]
}

# reported values here
beta_s_est <- 0.078
beta_s_se  <- 0.138
beta_d_est <- 0.206
beta_d_se  <- 0.118

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.127)
beta_diff_se <- 0.18

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
