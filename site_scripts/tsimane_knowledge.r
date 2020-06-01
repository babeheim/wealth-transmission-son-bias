
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "tsimane_knowledge"

d0 <- read.csv("./data/tsimane_knowledge.csv", stringsAsFactors = FALSE)

#all folks
d <- data.frame(
  c(d0$midpid, d0$dads_midpid, d0$moms_midpid),
  c(d0$percpossskills, d0$dads_percpossskills, d0$moms_percpossskills),
  c(d0$age, d0$agef, d0$agem)/10,
  c(d0$male, rep(1, nrow(d0)), rep(0, nrow(d0))),
  c(d0$comunidad, d0$dads_comunidad, d0$moms_comunidad))
names(d) <- c("id", "skill", "age", "male", "com")
drop <- which(duplicated(d$id) | is.na(d$id) | is.na(d$skill * d$age))
d <- d[-drop, ]
d$com <- as.numeric(factor(d$com))

# folks and their parents
drop <- which(
  duplicated(d0$midpid) |
  is.na(d0$midpid) |
  is.na(d0$age * d0$agef * d0$agem * d0$percpossskills)
)
d1 <- d0[-drop, ]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
d1$agem <- d1$agem/10
names(d1) <- c(
  "id", "age", "com", "skill",
  "idf", "agef", "comf", "skillf",
  "idm", "agem", "comm", "skillm",
  "male")
d1$com <- as.numeric(factor(d1$com))
d1$comf <- as.numeric(factor(d1$comf))
d1$comm <- as.numeric(factor(d1$comm))

(mean(d1$skill[which(d1$male == 1)])-mean(d1$skill[which(d1$male == 0)]))/sd(d1$skill)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  n = nrow(d),
  n1 = nrow(d1),
  skill = d$skill,
  age = d$age,
  male = d$male,
  com = d$com,
  skill_f = d1$skillf,
  age_f = d1$agef,
  com_f = d1$comf,
  skill_m = d1$skillm,
  age_m = d1$agem,
  com_m = d1$comm,
  skill_o = d1$skill,
  age_o = d1$age,
  male_o = d1$male,
  com_o = d1$com
)

expect_true(data_list$n == 267)
expect_true(data_list$n1 == 117)



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

for(i in 1:n_sim) {

  #find parent correlation
  dparents <- data.frame(post$Dev_f[i, ], post$Dev_m[i, ])
  dparents <- dparents[-which(duplicated(data.frame(d1$idf, d1$idm))), ]
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
beta_s_est <- 0.027
beta_s_se  <- 0.148
beta_d_est <- (-0.024)
beta_d_se  <- 0.135

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- 0.051
beta_diff_se <- 0.2

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
