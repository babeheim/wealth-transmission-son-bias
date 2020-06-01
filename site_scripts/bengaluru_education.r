
rm(list = setdiff(ls(), "analyses_to_run"))

source("./project_support.r")

measure <- "bengaluru_education"

d0 <- read.csv("./data/bengaluru_education.csv", stringsAsFactors = FALSE)

#just the kids
age <- c(2002-d0$child.year.birth)
soned <- d0$son.education.years
soned[which(is.na(soned))] <- 0
daughtered <- d0$daughter.education.years
daughtered[which(is.na(daughtered))] <- 0
ed <- soned + daughtered
ed[which(ed == 0)] <- NA
male <- rep(NA, nrow(d0))
male[which(!is.na(d0$son.RS))] <- 1
male[which(!is.na(d0$daughter.RS))] <- 0
do <- data.frame(male, age, ed)
drop <- which(is.na(do$male * do$age * do$ed) | do$age<18)
do <- do[-drop, ]
do$age <- do$age/10
D_mu_o <- cbind(1-do$male, (1-do$male) * do$age, (1-do$male) * do$age^2, do$male, do$male * do$age)
#go ahead and ignore the exponential; it works fine.

#now just the parents
age <- c(2002-d0$father.year.birth, 2002-d0$mother.year.birth)
id <- c(d0$RESPONDENT.ID+1000, d0$RESPONDENT.ID)
ed <- c(d0$father.education.years, d0$mother.education.years)
male <- c(rep(1, nrow(d0)), rep(0, nrow(d0)))
dp <- data.frame(id, male, age, ed)
drop <- which(is.na(dp$male * dp$age * dp$ed) | duplicated(dp$id))
dp <- dp[-drop, ]
dp$age <- dp$age/10
D_mu_p <- cbind(1-dp$male, dp$male)

#now folks and their parents
age <- c(2002-d0$child.year.birth)
soned <- d0$son.education.years
soned[which(is.na(soned))] <- 0
daughtered <- d0$daughter.education.years
daughtered[which(is.na(daughtered))] <- 0
ed <- soned + daughtered
ed[which(ed == 0)] <- NA
male <- rep(NA, nrow(d0))
male[which(!is.na(d0$son.RS))] <- 1
male[which(!is.na(d0$daughter.RS))] <- 0
edf <- d0$father.education.years
edm <- d0$mother.education.years
agef <- 2002-d0$father.year.birth
agem <- 2002-d0$mother.year.birth
idf <- d0$RESPONDENT.ID+1000
idm <- d0$RESPONDENT.ID
d1 <- data.frame(male, age, ed, agef, edf, agem, edm, idf, idm)
drop <- which(is.na(d1$male * d1$age * d1$ed * d1$edf * d1$edm * d1$agef * d1$agem) | d1$age<18)
d1 <- d1[-drop, ]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
d1$agem <- d1$agem/10
D_mu_o1 <- cbind(1-d1$male, (1-d1$male) * d1$age, (1-d1$male) * d1$age^2, d1$male, d1$male * d1$age)
D_mu_f1 <- cbind(rep(0, nrow(d1)), rep(1, nrow(d1)))
D_mu_m1 <- cbind(rep(1, nrow(d1)), rep(0, nrow(d1)))

#make list for stan
data_list <- list(
  no = nrow(do), 
  np = nrow(dp), 
  n1 = nrow(d1), 
  ed_o = do$ed, 
  ed_p = dp$ed, 
  ed_o1 = d1$ed, 
  ed_f1 = d1$edf, 
  ed_m1 = d1$edm, 
  D_mu_o = D_mu_o, 
  D_mu_p = D_mu_p, 
  D_mu_o1 = D_mu_o1, 
  D_mu_f1 = D_mu_f1, 
  D_mu_m1 = D_mu_m1
  )

expect_true(data_list$no == 588)
expect_true(data_list$np == 767)
expect_true(data_list$n1 == 562)



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
beta_s_est <- 0.668
beta_s_se  <- 0.034
beta_d_est <- 0.71
beta_d_se  <- 0.03

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.042)
beta_diff_se <- 0.044

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
