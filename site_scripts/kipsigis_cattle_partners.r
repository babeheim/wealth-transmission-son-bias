
rm(list = ls())

source("./project_support.r")

measure <- "kipsigis_cattle_partners"

d0 <- read.csv("./data/kipsigis_cattle_partners.csv", stringsAsFactors = FALSE)

#just fathers
df <- data.frame(d0$fid, d0$agef, d0$partnersf)
names(df) <- c("id", "age", "partners")
drop <- which(duplicated(df$id) | is.na(df$id) | is.na(df$age * df$partners))
df <- df[-drop, ]
df$age <- df$age/10

#the kids
drop <- which(duplicated(d0$id) | is.na(d0$age * d0$partners))
drop <- unique(
  c(drop, which(
    d0$male == 1 &
    d0$land %in% d0$landf &
    d0$livestock %in% d0$livestockf &
    d0$partners %in% d0$partnersf &
    d0$RS %in% d0$RSf &
    d0$has_parent_by)
  )
)
expect_true(length(drop) == 345)

# the latter removed known fathers
do <- d0[-drop, ]
do$age <- do$age/10
D_o <- cbind(rep(1, nrow(do)), do$age)

#those with parents
drop <- which(duplicated(d0$id) | is.na(d0$age/1000 * d0$partners * d0$partnersf))
d1 <- d0[-drop, ]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
D_o1 <- cbind(rep(1, nrow(d1)), d1$age)

(mean(d1$partners[which(d1$male == 1)])-mean(d1$partners[which(d1$male == 0)]))/sd(d1$partners)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  nf = nrow(df), 
  no = nrow(do), 
  n1 = nrow(d1), 
  partners_f = df$partners, 
  partners_o = do$partners, 
  partners_f1 = d1$partnersf, 
  partners_o1 = d1$partners, 
  D_o = D_o, 
  D_o1 = D_o1
)

expect_true(data_list$nf == 25)
expect_true(data_list$no == 130)
expect_true(data_list$n1 == 102)



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
beta_s_est <- (-0.02)
beta_s_se  <- 0.148
beta_d_est <- 0.161
beta_d_se  <- 0.206

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.181)
beta_diff_se <- 0.249

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
