
rm(list = ls())

source("./project_support.r")

measure <- "datoga_livestock"

d0 <- read.csv("./data/datoga_livestock.csv", stringsAsFactors = FALSE)

d0$male <- d0$sex
d0$male[which(d0$male == "s")] <-1
d0$male[which(d0$male == "d")] <-0
d0$male <- as.numeric(d0$male)

#clear difference in cattle holdings between fathers and offspring.
#gonna work in logs

#df: fathers only
df <- data.frame(d0$fid, d0$agef, d0$livestockf)
names(df) <- c("id", "age", "livestock")
drop <- which(duplicated(df$id))
df <- df[-drop, ]
df$age <- df$age/10
df$livestock <- log(df$livestock)
D_f <- cbind(rep(1, nrow(df)), df$age, df$age^2)

#do: the offspring
drop <- which(is.na(d0$livestock) | duplicated(d0$id))
# going to include the 17.5'ers, for now.
# There are 19 of them! Not much variance, though...
do <- d0[-drop, ]
do$age <- do$age/10
do$livestock <- log(do$livestock)
D_o <- cbind(rep(1, nrow(do)), do$age, do$age^2)

#d1: folks with parents
drop <- which(is.na(d0$livestock * d0$livestockf) | duplicated(d0$id)) 
d1 <- d0[-drop, ]
d1$age <- d1$age/10
d1$agef <- d1$agef/10
d1$livestock <- log(d1$livestock)
d1$livestockf <- log(d1$livestockf)
D_f1 <- cbind(rep(1, nrow(d1)), d1$agef, d1$agef^2)
D_o1 <- cbind(rep(1, nrow(d1)), d1$age, d1$age^2)

(mean(d1$livestock[which(d1$male == 1)])-mean(d1$livestock[which(d1$male == 0)]))/sd(d1$livestock)
mean(d1$age * 10)

#make list for stan
data_list <- list(
  nf = nrow(df), 
  no = nrow(do), 
  n1 = nrow(d1), 
  livestock_f = df$livestock, 
  livestock_o = do$livestock, 
  livestock_f1 = d1$livestockf, 
  livestock_o1 = d1$livestock, 
  D_f = D_f, 
  D_o = D_o, 
  D_f1 = D_f1, 
  D_o1 = D_o1
)

# add checks for the reported counts
expect_true(sum(d1$male == 1) == 95)
expect_true(sum(d1$male == 0) == 40)
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
beta_s_est <- 0.537
beta_s_se  <- 0.08
beta_d_est <- 0.544
beta_d_se  <- 0.118

expect_true(abs(mean(beta_s) - beta_s_est) < beta_tol)
expect_true(abs(sd(beta_s) - beta_s_se) < beta_tol)
expect_true(abs(mean(beta_d) - beta_d_est) < beta_tol)
expect_true(abs(sd(beta_d) - beta_d_se) < beta_tol)

beta_diff_est <- (-0.007)
beta_diff_se <- 0.132

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
