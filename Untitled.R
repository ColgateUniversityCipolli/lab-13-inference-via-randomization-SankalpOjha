library(tidyverse)
library(xtable)
library(VGAM)
library(e1071)
library(pwr)
library(effectsize)

data <- read_csv("zebrafinches.csv")


################################################################################
## Question 1A
################################################################################
library(tidyverse)
library(xtable)
library(VGAM)
library(e1071)
library(pwr)
library(effectsize)

data <- read_csv("zebrafinches.csv")

mu0 <- 0
t.stat <- t.test(data$further,
                    mu = mu0,
                    alternative = "less")
t.further <- t.stat$statistic

n <- length(data$further)

skewness <- skewness(data$further)
fz <- dnorm(t.further)
Fz <- pnorm(t.further)

(edgeworth.error <- (skewness/sqrt(n)) * (((2*(t.further)^2 + 1)/6)*(fz)))

(probability <- Fz + edgeworth.error)
################################################################################
## Question 1B
################################################################################
t.vals <- seq(-10,10, length = 1000)
fzb <- dnorm(t.vals)

error.vals <- (skewness/sqrt(n)) * (((2*(t.vals)^2 + 1)/6)*(fzb))

tvals.error <- tibble(
  t = t.vals,
  error = error.vals
)

ggplot()+
  geom_line(data = tvals.error,
            aes(x = t, y = error),
            color = "red")+
  theme_bw()+
  labs(title = "Edgeworth Approximation For Error",
       x = "T Values (-10,10)",
       y = "Error")
################################################################################
## Question 1C
################################################################################
val.t <- qnorm(0.05)
fzc <- dnorm(val.t)

(n <- ((skewness/(6*(0.10*0.05))) * (2*val.t^2 + 1) *fzc)^2)
################################################################################
## Question 2A
################################################################################
library(boot)

R <- 10000
sd.closer <- sd(data$closer)
num.closer <- length(data$closer)
resamples.closer <-  tibble(tstat=rep(NA, R),
                            xbar=rep(NA, R))
for(i in 1:R){
  curr.resample <- sample(x = data$closer,
                          size = num.closer,
                          replace = T)
  resamples.closer$tstat[i] <- (mean(curr.resample)-0)/(sd.closer/sqrt(num.closer))
  resamples.closer$xbar[i] <- mean(curr.resample)
}


closer.delta <- mean(resamples.closer$tstat) - 0

resamples.null.closer <- resamples.closer |>
  mutate(shifted = tstat - closer.delta)

R <- 10000
sd.further <- sd(data$further)
num.further <- length(data$further)
resamples.further <- tibble(tstat=rep(NA, R),
                            xbar=rep(NA, R))
for(i in 1:R){
  curr.resample <- sample(x = data$further,
                          size = num.further,
                          replace = T)
  resamples.further$tstat[i] <- (mean(curr.resample)-0)/(sd.further/sqrt(num.further))
  resamples.further$xbar[i] <- mean(curr.resample)
}


further.delta <- mean(resamples.further$tstat) - 0

resamples.null.further <- resamples.further |>
  mutate(shifted = tstat - further.delta)

R <- 10000
sd.diff <- sd(data$diff)
num.diff <- length(data$diff)
resamples.diff <-  tibble(tstat=rep(NA, R),
                          xbar=rep(NA, R))
for(i in 1:R){
  curr.resample <- sample(x = data$diff,
                          size = num.diff,
                          replace = T)
  resamples.diff$tstat[i] <- (mean(curr.resample)-0)/(sd.diff/sqrt(num.diff))
  resamples.diff$xbar[i] <- mean(curr.resample)
}

diff.delta <- mean(resamples.diff$tstat) - 0
resamples.null.diff <- resamples.diff |>
  mutate(shifted = tstat - diff.delta)
################################################################################
## Question 2B
################################################################################
# T-Test P-Value
t.p.close <- t.test(data$closer,
                     mu = 0,
                     alternative = "greater")$p.value

t.p.further <- t.test(data$further,
                      mu = 0,
                      alternative = "less")$p.value

t.p.diff <- t.test(data$diff,
                   mu = 0,
                   alternative = "two.sided")$p.value

# Bootstrap P-Value
boot.p.closer <- mean(resamples.null.closer$shifted >= closer.delta)

boot.p.further <- mean(resamples.null.further$shifted <= further.delta)

low <- 0 - diff.delta
high <- 0 + diff.delta
p.low <- mean(resamples.null.diff$shifted <= low)
p.high <- mean(resamples.null.diff$shifted >= high)
boot.p.diff <- p.low + p.high

(t.p.close)
(boot.p.closer)
(t.p.further)
(boot.p.further)
(t.p.diff)
(boot.p.diff)
################################################################################
## Question 2C
################################################################################
boot.closer.5th <- quantile(resamples.null.closer$shifted, 0.05)
t.closer.5th <- qt(0.05, df = num.closer -1)

boot.further.5th <- quantile(resamples.null.further$shifted, 0.05)
t.further.5th <- qt(0.05, df = num.further -1)

boot.diff.5th <- quantile(resamples.null.diff$shifted, 0.05)
t.diff.5th <- qt(0.05, df = num.diff -1)

(boot.closer.5th)
(t.closer.5th)
(boot.further.5th)
(t.further.5th)
(boot.diff.5th)
(t.diff.5th)
################################################################################
## Question 2D
################################################################################
boot.CI.closer <- quantile(resamples.null.closer$xbar, c(0.025, 0.975))
t.CI.closer <- t.test(x=data$closer, 
                      mu = 0, 
                      conf.level = 0.95, 
                      alternative = "two.sided")$conf.int

boot.CI.further <- quantile(resamples.null.further$xbar, c(0.025, 0.975))
t.CI.further <- t.test(x=data$further, 
                       mu = 0, 
                       conf.level = 0.95, 
                       alternative = "two.sided")$conf.int

boot.CI.diff <- quantile(resamples.null.diff$xbar, c(0.025, 0.975))
t.CI.diff <- t.test(x=data$diff, 
                    mu = 0, 
                    conf.level = 0.95, 
                    alternative = "two.sided")$conf.int

(boot.CI.closer)
(t.CI.closer)
(boot.CI.further)
(t.CI.further)
(boot.CI.diff)
(t.CI.diff)
################################################################################
## Question 3A
################################################################################
randomization.dist <- function(data, mu0, R = 10000) {
  random <- tibble(means = rep(NA, R))
  
  shifted.data <- data - mu0
  
  for(i in 1:R) {
    curr.random <- shifted.data * sample(x = c(-1, 1), 
                                         size = length(shifted.data), 
                                         replace = TRUE)
    random$means[i] <- mean(curr.random)
  }
  
  random <- random |>
    mutate(means = means + mu0)
  
  return(random)
}

closer.rand <- randomization.dist(data$closer, 0)
farther.rand <- randomization.dist(data$further, 0)
diff.rand <- randomization.dist(data$diff, 0)

(mean(closer.rand$means))
(mean(farther.rand$means))
(mean(diff.rand$means))
################################################################################
## Question 3B
################################################################################
closer.mean <- mean(data$closer)
closer.p.random <- mean(closer.rand$means >= closer.mean)

further.mean <- mean(data$further)
further.p.random <- mean(farther.rand$means <= further.mean)

diff.delta <- abs(mean(data$diff) - mu0)
mirror <- mu0 - diff.delta
xbar <- mu0 + diff.delta

diff.p.random <- mean(diff.rand$means <= mirror) + 
                 mean(diff.rand$means >= xbar)

(closer.p.random)
(further.p.random)
(diff.p.random)
################################################################################
## Question 3C
################################################################################
CI.finder <- function(data, R = 10000){
  
  mu0.iter <- 0.01
  mu0.lower <- mean(data)
  mu0.higher <- mean(data)
  
  repeat{
    random <- tibble(means = rep(NA, R))
    
    shifted.data <- data - mu0.lower
    
    for(i in 1:R) {
      curr.random <- shifted.data * sample(x = c(-1, 1), 
                                           size = length(shifted.data), 
                                           replace = TRUE)
      random$means[i] <- mean(curr.random)
    }
    
    random <- random |>
      mutate(means = means + mu0.lower)
    
    delta <- abs(mean(data) - mu0.lower)
    low <- mu0.lower - delta 
    high <- mu0.lower + delta   
    p.val <- mean(random$means <= low) +
             mean(random$means >= high)
    
    if(p.val < 0.05){
      break
    }else{
      mu0.lower <- mu0.lower - mu0.iter
    }
  }
  
  repeat{
    random <- tibble(means = rep(NA, R))
    
    shifted.data <- data - mu0.higher
    
    for(i in 1:R) {
      curr.random <- shifted.data * sample(x = c(-1, 1), 
                                           size = length(shifted.data), 
                                           replace = TRUE)
      random$means[i] <- mean(curr.random)
    }
    
    random <- random |>
      mutate(means = means + mu0.higher)
    
    delta <- abs(mean(data) - mu0.higher)
    low <- mu0.higher - delta 
    high <- mu0.higher + delta   
    p.val <- mean(random$means <= low) +
      mean(random$means >= high)
    
    if(p.val < 0.05){
      break
    }else{
      mu0.higher <- mu0.higher + mu0.iter
    }
  }
  
  return(c(mu0.lower, mu0.higher))
}

(close.CI.randomization <- CI.finder(data$closer))
(further.CI.randomization <- CI.finder(data$further))
(diff.CI.randomization <- CI.finder(data$diff))

