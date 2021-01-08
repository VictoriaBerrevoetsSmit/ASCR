# Packages for data processing and nice functions
library(copula)
library(tidyverse)
library(dplyr)
library(knitr)
library(parallel)
numCores <- detectCores()
# Packages for nice plots
library(ggplot2)
library(gridExtra)
library(ggthemes)


data.original <- read.csv('insurance.csv')
data.original$ID <- as.factor(data.original$ID)
set.seed(18101996)

ggplot(data.original, aes(WC, PLI)) + geom_bin2d()  + theme_void() + scale_color_colorblind() + theme(legend.position = "none")

ggplot(data.original) + geom_point(aes(WC, PLI))  + annotate("text", x=75, y=52, label="Large Claim") +theme_classic()

ppli <- ggplot(data.original) + geom_histogram(aes(PLI), bins=50) + theme_classic()
pwc <- ggplot(data.original) + geom_histogram(aes(WC), bins=50) + theme_classic()
grid.arrange(ppli, pwc, nrow=1)

sum.stats.pli <- c(minimum = min(data.original[,c("PLI")]), quantile(data.original[,c("PLI")], 0.25), median = median(data.original[,c("PLI")]), quantile(data.original[,c("PLI")], 0.75), maximum = max(data.original[,c("PLI")]), mean = mean(data.original[,c("PLI")]), sd = sd(data.original[,c("PLI")]))
sum.stats.wc <- c(minimum = min(data.original[,c("WC")]), quantile(data.original[,c("WC")], 0.25), median = median(data.original[,c("WC")]), quantile(data.original[,c("WC")], 0.75), maximum = max(data.original[,c("WC")]), mean = mean(data.original[,c("WC")]), sd = sd(data.original[,c("WC")]))
sum.stats <- data.frame(PLI = sum.stats.pli, WC = sum.stats.wc)
kable(t(data.frame(sum.stats)), caption = 'Summary Statistics of Original Data')

correlation.test <- cor.test(data.original[,c("WC")], data.original[,c("PLI")])
corrtest <- data.frame(estimate = correlation.test$estimate, statistic = correlation.test$statistic, df = correlation.test$parameter, p.value = correlation.test$p.value)

kable((corrtest), caption = 'Correlation test of Original Data')


estimate.mu <- function(mu, x){
  (sum(log(x))/length(x)) - mu}

estimate.sigma <- function(sigma, x){
  sum((log(x)-(sum(log(x))/length(x)))^2)/length(x) - sigma^2}


estimate.params <- function(data){
  #In: a dataframe
  #Out: a list with parameters
  
  X1 <- data[,c("PLI")]
  X2 <- data[,c("WC")]
  
  mu1 <- uniroot(estimate.mu, c(min(log(X1)), max(log(X1))), x=X1)$root
  sigma1 <- uniroot(estimate.sigma, c(0.1, 4*sd(log(X1))), x=X1)$root
  mu2 <- uniroot(estimate.mu, c(min(log(X2)), max(log(X2))), x=X2)$root
  sigma2 <- uniroot(estimate.sigma, c(0.1, 4*sd(log(X2))), x=X2)$root
  
  u <- pnorm((log(X1) - mu1)/sigma1)
  v <- pnorm((log(X2) - mu2)/sigma2)
  inv.data <- cbind(u, v)
  starting.params <- c(1.01, 4)
  theta1 = optimise(function(x){(sum(-dCopula(inv.data,  joeCopula(x), log = TRUE)))}, interval = starting.params)$minimum
  return(list(mu1= mu1, mu2=mu2, sigma1=sigma1, sigma2=sigma2, theta=theta1))
}

rjoint <- function(n, mu1, mu2, sigma1, sigma2, theta){
  #In: number of observations
  #In: parameters of the model
  #Out: dataframe with observations
  
  dist <-  mvdc(copula=joeCopula(theta,dim=2), margins=c("lnorm","lnorm"),
                paramMargins=list(list(mean=mu1, sd=sigma1),
                                  list(mean=mu2, sd=sigma2)))
  stot <-  rMvdc(n, dist)
  return(data.frame(PLI = stot[,1], WC = stot[,2]))}


data.estimates <- estimate.params(data.original)

set.seed(18101996)
data.generated <- rjoint(nrow(data.original), mu1= data.estimates$mu1, mu2 = data.estimates$mu2, sigma1 = data.estimates$sigma1, sigma2 = data.estimates$sigma2, theta= data.estimates$theta)

data.original.plot6 <- data.original[,2:3]
data.original.plot6$origin <- "Original"
data.generated$origin <- "Simulated"

data.plot6 <- rbind(data.original.plot6, data.generated)

p1 <- ggplot(data.plot6) + geom_point(aes(PLI, WC, colour=origin))+ xlim(0, 60) + ylim(0, 80)+ theme_classic() + scale_color_colorblind()
p1 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(p1)

# If mus higher then shift along x-axis and y-axis, respectively
set.seed(18101996)
data.generated2 <- rjoint(nrow(data.original), 3, 4, data.estimates$sigma1, data.estimates$sigma2, data.estimates$theta)
data.generated2$origin <- "Simulated"

data.plot62 <- rbind(data.original.plot6, data.generated2)

p2 <- ggplot(data.plot62) + geom_point(aes(PLI, WC, colour=origin)) + xlim(0, 60) + ylim(0, 80)+ theme_classic()+ scale_color_colorblind()

suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                          p2 + theme(legend.position="none"),
                                          nrow=1),
                              mylegend,nrow=2,heights=c(8, 2)))




set.seed(18101996)
data.generated3 <- rjoint(nrow(data.original), data.estimates$mu1, data.estimates$mu2, 1, 0.6, data.estimates$theta)
data.generated3$origin <- "Simulated"

data.plot63 <- rbind(data.original.plot6, data.generated3)


p3 <- ggplot(data.plot63) + geom_point(aes(PLI, WC, colour=origin)) + xlim(0, 60) + ylim(0, 80)+ theme_classic() + scale_color_colorblind()

suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                          p3 + theme(legend.position="none"),
                                          nrow=1),
                              mylegend,nrow=2,heights=c(8, 2)))


set.seed(18101996)
data.generated4 <- rjoint(nrow(data.original), data.estimates$mu1, data.estimates$mu2,  data.estimates$sigma1, data.estimates$sigma2, 2*data.estimates$theta)
data.generated4$origin <- "Simulated"

data.plot64 <- rbind(data.original.plot6, data.generated4)

p4 <- ggplot(data.plot64) + geom_point(aes(PLI, WC, colour=origin)) + xlim(0, 60) + ylim(0, 80)+theme_classic() + scale_color_colorblind()

suppressWarnings(grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                          p4 + theme(legend.position="none"),
                                          nrow=1),
                              mylegend,nrow=2,heights=c(8, 2)))


rm(data.plot6)
rm(data.plot62)
rm(data.plot63)
rm(data.plot64)
rm(data.generated)
rm(data.generated2)
rm(data.generated3)
rm(data.generated4)
rm(p1)
rm(p2)
rm(p3)
rm(p4)

estimate.paramstime <- function(data){
  #In: a dataframe
  #Out: a list with parameters
  
  X1 <- data[,c("PLI")]
  X2 <- data[,c("WC")]
  
  mu1 <- uniroot(estimate.mu, c(min(log(X1)), max(log(X1))), x=X1)$root
  sigma1 <- uniroot(estimate.sigma, c(0.1, 4*sd(log(X1))), x=X1)$root
  mu2 <- uniroot(estimate.mu, c(min(log(X2)), max(log(X2))), x=X2)$root
  sigma2 <- uniroot(estimate.sigma, c(0.1, 4*sd(log(X2))), x=X2)$root
  
  u <- pnorm((log(X1) - mu1)/sigma1)
  v <- pnorm((log(X2) - mu2)/sigma2)
  inv.data <- cbind(u, v)
  starting.params <- c(1.01, 4)
  theta1 = optimise(function(x){(sum(-dCopula(inv.data,  joeCopula(x), log = TRUE)))}, interval = starting.params)$minimum
  
  mu1.time <- system.time(uniroot(estimate.mu, c(min(log(X1)), max(log(X1))), x=X1)$root)[3]
  mu2.time <- system.time(uniroot(estimate.mu, c(min(log(X2)), max(log(X2))), x=X2)$root)[3]
  sigma1.time <- system.time(uniroot(estimate.sigma, c(0.1, 4*sd(log(X1))), x=X1)$root)[3]
  sigma2.time <- system.time(uniroot(estimate.sigma, c(0.1, 4*sd(log(X2))), x=X2)$root)[3]
  theta1.time <- system.time(optimise(function(x){(sum(-dCopula(inv.data,  joeCopula(x), log = TRUE)))}, interval = starting.params)$minimum)[3]
  
  return(list(mu1 = mu1, 
              mu2 = mu2, 
              sigma1 = sigma1, 
              sigma2 = sigma2, 
              theta=theta1, 
              mu1.time = unname(mu1.time), 
              mu2.time = unname(mu2.time), 
              sigma1.time= unname(sigma1.time),
              sigma2.time = unname(sigma2.time), 
              theta.time = unname(theta1.time)))
}
#

# Generate r datasets of size n

# For TA's: USES PARALLEL
simulator <- function(n, r, mu1, mu2, sigma1, sigma2, thetas){
  sim.dat <- replicate(r, rjoint(n, mu1, mu2, sigma1, sigma2, thetas), simplify = FALSE)
  sim.dat.est <- mclapply(sim.dat, estimate.paramstime)
  return(sim.dat.est)
}

r <- 100
n <- c(100, 200, 500, nrow(data.original), 1000)
mu1 <- 1
mu2 <- 3
sigma1 <- 2
sigma2 <- 0.5

thetas <- 2

set.seed(18101996)
q8.res <- mclapply(n, function(z) simulator(z, r, mu1, mu2, sigma1, sigma2, thetas))
q8.res.df <- data.frame(matrix(unlist(q8.res), nrow = (r*length(n)), byrow=T))
q8.res.df$n <- rep(n, each=r)
names(q8.res.df) <- c("mu1", "mu2", "sigma1", "sigma2", "theta", "time.mu1", "time.mu2", "time.sigma1", "time.sigma2", "time.theta", "n")

Indicator <- function(data, price){
  (rowSums(data) > price) * (rowSums(data))
}

rjoint2 <- function(n, mu1, mu2, sigma1, sigma2, theta){
  #In: number of observations
  #In: parameters of the model
  #Out: dataframe with observations
  
  dist <-  mvdc(copula=joeCopula(theta,dim=2), margins=c("lnorm","lnorm"),
                paramMargins=list(list(mean=mu1, sd=sigma1),
                                  list(mean=mu2, sd=sigma2)))
  stot <-  rMvdc(n, dist)
  return(stot)}


mu1 <- data.estimates$mu1
mu2 <- data.estimates$mu2
sigma1 <- data.estimates$sigma1
sigma2 <- data.estimates$sigma2
theta <-  data.estimates$theta


mc.plain <- function(n,B, price){
  datasets = replicate(B, rjoint2(n, mu1, mu2, sigma1, sigma2,theta), simplify=FALSE)
  mclapply(datasets,function(y) mean(Indicator(y, price)))
}

prices <- seq(100, 200, by=10)

mc.estimates <- mclapply(prices, function(z) mc.plain(1e5, 100, z))
mc.estimates.df <- data.frame(matrix(unlist(mc.estimates), ncol=11))


djoint <- function(x, mu1, mu2, sigma1, sigma2, theta){
  #In: a dataset
  #In: parameters of the model
  #Out: densities
  
  dist <-  mvdc(copula=joeCopula(theta,dim=2), margins=c("lnorm","lnorm"),
                paramMargins=list(list(mean=mu1, sd=sigma1),
                                  list(mean=mu2, sd=sigma2)))
  stot <-  dMvdc(x, dist)
  return(stot)
}

mu1.new <- 1.5*mu1

mc.is <- function(n,B, price){
  datasets = replicate(B, rjoint2(n, mu1.new, mu2, sigma1, sigma2, theta), simplify=FALSE)
  weights = mclapply(datasets, function(z) (djoint(z, mu1, mu2, sigma1, sigma2, theta) /
                                              djoint(z, mu1.new, mu2, sigma1, sigma2, theta)))
  indicators = mclapply(datasets, function(y) Indicator(y, price))
  weighted.inds = mcMap('*', indicators, weights)
  mclapply(weighted.inds, mean) 
}
set.seed(1810199)
mcis.estimates <- mclapply(prices, function(z) mc.is(1e5, 1, z), mc.set.seed = 10)


bootstrapper <- function(B, price){
  # Makes a large datalist with replicates of the data
  data <- replicate(B, data.original[sample(nrow(data.original), nrow(data.original), replace=T),2:3] , simplify=FALSE)
  
  # estimates parameters for that data
  data.par <- mclapply(data, function(z) estimate.params(z))
  
  # Generates new data
  data.gen <- mclapply(data.par, function(z) rjoint2(1e5, 1.5*z$mu1, z$mu2, z$sigma1, z$sigma2, z$theta))
  # Applies indicator to new data
  data.ind <- mclapply(data.gen, function(t) Indicator(t, price))
  # Bottleneck: finds the weights
  data.weights <- mcmapply(function(a,b) (djoint(b, a$mu1, a$mu2, a$sigma1, a$sigma2, a$theta) 
                                          /djoint(b, 1.5*a$mu1, a$mu2, a$sigma1, a$sigma2, a$theta)), a=data.par, b=data.gen, SIMPLIFY = FALSE)
  # Makes a weighted sample
  weighted.samp <- mcMap('*', data.ind, data.weights)
  # Calculates VT
  weighted.m <- mclapply(weighted.samp, mean)
  # Finds the quantiles
  quantile(unlist(weighted.m), c(0.2, 0.8))
}


data.boots <- mclapply(prices, function(v) bootstrapper(1e3, v))


q8.res.df.long.est <- gather(q8.res.df[,c(1:5, 11)], variable, estimate, 1:5,)
ggplot(q8.res.df.long.est %>% group_by(n, variable)) + geom_boxplot(aes(x=n, group=n, y=estimate)) + facet_wrap(~ variable) + theme_classic()

q8.res.df.long <- gather(q8.res.df[,6:11], variable, estimate, 1:5,)

times <- suppressWarnings(q8.res.df.long %>% group_by(n, variable) %>% summarise(mean.time = mean(estimate), .groups = "keep"))
suppressWarnings(ggplot(times) + geom_point(aes(x=n, group=n, y=mean.time)) + facet_wrap(~ variable) + theme_classic())

RMSE <- suppressWarnings(q8.res.df %>% dplyr::group_by(n) %>% summarise(RMSEmu1 = mean((mu1 - 1)^2), 
                                                                        RMSEmu2 = mean((mu2 - 3)^2), 
                                                                        RMSEsigma1 = mean((sigma1 - 2)^2), 
                                                                        RMSEsigma2 = mean((sigma2 - 0.5)^2), 
                                                                        RMSEtheta = mean((theta - 2)^2), .groups = "keep") )

RMSE.long <- gather(RMSE, variable, rmse, 2:6)
ggplot(RMSE.long) + geom_path(aes(n, sqrt(rmse))) + facet_wrap( ~ variable) + theme_classic() + ylab("RMSE")


price <- function(t){
  40000*exp(-t/7)
}

mc.estimates.long <- gather(mc.estimates.df, simulation, value, 1:11)
mc.estimates.long$simulation <- rep(prices, each=100)
total.picture <- mc.estimates.long %>% group_by(simulation) %>%  summarise(value=mean(value), .groups = "keep")
total.picture$func <- "V(t)"

prices.df <- data.frame(simulation=prices, value=price(prices))
prices.df$func <- "P(t)"
tot.pic <- rbind(total.picture, prices.df)

ggplot(tot.pic) + geom_path(aes(simulation, value, color=func)) + theme_classic() + scale_color_colorblind()