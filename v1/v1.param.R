#parameters and settings;
#seed
ran.seed <- 1
#Number of new infections per day within the whole ward 58 & 59 is Poisson distributed;
# lambda <- 100 #lambda for the peak season, now assume no seasonality;

n.pop <- 160000
n.household <-30000 #ask the student for details;
n.person.per.house <- 16/3
n.house.per.latrine <- 10
n.latrine <- n.pop/n.house.per.latrine/n.person.per.house
n.days <- 1000

#parameters of simulation sewage network and latrine points
num.lines <- 1000 #The sewage lines will be set to 1000, 2500, 5000.
num.points <- n.latrine
total.vol <- 2e+10 #mL from Peter and Andrew's rough calculation;
LLOD.test <- 2 #number 10^3 of DNA per 500 mL;
LLOD.uf.test <- 0.2 #number 10^3 of DNA per 5L ultrafiltration samples;
# LLOD.test <- 2/1000
# LLOD.uf.test <- 0.2/1000


dir.store <- paste0("./v1.SimIterative",Sys.Date(),".ssn")

#parameters of transmission;
mu.shed <- 10^8
sigma.shed <- 1
n.days.shed <- 14

#decay and lost parameters;
gamma.shape <- 2.5
gamma.rate <- 0.2

#pooling samples;
n.sample <- 10
n.sample.pooling <- 5
percent.x <- 0.02
percent.y <- 0.02
n.minimum.latrine <- 10

#adaptive sampling settings;
n.week.per.update <- 12
n.drop.add <- 1
n.update <- 100
n.site <- n.latrine
