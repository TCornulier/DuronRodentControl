### Duron et al. Neobiota ###
### Simulated scenarios of control effort ###
### November 2020 ###

library("sp")

# function to capture probability per day at trap location
p.trap<- function(trap.x, trap.y, HR.x, HR.y, HR.sigma, g0){
  p0<- dnorm(0, 0, HR.sigma)^2
  1 - prod(1 - dnorm(trap.x, HR.x, HR.sigma) * dnorm(trap.y, HR.y, HR.sigma) * g0 / p0)
}

# function to simulate trapping grid
grid.gen<- function(dist.within.transect, dist.between.transect) expand.grid(x= seq(0, 2000, by= dist.within.transect), y= seq(0, 1000, by= dist.between.transect))
trap.grid<- grid.gen(25, 90)

plot(trap.grid, asp= 1)

# function to compute capture probability per rat on a grid at day 1, 2, ... , 15
pcapt.ses<- function(dist.traps, dist.transects, sigma, g0, nsim= 10000, max.days= 15){
  trap.grid<- grid.gen(dist.traps, dist.transects)
  rat.x<- runif(nsim, 0, 2000)
  rat.y<- runif(nsim, 0, 1000)
  p.rat<- matrix(NA, nsim, max.days+1)
  for(i in 1:nsim){
    p.rat[i, 1:(max.days+1)]<- 1 - (1 - p.trap(trap.grid$x, trap.grid$y, HR.x= rat.x[i], HR.y= rat.y[i], HR.sigma= sigma, g0= g0))^(c(100, 1:max.days)) # cumulative probability of rat being caught on day t given location of home range
  }
  out<- apply(p.rat, 2, mean)
  names(out)<- paste("D", c(100, 1:max.days), sep= "")
  out
}

system.time(pcapt.ses(25, 25, 14.6, 0.09, 10000, 15))

# data frame for simulation parameters
simulation1<- expand.grid(dist.transect= seq(25, 100, by= 25), dist.trap= seq(25, 100, by= 25))
simulation1<- simulation1[!(simulation1$dist.trap > simulation1$dist.transect), ]
simulation1<- rbind(c(15, 15), simulation1)
row.names(simulation1)<- 1:nrow(simulation1)
simulation1<- simulation1[-11, ]

# nb hours required for 1 complete coverage of the trapping grid
simulation1$hours<- round((((1000*(2000/simulation1$dist.transect)+simulation1$dist.transect*((2000/simulation1$dist.transect))-1)*10)/100)/60+(((1000/simulation1$dist.trap)*(2000/simulation1$dist.transect))*37)/3600, 1) #37sec is the arbitrary average time spent per trap for a given day over the entire session
simulation1$numb.pers.day<- round(simulation1$hours/4, 1) # max working hours per day = 4

# simulation1
system.time(simulation1<- cbind(simulation1, t(apply(simulation1, 1, function(x){pcapt.ses(x[1], x[2], sigma= 14.6, g0= 0.09, nsim= 10000, max.days= 30)}))))

simulation1[, 1:6]

# plot cumulative capture probability over time (15 days)
matplot(t(as.matrix(simulation1[, c(6:20)])), type= "b", ylab= "Cumulative capture probability", xlab= "Trap nights", main= "Cumulative rat capture probability over time") #grid numbers indicated on plot

matplot(t(as.matrix(simulation1[, c(6:20)])), type= "l", lty= 1, lwd= 2, bpy.colors(10), ylab= "Cumulative capture probability", xlab= "Trap nights")
abline(h= 0.5, lty= 3)
abline(h= 0.8, lty= 5)

# plot cumulative capture probability over time (100 days)
matplot(t(as.matrix(simulation1[, c(6:20, 5)])), type= "b", ylab= "Cumulative capture probability", xlab= "Trap nights", main= "Cumulative rat capture probability over time")

# calculate difference in capture proba between successive days
apply(simulation1[, c(6:20)], 1, diff)

# plot difference over time
matplot(apply(simulation1[, c(6:20)], 1, diff), type= "b", ylab= "Difference in capture probability", xlab= "Days", main= "Difference in rat capture probability between successive days")

gain<- simulation1[, c(6:35)]
expenditure<- outer(simulation1$hours, 1:15) + simulation1$hours # add 1 extra day for initial baiting (extra baiting day not included in day count - columns represent trapping days only)


#### showing what is achievable under sets of social/financial constraints 
## maximum budget
wage.per.hour<- 10
nb.people<- 10
hour.per.day.person<- 4
nb.day.max<- 30
provisioning.cost<- 6000 # includes 1000 euros for food (10 pers/15 days) and one helicopter trip of 5000 euros
provisioning.recur<- 15 # provisioning interval in days

nb.pers.day<- round(simulation1$hours/4, 1) # nb people required to complete grid in 1 day
splits<- ceiling(nb.pers.day / nb.people) # nb splits required to complete grid with #nb.people
nb.people.for.last.split<- ceiling(nb.pers.day %% nb.people) # nb people required to complete the last split (modulus)

baiting.days<- splits

session.length.max<- floor((nb.day.max - baiting.days) / splits)

cum.cost<- matrix(NA, 10, max(session.length.max))
cum.capt<- matrix(NA, 10, max(session.length.max))

for(i in 1:10){
    for(j in 1:max(session.length.max)){
        npd<- nb.people * (splits[i] - 1) + nb.people.for.last.split[i] # nb people per unit session length
        prov.cost<- provisioning.cost * ceiling((baiting.days[i] + splits[i] * j) / 10)
        cum.cost[i, j]<- wage.per.hour * hour.per.day.person * npd * (j + 1) + prov.cost
        cum.capt[i, j]<- gain[i, j]
        if((baiting.days[i] + splits[i] * j) > nb.day.max){
            cum.cost[i, j]<- NA
            cum.capt[i, j]<- NA
        }
    }
}

# Fig.4
png("./Fig.4.png", 3600, 1600, type = "cairo-png", pointsize= 60)
par(mfrow= c(1, 2), mar= c(4.1, 4.1, 3.1, 1))
cex.lwd<- 2.5
cex.axis<- 1
cex.lab<- 1.1
matplot(t(as.matrix(simulation1[, c(6:20)])), type= "l", lty= 1, col= bpy.colors(10), ylab= "Cumulative capture probability", xlab= "Trap nights", cex.lab= cex.lab, cex.axis= cex.axis, lwd= cex.lwd)
matplot(t(cum.cost), t(cum.capt), col= bpy.colors(10), ylim= c(0,1), ylab= "Cumulative capture probability", xlab= "Cumulative cost in euros", type= "l", lty= 1, lwd= cex.lwd, cex.lab= cex.lab, cex.axis= cex.axis)
legend("topleft", inset=.03, cex= 0.8, lty= 1, ncol= 5, x.intersp = 0.5, legend=row.names(simulation1), col= bpy.colors(10), lwd= cex.lwd)
abline(h= 0.5, lty= 3, lwd= cex.lwd)
abline(h= 0.8, lty= 5, lwd= cex.lwd)

# add points on plot to indicate what is achievable in 15 days and at what cost
simulation1$cum.cost.15d<- c(NA, 17160, 16400, 16200, 15360, 17600, 15920, 15080, 15640, 18400) # cumulative costs calculated as cost of trapping n occasions; n = (15 days / (total number of spilts - 1 baiting split); e.g. grid 2: 15/(5-1)=2.5 rounded to 2 full occasions
simulation1$cum.capt15d<- c(NA, 0.32529284, 0.32165323, 0.29971891, 0.23025682, 0.25187409, 0.16966158, 0.12785644, 0.11547691, 0.1785705) # cumulative capture proba achieved after 15 days
points(simulation1$cum.cost.15d, simulation1$cum.capt15d, col= bpy.colors(10), pch= 16, cex= 1)
dev.off()