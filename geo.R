library("geoR")
library(ggplot2)
# load the data
da = read.table("arsenic.txt",header=TRUE)
# There seems to be some duplicate coordinates. So we add a little jitter
newcoords = jitterDupCoords(cbind(da$latitude,da$longitude),max=1e-4) 
da$latitude = newcoords[,1]
da$longitude = newcoords[,2]

pdf("arsenic.pdf",width = 5,height = 5)
ggplot(data=da,aes(x=longitude,latitude,color=log(arsenic)))+geom_point(size=0.5)+ scale_color_gradient2(low = "blue", midpoint = median(log(da$arsenic)), high = "red")+theme_bw()
dev.off()
pdf("arsenic1.pdf",width = 5,height = 5)
ggplot(data=da,aes(x=longitude,latitude,color=log(arsenic)))+geom_point(size=0.5)+theme_bw()
dev.off()
# make the data from a frame to a geodata object.
arsenic = as.geodata(da,covar.col = 4:7)


d = dist(newcoords);
maxdist = max(d)

# Compute the variogram under isotropic assumption
vg = variog(arsenic,max.dist = maxdist/2,option="bin",breaks = seq(0.01,maxdist/2,by=0.02),lambda=0)
#lambda = 0 means log transformation of the arsenic concentration.
plot(vg)
# The variogram seems to have a sill around 6, nugget around 2. 
# But the range is not clear. Let's say 3.5


# Check for anisotropy. Compute the variogram along 4 directions
vg4 = variog4(arsenic,option = "bin",max.dist = maxdist/2,lambda=0)
plot(vg4)


# fit Matern variogram to the data: partial sill estimate = 6-2 = 4.
# weights = "cressie" uses the wighted least squares. weights = "npair" means OLS.
fitV = variofit(vg,ini.cov.pars = c(4,3.5),cov.model = "matern",
                fix.nugget = FALSE,nugget = 2,kappa = 0.5,fix.kappa = FALSE,
                weights = "cressie")

sill = fitV$nugget + fitV$cov.pars[1]

fitV$nugget
fitV$cov.pars

fitV.ols = variofit(vg,ini.cov.pars = c(4,3.5),cov.model = "matern",
                fix.nugget = FALSE,nugget = 2,kappa = 0.5,fix.kappa = FALSE,
                weights = "npairs")
fitV.ols$nugget 
fitV.ols$cov.pars

# compute the covariance functions at the required lags given in vg$uvec
fit.cov = cov.spatial(obj= vg$uvec,cov.model=fitV$cov.model, cov.pars = fitV$cov.pars,kappa = fitV$kappa)

# compute teh variogram
fit.variog = fitV$cov.pars[1] - fit.cov + fitV$nugget

residuals = sqrt(vg$n)*(vg$v - fit.variog)/fit.variog 
  
plot(vg$u,residuals,xlab="lag distance")
  
  
# compute teh variogram for ols
fit.cov.ols = cov.spatial(obj= vg$uvec,cov.model=fitV.ols$cov.model, cov.pars = fitV.ols$cov.pars,kappa = fitV.ols$kappa)


fit.variog.ols = fitV.ols$cov.pars[1] - fit.cov.ols + fitV.ols$nugget

residuals.ols = sqrt(vg$n)*(vg$v - fit.variog.ols)/fit.variog.ols

plot(vg$u,residuals.ols,xlab="lag distance")

### likelihood 
res.lik <- likfit(geodata = arsenic,fix.nugget = FALSE,nugget = 3.5, fix.kappa=FALSE,lambda = 0,cov.model = "matern",ini.cov.pars = c(2,0.5),limits =pars.limits(phi=c(lower=0,upper=7)))

intermat <- read.table("interior500x300.txt")
dim(intermat)

fun1 <- function(ind)
{
  x <- which(intermat[ind,]==1)
  if(length(x) == 0){
    return(NULL)
  }else{
    return(cbind(ind,x))
  }
}

locind <- do.call("rbind",lapply(1:nrow(intermat),fun1))
loc <- data.frame(latitude = seq(20,27,l=500)[locind[,1]],longitude = seq(88,93,l=300)[locind[,2]])

kg = krige.conv(arsenic, loc = loc, krige = krige.control(obj.m = res.lik,lambda=0))
save(res.lik,kg,file="lik_krige.RData")

df.pred <- data.frame(cbind(loc),arsenic =kg$predict)

#krige.conv(arsenic, loc = loc[c(1,2000),], krige = krige.control(obj.model = res.lik))

pdf("arsenic_map.pdf",width = 7,height = 5)
ggplot(data = df.pred,aes(x= latitude,y=longitude,colour=log(arsenic)))+geom_point()+theme_bw()
dev.off()

