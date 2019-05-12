## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
	fig.align = "center",
	fig.width = 6, fig.height = 6,
	cache = FALSE,
	collapse = TRUE,
	comment = "#>",
	highlight = TRUE
)

## ----Willow-tit, echo=FALSE, out.width="\\textwidth", fig.cap="(ref:cap-Willow-tit)"----
knitr::include_graphics("figures/Poecile_montanus.jpg")

## ----libraries-----------------------------------------------------------
# Load libraries
library(sp)
library(raster)
library(hSDM)

## ----Kery2010-data-------------------------------------------------------
# Load Kéry et al. 2010 data
data(data.Kery2010, package="hSDM")
head(data.Kery2010)

## ----normalizing---------------------------------------------------------
# Normalized variables
elev.mean <- mean(data.Kery2010$elevation)
elev.sd <- sd(data.Kery2010$elevation)
juldate.mean <- mean(c(data.Kery2010$juldate1,
                     data.Kery2010$juldate2,
                     data.Kery2010$juldate3),na.rm=TRUE)
juldate.sd <- sd(c(data.Kery2010$juldate1,
                 data.Kery2010$juldate2,
                 data.Kery2010$juldate3),na.rm=TRUE)
data.Kery2010$elevation <- (data.Kery2010$elevation-elev.mean)/elev.sd
data.Kery2010$juldate1 <- (data.Kery2010$juldate1-juldate.mean)/juldate.sd
data.Kery2010$juldate2 <- (data.Kery2010$juldate2-juldate.mean)/juldate.sd
data.Kery2010$juldate3 <- (data.Kery2010$juldate3-juldate.mean)/juldate.sd

## ----sites, out.width="\\textwidth", fig.cap="(ref:cap-sites)", fig.width=9, fig.height=6----
# Landscape and observation sites
sites.sp <- SpatialPointsDataFrame(coords=data.Kery2010[c("coordx","coordy")],
                                   data=data.Kery2010[,-c(1,2)])
xmin <- min(data.Kery2010$coordx)
xmax <- max(data.Kery2010$coordx)
ymin <- min(data.Kery2010$coordy)
ymax <- max(data.Kery2010$coordy)
res_spatial_cell <- 10
ncol <- ceiling((xmax-xmin)/res_spatial_cell)
nrow <- ceiling((ymax-ymin)/res_spatial_cell)
Xmax <- xmin+ncol*res_spatial_cell
Ymax <- ymin+nrow*res_spatial_cell
ext <- extent(c(xmin, Xmax,
                ymin, Ymax))
landscape <- raster(ncols=ncol,nrows=nrow,ext)
values(landscape) <- runif(ncell(landscape),0,1)
landscape.po <- rasterToPolygons(landscape)
par(mar=c(0.1,0.1,0.1,0.1))
plot(landscape.po)
points(sites.sp, pch=16, cex=1, col="black")

## ----Neighborhood--------------------------------------------------------
# Neighborhood
# Rasters must be projected to correctly compute the neighborhood
crs(landscape) <- '+proj=utm +zone=1'
# Cell for each site
cells <- extract(landscape,sites.sp,cell=TRUE)[,1]
# Neighborhood matrix
ncells <- ncell(landscape)
neighbors.mat <- adjacent(landscape, cells=c(1:ncells), directions=8,
                          pairs=TRUE, sorted=TRUE)
# Number of neighbors by cell
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
# Adjacent cells
adj <- neighbors.mat[,2]

## ----Arranging-data------------------------------------------------------
# Arranging data
# data.obs
nsite <- length(data.Kery2010$coordx)
count <- c(data.Kery2010$count1,data.Kery2010$count2,data.Kery2010$count3)
juldate <- c(data.Kery2010$juldate1,data.Kery2010$juldate2,
             data.Kery2010$juldate3)
site <- rep(1:nsite,3)
data.obs <- data.frame(count,juldate,site)
data.obs <- data.obs[!is.na(data.obs$juldate),]
# data.suit
data.suit <- data.Kery2010[c("coordx","coordy","elevation")]
data.suit$cells <- cells
data.suit <- data.suit[-139,] # Removing site 139 with no juldate

## ----Kery2010-pois-mod, cache=TRUE---------------------------------------
# hSDM.poisson
data.pois <- data.obs
data.pois$elevation <- data.suit$elevation[as.numeric(as.factor(data.obs$site))]
mod.Kery2010.pois <- hSDM.poisson(counts=data.pois$count,
                                  suitability=~elevation+I(elevation^2),
                                  data=data.pois, beta.start=0)

## ----Kery2010-pois-mcmc, results="markup"--------------------------------
# Outputs
summary(mod.Kery2010.pois$mcmc)

## ----Kery2010-pois-pred--------------------------------------------------
# Predictions
npred <- 100
nsamp <- dim(mod.Kery2010.pois$mcmc)[1]
# Abundance-elevation
elev.seq <- seq(500,3000,length.out=npred)
elev.seq.n <- (elev.seq-elev.mean)/elev.sd
beta <- as.matrix(mod.Kery2010.pois$mcmc[,1:3])
tbeta <- t(beta)
X <- matrix(c(rep(1,npred),elev.seq.n,elev.seq.n^2),ncol=3)
N <- matrix(NA,nrow=nsamp,ncol=npred)
for (i in 1:npred) {
    N[,i] <- exp(X[i,] %*% tbeta)
}
N.est.pois <- apply(N,2,mean)
N.q1.pois <- apply(N,2,quantile,0.025)
N.q2.pois <- apply(N,2,quantile,0.975)

## ----Kery2010-Nmix-mod, cache=TRUE---------------------------------------
# hSDM.Nmixture
mod.Kery2010.Nmix <- hSDM.Nmixture(# Observations
                     counts=data.obs$count,
                     observability=~juldate+I(juldate^2),
                     site=data.obs$site,
                     data.observability=data.obs,
                     # Habitat
                     suitability=~elevation+I(elevation^2),
                     data.suitability=data.suit,
                     # Predictions
                     suitability.pred=NULL,
                     # Chains
                     burnin=10000, mcmc=5000, thin=5,
                     # Starting values
                     beta.start=0,
                     gamma.start=0,
                     # Priors
                     mubeta=0, Vbeta=1.0E6,
                     mugamma=0, Vgamma=1.0E6,
                     # Various
                     seed=1234, verbose=1,
                     save.p=0, save.N=0)

## ----Kery2010-Nmix-mcmc, results="markup"--------------------------------
# Outputs
summary(mod.Kery2010.Nmix$mcmc)

## ----Kery2010-Nmix-pred--------------------------------------------------
# Predictions
nsamp <- dim(mod.Kery2010.Nmix$mcmc)[1]
# Abundance-elevation
beta <- as.matrix(mod.Kery2010.Nmix$mcmc[,1:3])
tbeta <- t(beta)
N <- matrix(NA,nrow=nsamp,ncol=npred)
for (i in 1:npred) {
    N[,i] <- exp(X[i,] %*% tbeta)
}
N.est.Nmix <- apply(N,2,mean)
N.q1.Nmix <- apply(N,2,quantile,0.025)
N.q2.Nmix <- apply(N,2,quantile,0.975)
# Detection-Julian date
juldate.seq <- seq(100,200,length.out=npred)
juldate.seq.n <- (juldate.seq-juldate.mean)/juldate.sd
gamma <- as.matrix(mod.Kery2010.Nmix$mcmc[,4:6])
tgamma <- t(gamma)
W <- matrix(c(rep(1,npred),juldate.seq.n,juldate.seq.n^2),ncol=3)
delta <- matrix(NA,nrow=nsamp,ncol=npred)
for (i in 1:npred) {
    delta[,i] <- inv.logit(X[i,] %*% tgamma)
}
delta.est.Nmix <- apply(delta,2,mean)
delta.q1.Nmix <- apply(delta,2,quantile,0.025)
delta.q2.Nmix <- apply(delta,2,quantile,0.975)

## ----Kery2010-Nmix-iCAR-mod, cache=TRUE----------------------------------
# hSDM.Nmixture.iCAR
mod.Kery2010.Nmix.iCAR <- hSDM.Nmixture.iCAR(# Observations
                            counts=data.obs$count,
                            observability=~juldate+I(juldate^2),
                            site=data.obs$site,
                            data.observability=data.obs,
                            # Habitat
                            suitability=~elevation+I(elevation^2),
                            data.suitability=data.suit,
                            # Spatial structure
                            spatial.entity=data.suit$cells,
                            n.neighbors=n.neighbors, neighbors=adj,
                            # Chains
                            burnin=20000, mcmc=10000, thin=10,
                            # Starting values
                            beta.start=0,
                            gamma.start=0,
                            Vrho.start=1,
                            # Priors
                            mubeta=0, Vbeta=1.0E6,
                            mugamma=0, Vgamma=1.0E6,
                            priorVrho="1/Gamma",
                            shape=1, rate=1,
                            # Various
                            seed=1234, verbose=1,
                            save.rho=0, save.p=0, save.N=0)

## ----Kery2010-Nmix-iCAR-mcmc,results="markup"----------------------------
summary(mod.Kery2010.Nmix.iCAR$mcmc)

## ----Kery2010-Nmix-iCAR-spatial-effects, out.width="\\textwidth", fig.width=9, fig.height=6, fig.cap="(ref:cap-spatial)"----
# Spatial random effects
rho.pred <- mod.Kery2010.Nmix.iCAR$rho.pred
r.rho.pred <- rasterFromXYZ(cbind(coordinates(landscape),rho.pred))
plot(r.rho.pred)
# Mean abundance by site
ma <- apply(sites.sp@data[,3:5],1,mean,na.rm=TRUE)
points(sites.sp,pch=16,cex=0.5)
points(sites.sp,pch=1,cex=ma/2)

## ----Kery2010-Nmix-iCAR-pred---------------------------------------------
# Predictions
nsamp <- dim(mod.Kery2010.Nmix.iCAR$mcmc)[1]
# Abundance-elevation
beta <- as.matrix(mod.Kery2010.Nmix.iCAR$mcmc[,1:3])
tbeta <- t(beta)
N <- matrix(NA,nrow=nsamp,ncol=npred)
# Simplified way of obtaining samples for rho
rho.samp <- sample(rho.pred,nsamp,replace=TRUE)
for (i in 1:npred) {
    N[,i] <- exp(X[i,] %*% tbeta + rho.samp)
}
N.est.Nmix.iCAR <- apply(N,2,mean)
N.q1.Nmix.iCAR <- apply(N,2,quantile,0.025)
N.q2.Nmix.iCAR <- apply(N,2,quantile,0.975)

# Detection-Julian date
gamma <- as.matrix(mod.Kery2010.Nmix.iCAR$mcmc[,4:6])
tgamma <- t(gamma)
delta <- matrix(NA,nrow=nsamp,ncol=npred)
for (i in 1:npred) {
    delta[,i] <- inv.logit(X[i,] %*% tgamma)
}
delta.est.Nmix.iCAR <- apply(delta,2,mean)
delta.q1.Nmix.iCAR <- apply(delta,2,quantile,0.025)
delta.q2.Nmix.iCAR <- apply(delta,2,quantile,0.975)

## ----Kery2010-comp-abundance, fig.cap="(ref:cap-abundance)"--------------
# Expected abundance - Elevation
par(mar=c(4,4,1,1),cex=1.4,tcl=+0.5)
plot(elev.seq,N.est.pois,type="l",
     xlim=c(500,3000),
     ylim=c(0,10),
     lwd=2,
     xlab="Elevation (m a.s.l.)",
     ylab="Expected abundance",
     axes=FALSE)
#lines(elev.seq,N.q1.pois,lty=3,lwd=1)
#lines(elev.seq,N.q2.pois,lty=3,lwd=1)
axis(1,at=seq(500,3000,by=500),labels=seq(500,3000,by=500))
axis(2,at=seq(0,10,by=2),labels=seq(0,10,by=2))
# Nmix
lines(elev.seq,N.est.Nmix,lwd=2,col="red")
#lines(elev.seq,N.q1.Nmix,lty=3,lwd=1,col="red")
#lines(elev.seq,N.q2.Nmix,lty=3,lwd=1,col="red")
# Nmix.iCAR
lines(elev.seq,N.est.Nmix.iCAR,lwd=2,col="dark green")
#lines(elev.seq,N.q1.Nmix.iCAR,lty=3,lwd=1,col="dark green")
#lines(elev.seq,N.q2.Nmix.iCAR,lty=3,lwd=1,col="dark green")

## ----Kery2010-comp-detection, fig.cap="(ref:cap-detection)"--------------
# Detection probability - Julian date
par(mar=c(4,4,1,1),cex=1.4,tcl=+0.5)
plot(juldate.seq,delta.est.Nmix,type="l",
     xlim=c(100,200),
     ylim=c(0,1),
     lwd=2,
     col="red",
     xlab="Julian date",
     ylab="Detection probability",
     axes=FALSE)
lines(juldate.seq,delta.q1.Nmix,lty=3,lwd=1,col="red")
lines(juldate.seq,delta.q2.Nmix,lty=3,lwd=1,col="red")
axis(1,at=seq(100,200,by=20),labels=seq(100,200,by=20))
axis(2,at=seq(0,1,by=0.2),labels=seq(0,1,by=0.2))
# Nmix.iCAR
lines(juldate.seq,delta.est.Nmix.iCAR,lwd=2,col="dark green")
lines(juldate.seq,delta.q1.Nmix.iCAR,lty=3,lwd=1,col="dark green")
lines(juldate.seq,delta.q2.Nmix.iCAR,lty=3,lwd=1,col="dark green")

