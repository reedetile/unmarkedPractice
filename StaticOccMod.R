#Description-----------------------------------------
#Practice of a single season occupancy model in unmarked
#  06 Nov 2024
#RCS

#Initialize -----------------------------------------
library(unmarked)

# Load functions--------------------------------------


# Global Variables-------------------------------------
getwd()
wt <- read.csv(system.file("csv","widewt.csv", package="unmarked")) #full data set

# Program Body------------------------------------------

#prepping data in an unmarked frame

PA <- wt[,2:4] #just the presence absence of the spp
siteCovs <- wt[,c('elev','forest','length')] #pulling out site lvl covariates
obsCovs <- list(date = wt[,c('date.1','date.2','date.3')],
                ivel = wt[,c('ivel.1','ivel.2','ivel.3')]) #pulling out survey lvl covars
data <- unmarkedFrameOccu(y = PA,
                          siteCovs = siteCovs,
                          obsCovs = obsCovs)

summary(data)
obsCovs(data) <- scale(obsCovs(data)) #remember that the data needs to be scaled!
### fitting models ###
# Detection covariates follow first tilde, then occupancy covariates
fm1 <- occu(~1 ~1, data = data)
fm2 <- occu(~elev ~1, data = data)
fm3 <- occu(~1 ~ elev, data = data)
fm4 <- occu(~date ~1, data = data)


#back transforming parameter estimates
backTransform(fm2, 'state')
backTransform(fm2, 'det') #won't work because covariates are present

#predict might be better
newData <- data.frame(ivel = 0, date = seq(from = -1.8, to = 2.2, by  = 0.1))
output <- predict(fm2, type = 'det',newData = newData, appendData = T)

#get confidence intervals
confint(fm2, type = 'det')
confint(fm2, type = 'det', method = 'profile') #not sure what profiling does


#Model selection + model fit
fms <- fitList('psi(.)p(.)' = fm1,
               'psi(.)p(date)' = fm2)

modSel(fms) #date strongly affects detection!

#using a parametric bootstrap
chisq <- function(fm) {
  umf <- fm@data
  y <- umf@y
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}

(pb <- parboot(fm2, statistic=chisq, nsim=100, parallel=FALSE))

ranef(fm2) #does this provide the conditional dist of occurrence at each site?
re <- ranef(fm2)
#best unbiarsed predictor
EBUP <- bup(re, stat = "mode")
EBUP #okay so this is the best estimated best predictor at each site?

sum(EBUP)/numSites(data)
