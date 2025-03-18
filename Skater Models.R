library(spdep)
library(surveyPrev)
library(SUMMER)
library(tidyverse)
library(INLA)
library(ggpubr)
library(R.utils)
source("/Users/alanamcgovern/Desktop/Research/my_helpers.R")

country <- c('Nigeria','Zambia')[1]

# load polygons COUNTRY SPECIFIC ----
# run script to load polygon with correct names
if(country=='Zambia'){
  poly.path <- paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_ZMB_shp")
  poly.adm2 <- st_read(dsn = poly.path, layer = "gadm41_ZMB_2", options = "ENCODING=UTF-8")
}else if(country=='Nigeria'){
  poly.path <- paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_NGA_shp")
  poly.adm2 <- st_read(dsn = poly.path, layer = "gadm41_NGA_2", options = "ENCODING=UTF-8")
  poly.adm2 <- poly.adm2[poly.adm2$ENGTYPE_2=='Local Authority',]
}else{
  stop('No country info available')
}

poly.adm2$admin2 <- 1:nrow(poly.adm2)
poly.adm2$admin2.char <- paste0('admin2_',1:nrow(poly.adm2))
n_admin2 <- nrow(poly.adm2)
poly.adm2$admin1.name <- poly.adm2$NAME_1
poly.adm2$admin2.name <- poly.adm2$NAME_2
admin2.mat <- nb2mat(poly2nb(poly.adm2), zero.policy = TRUE)
colnames(admin2.mat) <- rownames(admin2.mat) <- poly.adm2$admin2.char

poly.adm1 <- poly.adm2 %>% group_by(admin1.name) %>% summarise(geometry = st_union(geometry))
poly.adm1$admin1 <- 1:nrow(poly.adm1)
poly.adm1$admin1.char <- paste0('admin1_',1:nrow(poly.adm1))
poly.adm1 <- sf::st_as_sf(poly.adm1)
n_admin1 <- nrow(poly.adm1)
admin1.mat <- nb2mat(poly2nb(poly.adm1), zero.policy = TRUE)
colnames(admin1.mat) <- rownames(admin1.mat) <- poly.adm1$admin1.char

admin.key <- as.data.frame(poly.adm2[,c('admin2','admin2.char','admin2.name','admin1.name')])
admin.key <- merge(admin.key,poly.adm1,by='admin1.name') %>% arrange(admin2)
admin.key <- admin.key %>% dplyr::select(-geometry.x,-geometry.y)

# load survey data and cluster info w surveyPrev COUNTRY SPECIFIC -----
#  geo <- getDHSgeo(country='Zambia',year=2013)
# cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, by.adm1 = "admin1.name",by.adm2 = "admin2.name")
# cluster.info$data <- left_join(cluster.info$data,poly.adm2[c('admin1.name','admin2.name','admin2','admin2.char')],
#                                by=c('admin1.name','admin2.name'))
# cluster.info$data <- left_join(cluster.info$data,poly.adm1[c('admin1.name','admin1','admin1.char')],by='admin1.name') %>% dplyr::select(-c(geometry.x,geometry.y,geometry))
# 
# dat.path <- getDHSdata(country = 'Nigeria',indicator='nmr',year=2013)
# dat.tmp <- getDHSindicator(dat.path,indicator='nmr',nmr.year = 10)
# 
# dat.tmp <- dat.tmp[(dat.tmp$cluster %in% cluster.info$data$cluster),]
# dat <- dat.tmp[!is.na(dat.tmp$value),]
# dat <- as.data.frame(merge(dat,cluster.info$data))
# 
# clusterdt.prev<- dat %>% group_by(admin1,admin2,v022,weight,cluster) %>% summarise(value=sum(value),N=n())
# 
# save(clusterdt.prev,file="/Users/alanamcgovern/Desktop/Research/Regionalization/ZMB2013NMR_clusterdt.rda")

# geo <- getDHSgeo(country='Zambia',year=2018)
# cluster.info <- clusterInfo(geo=geo, poly.adm1=poly.adm1, poly.adm2=poly.adm2, by.adm1 = "admin1.name",by.adm2 = "admin2.name")
# cluster.info$data <- left_join(cluster.info$data,poly.adm2[c('admin1.name','admin2.name','admin2','admin2.char')],
#                                by=c('admin1.name','admin2.name'))
# cluster.info$data <- left_join(cluster.info$data,poly.adm1[c('admin1.name','admin1','admin1.char')],by='admin1.name') %>% dplyr::select(-c(geometry.x,geometry.y,geometry))
# 
# dat.path <- getDHSdata(country = 'Nigeria',indicator='nmr',year=2018)
# dat.tmp <- getDHSindicator(dat.path,indicator='nmr',nmr.year = 10)
# 
# dat.tmp <- dat.tmp[(dat.tmp$cluster %in% cluster.info$data$cluster),]
# dat <- dat.tmp[!is.na(dat.tmp$value),]
# dat <- as.data.frame(merge(dat,cluster.info$data))
# 
# clusterdt<- dat %>% group_by(admin1,admin2,v022,weight,cluster) %>% summarise(value=sum(value),N=n())
# 
# save(clusterdt,file="/Users/alanamcgovern/Desktop/Research/Regionalization/ZMB2018NMR_clusterdt.rda")

if(country=='Zambia'){
  load(file="/Users/alanamcgovern/Desktop/Research/Regionalization/ZMB2013NMR_clusterdt.rda")
  load(file="/Users/alanamcgovern/Desktop/Research/Regionalization/ZMB2018NMR_clusterdt.rda")
  
}else if(country=='Nigeria'){
  load(file="/Users/alanamcgovern/Desktop/Research/Regionalization/NGA2013NMR_clusterdt.rda")
  load(file="/Users/alanamcgovern/Desktop/Research/Regionalization/NGA2018NMR_clusterdt.rda")
  
}else{
  stop('No country info available')
}


# load admin weights ---------

if(country=='Nigeria'){
  load(file='/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Nigeria/worldpop/adm2_weights_u1.rda')
}

weight.adm2.u1 <- weight.adm2.u1[weight.adm2.u1$year==2017,]
weight.adm2.u1 <- merge(weight.adm2.u1,admin.key)
weight.adm2.u1 <- merge(weight.adm2.u1,
                        weight.adm2.u1 %>% group_by(admin1.char) %>% summarise(admin1_prop = sum(proportion))) %>%
  arrange(admin2)
weight.adm2.u1$within_prop <- weight.adm2.u1$proportion/weight.adm2.u1$admin1_prop

# get admin1 HT estimates ----
options(survey.lonely.psu = "adjust")
CI = 0.95
my.svydesign <- survey::svydesign(ids = ~cluster, 
                                  strata = ~v022, nest = T, weights = ~weight, data = clusterdt)

admin1.HT.withNA <- function(svy.design,which.area) {
  admin1 <- NULL
  tmp <- subset(svy.design, (admin1 == as.character(which.area)))
  
  if (dim(tmp)[1] == 0) {
    return(rep(NA, 5))
  }else if (sum(tmp$variables$value) == 0) {
    message(paste0(which.area, " has no death, set to NA"))
    return(rep(NA, 5))
  }else {
    glm.ob <- survey::svyglm(value/N ~ 1, design = tmp, 
                             family = stats::quasibinomial, maxit = 50, 
                             weights = tmp$variables$N)
    var.est <- stats::vcov(glm.ob)
    mean <- expit(summary(glm.ob)$coefficient[1])
    lims <- expit(logit(mean) + stats::qnorm(c((1 - CI)/2, 
                                               1 - (1 - CI)/2)) * sqrt(c(var.est)))
    return(c(mean, lims, logit(mean), var.est))
  }
}

admin1.dir <- data.frame(admin1 = sort(as.numeric(unique(clusterdt$admin1))))
admin1.dir$logit.var <- admin1.dir$logit.est <- admin1.dir$upper <- admin1.dir$lower <- admin1.dir$mean <- NA
x <- sapply(admin1.dir$admin1,function(x){admin1.HT.withNA(svy.design=my.svydesign,which.area=x)})
admin1.dir[, 2:6] <- t(x)

# get admin2 HT estimates ----

admin2.HT.withNA <- function(svy.design,which.area) {
  admin2 <- NULL
  tmp <- subset(svy.design, (admin2 == as.character(which.area)))
  
  if (dim(tmp)[1] == 0) {
    return(rep(NA, 5))
  }else if (sum(tmp$variables$value) == 0) {
    message(paste0(which.area, " has no death, set to NA"))
    return(rep(NA, 5))
  }else {
    glm.ob <- survey::svyglm(value/N ~ 1, design = tmp, 
                             family = stats::quasibinomial, maxit = 50, 
                             weights = tmp$variables$N)
    var.est <- stats::vcov(glm.ob)
    mean <- expit(summary(glm.ob)$coefficient[1])
    lims <- expit(logit(mean) + stats::qnorm(c((1 - CI)/2, 
                                               1 - (1 - CI)/2)) * sqrt(c(var.est)))
    return(c(mean, lims, logit(mean), var.est))
  }
}

admin2.dir <- data.frame(admin2 = sort(as.numeric(unique(clusterdt$admin2))))
admin2.dir$logit.var <- admin2.dir$logit.est <- admin2.dir$upper <- admin2.dir$lower <- admin2.dir$mean <- NA
x <- sapply(admin2.dir$admin2,function(x){admin2.HT.withNA(svy.design=my.svydesign,which.area=x)})
admin2.dir[, 2:6] <- t(x)

# formula to get HT for skater regions -----
skater.HT.withNA <- function(svy.design,which.area) {
  skater <- NULL
  tmp <- subset(svy.design, (skater == as.character(which.area)))
  
  if (dim(tmp)[1] == 0) {
    return(rep(NA, 5))
  }else if (sum(tmp$variables$value) == 0) {
    message(paste0(which.area, " has no death, set to NA"))
    return(rep(NA, 5))
  }else {
    glm.ob <- survey::svyglm(value/N ~ 1, design = tmp, 
                             family = stats::quasibinomial, maxit = 50, 
                             weights = tmp$variables$N)
    var.est <- stats::vcov(glm.ob)
    mean <- expit(summary(glm.ob)$coefficient[1])
    lims <- expit(logit(mean) + stats::qnorm(c((1 - CI)/2, 
                                               1 - (1 - CI)/2)) * sqrt(c(var.est)))
    return(c(mean, lims, logit(mean), var.est))
  }
}


# get admin2 BYM2 estimates with previous survey ------

overdisp.prior <- list(rho = list(param = c(0, 0.1), initial = 0))
bym2.prior <- list(phi=list(prior="pc", param=c(0.5, 2/3)),
                   prec=list(prior="pc.prec",param=c(1,0.01)))

mod1 <- inla(value ~ 1 + f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T,
                           hyper=bym2.prior),
             control.family = list(hyper = overdisp.prior),
             family='betabinomial',
             control.compute = list(config=T),
             data=clusterdt.prev, Ntrials=N)
mod1.draws <- get_inla_samples_local(mod1,1000)

eta.samples <- mod1.draws[,1:n_admin2] + mod1.draws[,ncol(mod1.draws)]
res.prev <- data.frame(admin2 = 1:n_admin2,
                       inla.median = apply(expit(eta.samples),2,median),
                       inla.sd = apply(expit(eta.samples),2,sd),
                       inla.lower = apply(expit(eta.samples),2,quantile,probs=0.025),
                       inla.upper = apply(expit(eta.samples),2,quantile,probs=0.975))

# get admin2 BYM2 estimates ------

overdisp.prior <- list(rho = list(param = c(0, 0.1), initial = 0))
bym2.prior <- list(phi=list(prior="pc", param=c(0.5, 2/3)),
                   prec=list(prior="pc.prec",param=c(1,0.01)))

mod1 <- inla(value ~ 1 + f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T,
                           hyper=bym2.prior),
             control.family = list(hyper = overdisp.prior),
             family='betabinomial',
             control.compute = list(config=T),
             data=clusterdt, Ntrials=N)
mod1.draws <- get_inla_samples_local(mod1,1000)

eta.samples <- mod1.draws[,1:n_admin2] + mod1.draws[,ncol(mod1.draws)]
admin2.res <- data.frame(admin2 = 1:n_admin2,
                       model=1,
                       inla.median = apply(expit(eta.samples),2,median),
                       inla.sd = apply(expit(eta.samples),2,sd),
                       inla.lower = apply(expit(eta.samples),2,quantile,probs=0.025),
                       inla.upper = apply(expit(eta.samples),2,quantile,probs=0.975))

adm1_wt_mat <- matrix(0,n_admin2,n_admin1)
for(area in 1:n_admin1){
  adm2_areas_tmp <- admin.key[admin.key$admin1==area,]$admin2
  adm1_wt_mat[adm2_areas_tmp,area] <- weight.adm2.u1$within_prop[adm2_areas_tmp]
}

admin1.samples <- eta.samples %*% adm1_wt_mat
admin1.res <- data.frame(admin1 = 1:n_admin1,
                         model=1,
                         inla.median = apply(expit(admin1.samples),2,median),
                         inla.sd = apply(expit(admin1.samples),2,sd),
                         inla.lower = apply(expit(admin1.samples),2,quantile,probs=0.025),
                         inla.upper = apply(expit(admin1.samples),2,quantile,probs=0.975))


# get admin1 + admin2 BYM2 estimates -----

mod2 <- inla(value ~ 1 + factor(admin1) + 
               f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T,
                 hyper=bym2.prior),
             control.family = list(hyper = overdisp.prior),
             family='betabinomial',
             control.compute = list(config=T),
             data=clusterdt, Ntrials=N)
mod2.draws <- get_inla_samples_local(mod2,1000)

eta.samples <- matrix(0,ncol=n_admin2,nrow=1000)
for(area in 1:n_admin2){
  adm1_area_tmp <- admin.key[admin.key$admin2==area,]$admin1
  eta.samples[,area] <- mod2.draws[,2*n_admin2+1] + mod2.draws[,area]
  # add admin1 fixed effect if not in admin1_1
  if(adm1_area_tmp!=1){
    eta.samples[,area] <-  eta.samples[,area] + mod2.draws[,2*n_admin2+adm1_area_tmp]
  }
}

admin2.res <- rbind(admin2.res,
                    data.frame(admin2 = 1:n_admin2,
                               model=2,
                               inla.median = apply(expit(eta.samples),2,median),
                               inla.sd = apply(expit(eta.samples),2,sd),
                               inla.lower = apply(expit(eta.samples),2,quantile,probs=0.025),
                               inla.upper = apply(expit(eta.samples),2,quantile,probs=0.975)))

admin1.samples <- eta.samples %*% adm1_wt_mat
admin1.res <- rbind(admin1.res,
                    data.frame(admin1 = 1:n_admin1,
                         model=2,
                         inla.median = apply(expit(admin1.samples),2,median),
                         inla.sd = apply(expit(admin1.samples),2,sd),
                         inla.lower = apply(expit(admin1.samples),2,quantile,probs=0.025),
                         inla.upper = apply(expit(admin1.samples),2,quantile,probs=0.975)))

# get admin1 FE estimates --------

mod3 <- inla(value ~ factor(admin1) - 1,
             control.family = list(hyper = overdisp.prior),
             family='betabinomial',
             control.compute = list(config=T),
             data=clusterdt, Ntrials=N)
mod3.draws <- get_inla_samples_local(mod3,1000)

admin1.samples <- mod3.draws
admin1.res <- rbind(admin1.res,
                    data.frame(admin1 = 1:n_admin1,
                               model=3,
                               inla.median = apply(expit(admin1.samples),2,median),
                               inla.sd = apply(expit(admin1.samples),2,sd),
                               inla.lower = apply(expit(admin1.samples),2,quantile,probs=0.025),
                               inla.upper = apply(expit(admin1.samples),2,quantile,probs=0.975)))


# make MST based on BB estimates -------

# get list of neighbors
nb <- poly2nb(poly.adm2)

# get cost of each edge
costs <- spdep::nbcosts(nb,admin2.res[admin2.res$model==1,]$inla.median)
#costs <- spdep::nbcosts(nb,res.prev$inla.median)

# get a minimum spanning tree (only unique if you start with the same node)
mst <- spdep::mstree(nb2listw(nb,costs,style = 'B'),ini = which.max(admin2.res[admin2.res$model==1,]$inla.median))
#mst <- spdep::mstree(nb2listw(nb,costs,style = 'B'),ini = which.max(res.prev$inla.median))

#plot MST (as graph)
# plot(st_geometry(poly.adm2), border=gray(.5))
# pts <- st_coordinates(st_centroid(poly.adm2))
# plot(mst, pts, col=2, cex.lab=.6, cex.circles=0.035, fg="blue", add=TRUE)

# fit models with K skater regions ------
K_seq <- seq(5,85,15)
skatermod.res <- skatermod.res.agg <- NULL
for(K in K_seq){
  print(K)
  # get groups with skater, ncuts+1 regions, each must contain at least crit areas (deterministic)
  res <- spdep::skater(edges=mst[,1:2], data=admin2.res[admin2.res$model==1,]$inla.median, ncuts=K-1, crit=3)
 # res <- spdep::skater(edges=mst[,1:2], data=res.prev$inla.median, ncuts=K-1, crit=3)
  
  poly.adm2$skater <- res$groups
  admin.key.new <- merge(admin.key,poly.adm2[,c('admin2','skater')])
  clusterdt.new <- merge(clusterdt,poly.adm2[,c('admin2','skater')])
  
  # clean_map_theme +
  #   geom_sf(data=poly.adm2,aes(fill=factor(skater))) +
  #   geom_sf(fill = "transparent", color = "grey30", lwd=0.75, data = poly.adm2 %>% group_by(skater) %>% summarise()) +
  #   theme(legend.position = 'none') + ggtitle('SKATER Regions (min size=3)')
  
  ## get direct estimates for skater region
  # my.svydesign <- survey::svydesign(ids = ~cluster, 
  #                                   strata = ~v022, nest = T, weights = ~weight, data = clusterdt.new)
  # 
  # skater.dir <- data.frame(region = sort(as.numeric(unique(clusterdt.new$skater))))
  # skater.dir$logit.var <- skater.dir$logit.est <- skater.dir$upper <- skater.dir$lower <- skater.dir$mean <- NA
  # skater.dir[, 2:6] <- t(sapply(skater.dir$region,function(x){skater.HT.withNA(svy.design=my.svydesign,which.area=x)}))
  # 
  
  poly.skater <- poly.adm2 %>% group_by(skater) %>% summarise(geometry = st_union(geometry))
  poly.skater$region <- 1:nrow(poly.skater)
  poly.skater$region.char <- paste0('region_',1:nrow(poly.skater))
  poly.skater <- sf::st_as_sf(poly.skater)
  skater.mat <- nb2mat(poly2nb(poly.skater), zero.policy = TRUE)
  colnames(skater.mat) <- rownames(skater.mat) <- poly.skater$region.char
  # check how many regions have at least 1 death
  death_regions <- (clusterdt.new %>% group_by(skater) %>% summarise(deaths=sum(value)) %>% summarise(val = sum(deaths>0)))$val
  
  # FIXED EFFECTS -- should we be doing a fixed effects model when some regions have no observed deaths? skip if not all regions have death
  if(death_regions < K){
    stop('Too many regions')
  }
  
  skatermod <- inla(value ~ 1 + factor(skater) + 
                        f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T, hyper=bym2.prior),
                      control.family = list(hyper = overdisp.prior),
                      family='betabinomial',
                      control.compute = list(config=T),
                      data=clusterdt.new, Ntrials=N)
  
    skatermod.draws <- get_inla_samples_local(skatermod,1000)
    
    eta.samples <- matrix(0,ncol=n_admin2,nrow=1000)
    for(area in 1:n_admin2){
      skater_area_tmp <- admin.key.new[admin.key.new$admin2==area,]$skater
      eta.samples[,area] <- skatermod.draws[,2*n_admin2+1] + skatermod.draws[,area]
      # add fixed effect if not in region 1
      if(skater_area_tmp!=1){
        eta.samples[,area] <-  eta.samples[,area] + skatermod.draws[,2*n_admin2+skater_area_tmp]
      }
    }
    
    skatermod.res.tmp <- data.frame(admin2 = 1:n_admin2, K = K,
                                   # terms='fixed',
                                    inla.median = apply(expit(eta.samples),2,median),
                                    inla.sd = apply(expit(eta.samples),2,sd),
                                    inla.lower = apply(expit(eta.samples),2,quantile,probs=0.025),
                                    inla.upper = apply(expit(eta.samples),2,quantile,probs=0.975))
    
    skatermod.res <- rbind(skatermod.res,skatermod.res.tmp)
    
    admin1.samples <- eta.samples %*% adm1_wt_mat
    
    skatermod.res.agg.tmp <- data.frame(admin1 = 1:n_admin1, K = K,
                                        inla.median = apply(expit(admin1.samples),2,median),
                                        inla.sd = apply(expit(admin1.samples),2,sd),
                                        inla.lower = apply(expit(admin1.samples),2,quantile,probs=0.025),
                                        inla.upper = apply(expit(admin1.samples),2,quantile,probs=0.975))
  
    skatermod.res.agg <- rbind(skatermod.res.agg,skatermod.res.agg.tmp)
  
    {
  # RANDOM IID EFFECTS
  # skatermod2 <- inla(value ~ 1 + f(skater,model='iid') + 
  #                      f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T, hyper=bym2.prior),
  #                    control.family = list(hyper = overdisp.prior),
  #                    family='betabinomial',
  #                    control.compute = list(config=T),
  #                    data=clusterdt.new, Ntrials=N)
  # skatermod2.draws <- get_inla_samples_local(skatermod2,1000)
  # 
  # eta.samples <- skatermod2.draws[,K + 2*n_admin2 +1] + #intercept
  #   skatermod2.draws[,K + 1:n_admin2] + #admin2
  #   skatermod2.draws[,admin.key.new[admin.key.new$admin2==1:n_admin2,]$skater] # skater region 
  # 
  # skatermod.res.tmp <- data.frame(admin2 = 1:n_admin2, K= K,
  #                                 terms='random_IID',
  #                                 inla.median = apply(expit(eta.samples),2,median),
  #                                 inla.sd = apply(expit(eta.samples),2,sd),
  #                                 inla.lower = apply(expit(eta.samples),2,quantile,probs=0.025),
  #                                 inla.upper = apply(expit(eta.samples),2,quantile,probs=0.975))
  # 
  # skatermod.res <- rbind(skatermod.res,skatermod.res.tmp)
  
  # RANDOM BYM2 EFFECTS
  # skatermod3 <- inla(value ~ 1 + f(skater,model='bym2',graph=skater.mat, scale.model=T, constr=T, hyper=bym2.prior) + 
  #                      f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T, hyper=bym2.prior),
  #                    control.family = list(hyper = overdisp.prior),
  #                    family='betabinomial',
  #                    control.compute = list(config=T),
  #                    data=clusterdt.new, Ntrials=N)
  # skatermod3.draws <- get_inla_samples_local(skatermod3,1000)
  # 
  # eta.samples <- skatermod3.draws[,2*K + 2*n_admin2 +1] + #intercept
  #   skatermod3.draws[,2*K + 1:n_admin2] + #admin2
  #   skatermod3.draws[,admin.key.new[admin.key.new$admin2==1:n_admin2,]$skater] # skater region 
  # 
  # skatermod.res.tmp <- data.frame(admin2 = 1:n_admin2, K= K,
  #                                 terms='random_BYM2',
  #                                 inla.median = apply(expit(eta.samples),2,median),
  #                                 inla.sd = apply(expit(eta.samples),2,sd),
  #                                 inla.lower = apply(expit(eta.samples),2,quantile,probs=0.025),
  #                                 inla.upper = apply(expit(eta.samples),2,quantile,probs=0.975))
  # 
  # skatermod.res <- rbind(skatermod.res,skatermod.res.tmp)
    }

}


# compare skater region models to admin2 and admin1 + admin2 -------

fixed_col <- 'red'
iid_col <- 'green'
bym2_col <- 'blue'


# check spread of posterior medians
skatermod.res %>% group_by(K) %>% summarise(val = sd(inla.median)) %>% ggplot() + 
  geom_point(aes(K,val),size=2) + geom_line(aes(K,val)) +
  xlab('Number of regions') + ylab('SD of Admin 2 estimates') +
  # scale_color_manual(name='Skater region effect', 
  #                    labels=c('Fixed','IID Random','BYM2 Random'),
  #                    values=c(fixed_col,iid_col,bym2_col)) +
  annotate('point',x=1,y=sd(admin2.res[admin2.res$model==1,]$inla.median),pch=4,size=5) +
  annotate('point',x=n_admin1,y=sd(admin2.res[admin2.res$model==2,]$inla.median),col=fixed_col,pch=4,size=5)

# check median of posterior SDs
skatermod.res %>% group_by(K) %>% summarise(val = median(inla.sd)) %>% ggplot() + 
  geom_point(aes(K,val),size=2) + geom_line(aes(K,val)) +
  xlab('Number of regions') + ylab('Median of SDs of Admin 2 estimates') +
  # scale_color_manual(name='Skater region effect', 
  #                    labels=c('Fixed','IID Random','BYM2 Random'),
  #                    values=c(fixed_col,iid_col,bym2_col)) +
  annotate('point',x=1,y=median(admin2.res[admin2.res$model==1,]$inla.sd),pch=4,size=5) +
  annotate('point',x=n_admin1,y=median(admin2.res[admin2.res$model==2,]$inla.sd),col=fixed_col,pch=4,size=5)


plot_theme <- clean_map_theme + 
  scale_fill_continuous(name='Posterior median',
    lim=range(skatermod.res$inla.median),
    breaks=c(0.01,0.03,0.05,0.07),
   # trans='log',
    high = "#132B43", low = "#56B1F7")

plots <- list()
plots[[1]] <- plot_theme + ggtitle('Admin2 only') +
  geom_sf(data=merge(poly.adm2,mod1.res),aes(fill=inla.median))

plots[[2]] <- plot_theme + ggtitle('Admin1 FE (37 areas)') +
  geom_sf(data=merge(poly.adm2,mod2.res),aes(fill=inla.median))

for(k in K_seq[1:4]){
  plots[[length(plots)+1]] <- plot_theme + ggtitle(paste0('K = ',k)) +
    geom_sf(data = merge(poly.adm2, skatermod.res %>% filter(terms=='fixed',K==k)),aes(fill=inla.median)) 
  
}

ggarrange(plotlist = plots,common.legend = T)

plot_theme <- clean_map_theme + 
  scale_fill_continuous(name='Posterior SD',
                        lim=range(skatermod.res$inla.sd),
                        breaks=c(0.005,0.02,0.05),
                        # trans='log',
                        high = "#132B43", low = "#56B1F7")

plots <- list()
plots[[1]] <- plot_theme + ggtitle('Admin2 only') +
  geom_sf(data=merge(poly.adm2,mod1.res),aes(fill=inla.sd))

plots[[2]] <- plot_theme + ggtitle('Admin1 FE (37 areas)') +
  geom_sf(data=merge(poly.adm2,mod2.res),aes(fill=inla.sd))

for(k in K_seq[1:4]){
  plots[[length(plots)+1]] <- plot_theme + ggtitle(paste0('K = ',k)) +
    geom_sf(data = merge(poly.adm2, skatermod.res %>% filter(terms=='fixed',K==k)),aes(fill=inla.sd)) 
  
}

ggarrange(plotlist = plots,common.legend = T)



# MSE between admin2 direct and model-based
merge(admin2.dir, skatermod.res) %>% mutate(val = (mean-inla.median)^2) %>% group_by(K) %>% summarise(sum(val,na.rm=T))
merge(admin2.dir, admin2.res) %>% group_by(model) %>% mutate(val = (mean-inla.median)^2) %>% summarise(sum(val,na.rm=T))

# Scatter plots
plot(admin2.res[admin2.res$model==2,]$inla.median,skatermod.res[skatermod.res$K==30,]$inla.median,
     xlab='Admin1 FE + Admin2',ylab='30 Skaters + Admin2')
abline(0,1)

plot(admin1.dir$mean,skatermod.res.agg[skatermod.res.agg$K==30,]$inla.median,
     xlab='Admin1 Direct',ylab='30 Skaters + Admin2')
abline(0,1)

plot(admin1.res[admin1.res$model==2,]$inla.median,skatermod.res.agg[skatermod.res.agg$K==30,]$inla.median,
     xlab='Admin1 FE + Admin2',ylab='30 Skaters + Admin2')
abline(0,1)

plot(admin1.dir$mean,admin1.res[admin1.res$model==2,]$inla.median,xlab='Admin1 Direct',ylab='Admin1 FE + Admin2')
abline(0,1)

