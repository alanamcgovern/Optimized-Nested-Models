 library(spdep)
 library(tidyverse)
 library(INLA)
 library(ggpubr)
 library(parallel)
 inla.setOption(inla.timeout=200)
source("/Users/alanamcgovern/Desktop/Research/my_helpers.R")

country <- c('Nigeria','Zambia')[2]

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

# get list of neighbors
nb <- poly2nb(poly.adm2)

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

if(country=='Zambia'){
  load(file="ZMB2013NMR_clusterdt.rda")
  load(file="ZMB2018NMR_clusterdt.rda")
  
}else if(country=='Nigeria'){
  load(file="/Users/alanamcgovern/Desktop/Research/Regionalization/NGA2013NMR_clusterdt.rda")
  load(file="/Users/alanamcgovern/Desktop/Research/Regionalization/NGA2018NMR_clusterdt.rda")
  
}else{
  stop('No country info available')
}



# load admin weights COUNTRY SPECIFIC ---------

if(country=='Nigeria'){
  load(file='/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Nigeria/worldpop/adm2_weights_u1.rda')
  weight.adm2.u1 <- weight.adm2.u1[weight.adm2.u1$year==2017,]
}else if(country=='Zambia'){
  load(file='/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/Zambia/worldpop/adm2_weights_u1.rda')
  weight.adm2.u1 <- weight.adm2.u1[weight.adm2.u1$year==2017,]
}

weight.adm2.u1 <- merge(weight.adm2.u1,admin.key)
weight.adm2.u1 <- merge(weight.adm2.u1,
                        weight.adm2.u1 %>% group_by(admin1.char) %>% summarise(admin1_prop = sum(proportion))) %>%
  arrange(admin2)
weight.adm2.u1$within_prop <- weight.adm2.u1$proportion/weight.adm2.u1$admin1_prop

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
mod1.samples <- mod1.draws[,2*n_admin2+1] + mod1.draws[,1:n_admin2]

mod1.res <- data.frame(admin2 = 1:n_admin2, 
                       median = apply(expit(mod1.samples),2,median),
                       sd = apply(expit(mod1.samples),2,sd),
                       lower = apply(expit(mod1.samples),2,quantile,probs=0.05),
                       upper = apply(expit(mod1.samples),2,quantile,probs=0.95))

# get Skater model for each set of draws ------

weighted.post.samples <- skater.res <- n_models <- NULL
# how many models are averaged for each K
K_seq <- c(10,15,20,25)
K_seq <- c(10,30)

for(K in K_seq){
  print(K)
 # eta.stage2.samples <- mclapply(1:10,function(j){
  eta.stage2.samples <- list()
  mlik <- NULL
  for(j in 1:3){
    print(j)
    # get MST
    costs <- spdep::nbcosts(nb,expit(mod1.samples[j,]))
    mst <- spdep::mstree(nb2listw(nb,costs,style = 'B'),ini = which.max(expit(mod1.samples[j,])))
    
    # get skater regions and add to data
    skater_regions <- spdep::skater(edges=mst[,1:2], data=expit(mod1.samples[j,]), ncuts=K-1, crit=3)
    
    poly.adm2$skater <- skater_regions$groups
    admin.key.new <- merge(admin.key,poly.adm2[,c('admin2','skater')])
    clusterdt.new <- merge(clusterdt,poly.adm2[,c('admin2','skater')])
    
    # check that all regions have at least 1 death
    death_regions <- (clusterdt.new %>% group_by(skater) %>% summarise(deaths=sum(value)) %>% summarise(val = sum(deaths>0)))$val
    
   if(death_regions < K){
     message('Too many regions')
   }else{
      skatermod <- inla(value ~ factor(skater) -1 + 
                          f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T, hyper=bym2.prior),
                        control.family = list(hyper = overdisp.prior),
                        family='betabinomial',
                        control.compute = list(config=T),
                        data=clusterdt.new, Ntrials=N)
      
      skatermod.draws <- get_inla_samples_local(skatermod,1000)
      
      out <- skatermod.draws[,1:n_admin2] + skatermod.draws[,2*n_admin2 + admin.key.new$skater]
      eta.stage2.samples[[length(eta.stage2.samples)+1]] <- out
      mlik <- c(mlik,skatermod$mlik[1,1])
    }
    #return(out)
  }
 # },mc.cores=5)
  
  # get weights for posterior
  post_weights <- as.numeric(exp(mlik - (min(mlik) + log(sum(exp(mlik-min(mlik)))))))
  
  eta.samples <- lapply(1:length(eta.stage2.samples),function(m){data.frame(eta.stage2.samples[[m]],y=post_weights[m])})
  eta.samples.array <- do.call('rbind',eta.samples)
  
  #find quantile, q
  q <- c(0.05,0.5,0.95)
  target = q*sum(eta.samples.array$y)
  
  # for each each area, find that quantile from weighted posterior
  quants <- t(sapply(1:n_admin2,function(area){
    tmp <- eta.samples.array[,c(area,n_admin2+1)]
    tmp_ord <- tmp[order(tmp[,1]),]
    cdf <- cumsum(tmp_ord$y)
    out <- sapply(target,function(x){tmp_ord[which.min(cdf<x),][,1]})
    return(out)
  }))
  
  skater.res.tmp <- data.frame(admin2 = 1:n_admin2, K=K)
  skater.res.tmp[,3:5] <- expit(quants)
  colnames(skater.res.tmp[,3:5]) <- c('lower90','median','upper90')
  
  
  # record weighted posterior samples (in case we need other summaries later)
  weighted.post.samples[[length(weighted.post.samples)+1]] <- eta.samples.array
  # record summarized results
  skater.res <- rbind(skater.res,skater.res.tmp)
  # record how many models were averaged
  n_models[length(n_models)+1] <- length(mlik)
  
}
colnames(skater.res) <- c('admin2','K','lower90','median','upper90')

save(skater.res,file="/Users/alanamcgovern/Desktop/Research/Regionalization/ZMB2018_EmpBMA_250318.rda")

## other models to compare to -------

# admin1 direct
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


# admin1 FE + Admin2
mod2 <- inla(value ~ factor(admin1) -1 + f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T,
                                           hyper=bym2.prior),
             control.family = list(hyper = overdisp.prior),
             family='betabinomial',
             control.compute = list(config=T),
             data=clusterdt, Ntrials=N)
mod2.draws <- get_inla_samples_local(mod2,1000)

mod2.samples <- mod2.draws[,1:n_admin2] + mod2.draws[,2*n_admin2 + admin.key$admin1]
mod2.res <- data.frame(admin2 = 1:n_admin2, 
                       median = apply(expit(mod2.samples),2,median),
                       sd = apply(expit(mod2.samples),2,sd),
                       lower = apply(expit(mod2.samples),2,quantile,probs=0.05),
                       upper = apply(expit(mod2.samples),2,quantile,probs=0.95))

# aggregate models to admin1 (just medians right now) -----

mod1.res.agg <- merge(mod1.res,weight.adm2.u1) %>% group_by(admin1) %>% summarise(median = sum(median*within_prop))
mod2.res.agg <- merge(mod2.res,weight.adm2.u1) %>% group_by(admin1) %>% summarise(median = sum(median*within_prop))
skater.res.agg <- merge(skater.res,weight.adm2.u1) %>% group_by(K,admin1) %>% summarise(median = sum(median*within_prop))

# scatter plots ------

par(mfrow=c(2,2))

# compare admin2 medians -- shrinkage is reduced

## compare skater to no fixed effect
for(K in K_seq){
  plot(mod1.res$median,skater.res[skater.res$K==K,]$median,
       xlab='Admin2 model',ylab=paste0('Average of 100 Skater models with ', K, ' regions'))
  abline(0,1)
}

## compare skater to admin1 FE
for(K in K_seq){
  plot(mod2.res$median,skater.res[skater.res$K==K,]$median,
       xlab='Admin1 + Admin2 model',ylab=paste0(K, ' skater regions'))
  abline(0,1)
}


# compare admin2 SDs -- uncertainty increases
mod1.res$ci.width.scaled <- (mod1.res$upper-mod1.res$lower)/mod1.res$median
mod2.res$ci.width.scaled <- (mod2.res$upper-mod2.res$lower)/mod2.res$median
skater.res$ci.width.scaled <- (skater.res$upper90 - skater.res$lower90)/skater.res$median


## compare skater to no fixed effect
for(K in K_seq){
  plot(mod1.res$ci.width.scaled,skater.res[skater.res$K==K,]$ci.width.scaled,
       xlab='Admin2 model',ylab=paste0('Average of 10 Skater models with ', K, ' regions'),main='scaled CI width of NMR estimate')
  abline(0,1)
}

## compare skater to admin1 FE
for(K in K_seq){
  plot(mod2.res$ci.width.scaled,skater.res[skater.res$K==K,]$ci.width.scaled,
       xlab='Admin1 + Admin2 model',ylab=paste0('Average of 10 Skater models with ', K, ' regions'))
  abline(0,1)
}


# compare (aggregated) admin1 medians to admin1 direct

par(mfrow=c(3,2))
plot(admin1.dir$mean,mod1.res.agg$median,
     xlab='Admin1 Direct',ylab='Admin2 model',
     ylim=range(admin1.dir$mean))
abline(0,1)

plot(admin1.dir$mean,mod2.res.agg$median,
     xlab='Admin1 Direct',ylab='Admin1 + Admin2 model',
     ylim=range(admin1.dir$mean))
abline(0,1)

for(K in K_seq){
  plot(admin1.dir$mean,skater.res.agg[skater.res.agg$K==K,]$median,
       xlab='Admin1 Direct',ylab=paste0('', K, ' skater regions'),
       ylim=range(admin1.dir$mean))
  abline(0,1)
}

par(mfrow=c(2,2))
# compare (aggregated) admin1 medians from skater and admin1 FE models
for(K in K_seq){
  plot(mod2.res.agg$median,skater.res.agg[skater.res.agg$K==K,]$median,
       xlab='Admin1 + Admin2 model',ylab=paste0(K, ' skater regions'),
       ylim=range(mod2.res.agg$median))
  abline(0,1)
}


## map plots ----

maps <- list()
maps[[1]] <- clean_map_theme + geom_sf(data=merge(poly.adm2,mod1.res),aes(fill=median)) +
  scale_fill_viridis_c(name='NMR',lim=range(mod1.res$median,mod2.res$median,skater.res$median),
                       direction=-1,
                        breaks=c(0.03,0.04,0.05)) + 
  ggtitle('No FE')
maps[[2]] <- clean_map_theme + geom_sf(data=merge(poly.adm2,mod2.res),aes(fill=median)) +
  scale_fill_viridis_c(lim=range(mod1.res$median,mod2.res$median,skater.res$median),direction=-1) +
  theme(legend.position = 'none') +
  ggtitle('Admin1 FE (10 regions)')

for(K in K_seq){
  maps[[length(maps)+1]] <- clean_map_theme + geom_sf(data=merge(poly.adm2,skater.res[skater.res$K==K,]),aes(fill=median)) +
    scale_fill_viridis_c(lim=range(mod1.res$median,mod2.res$median,skater.res$median),direction=-1) +
    theme(legend.position = 'none') +
    ggtitle(paste0(K, ' skater regions'))
}

ggarrange(plotlist = maps,common.legend = T,nrow=2,ncol=3)


