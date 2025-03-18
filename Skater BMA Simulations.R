
library(spdep)
library(tidyverse)
library(INLA)
inla.setOption(inla.timeout=200)
source("/Users/alanamcgovern/Desktop/Research/my_helpers.R")

# load polygons (using Zambia for now) ----
# run script to load polygon with correct names
poly.path <- paste0("/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data/shapeFiles_gadm/gadm41_ZMB_shp")
poly.adm2 <- st_read(dsn = poly.path, layer = "gadm41_ZMB_2", options = "ENCODING=UTF-8")

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

Q.admin2 <- -admin2.mat
Q.admin2 <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})*Q.admin2
diag(Q.admin2) <- sapply(1:nrow(Q.admin2),function(i){sum(I(Q.admin2[i,]!=0))})
diag(Q.admin2)[diag(Q.admin2)==0] <- 1
Q_scaled <- INLA::inla.scale.model(Q.admin2, constr=list(A=matrix(1,nrow=1,ncol=n_admin2), e=0))
Q_scaled_inv <- INLA:::inla.ginv(as.matrix(Q_scaled))

# population surface parameters (similar to Zambia 2018 frame) -----

n_clusters_urban <- pmax(100,round(rnorm(n_admin1,750,500))) #number of urban clusters in each admin1 area
n_clusters_rural <- pmax(100,round(rnorm(n_admin1,1800,800))) #number of rural clusters in each admin1 area

n_hh_urban <- 150 #average number of HH per urban cluster
sd_n_hh_urban <- 20
n_hh_rural <-  100 #average number of HH per rural cluster
sd_n_hh_rural <- 10

n_births_urban <- 2 # average number of births per HH in urban
n_births_rural <- 2.5  # average number of births per HH in rural
sd_n_births <- 0.75

# simulate population surface -----
all_dat <- data.frame(N_hh = round(c(rnorm(sum(n_clusters_urban),n_hh_urban,sd_n_hh_urban),
                                         rnorm(sum(n_clusters_rural),n_hh_rural,sd_n_hh_rural))),
                          U = c(rep(1,sum(n_clusters_urban)), rep(0,sum(n_clusters_rural))),
                           # admin area of each cluster
                          admin1 = c(unlist(sapply(1:n_admin1,function(i){
                             rep(i,n_clusters_urban[i])})),unlist(sapply(1:n_admin1,function(i){rep(i,n_clusters_rural[i])}))))
all_dat$cluster <- 1:nrow(all_dat)
# number of births in each cluster
# make sure each cluster has nonzero births
all_dat$N <- sapply(1:nrow(all_dat),function(j){ 
  if(all_dat$U[j]==1){
    round(all_dat$N_hh[j]*pmax(0.5,rnorm(1,n_births_urban,sd_n_births)))
  }else{
    round(all_dat$N_hh[j]*pmax(0.5,rnorm(1,n_births_rural,sd_n_births)))
  }
})

#assign admin2 -- CAN CHANGE TO VARY SIZES OF ADMIN2 AREAS
adm2_probs <- pmax(1,rnorm(n_admin2,10,5))
all_dat$admin2 <- sapply(1:nrow(all_dat),function(x){
  adm2.areas.tmp <- admin.key[admin.key$admin1==all_dat$admin1[x],]$admin2
#  sample(adm2.areas.tmp,1,prob = rep(1,length(adm2.areas.tmp)))
  sample(adm2.areas.tmp,1,prob = adm2_probs[adm2.areas.tmp])
})

# calculate population weights -----
admin1_UR_weights <- all_dat %>% group_by(admin1) %>% summarise(urban_prop = sum(U*N)/sum(N))
admin2_UR_weights <- all_dat %>% group_by(admin2) %>% summarise(urban_prop = sum(U*N)/sum(N))
admin2_weights <- all_dat %>% group_by(admin2,admin1) %>% summarise(prop = sum(N)/sum(all_dat$N))
admin2_to_admin1_weights <- admin2_weights %>% group_by(admin1) %>% mutate(prop = prop/sum(prop))

# risk surface parameters -----
rho <- c(0.005,0.01)[2]
sigma <- c(0.15,0.25,0.5)[3]
phi <- c(0.25,0.8)[2]
alpha <- logit(0.035)

# simulate  risk surface ------
b <- as.vector(Rfast::rmvnorm(1,rep(0,n_admin2),sigma^2*(diag((1-phi),n_admin2) + phi*Q_scaled_inv)))

# draw deaths from beta binomial
mu <- sapply(1:nrow(all_dat),function(i){expit(alpha + b[all_dat$admin2[i]])})
alpha <- mu*(1-rho)/rho
beta <- (1-mu)*(1-rho)/rho
p <- sapply(1:nrow(all_dat),function(i){rbeta(1,alpha[i],beta[i])})
all_dat$Y <- sapply(1:nrow(all_dat),function(i){rbinom(1, size = all_dat$N[i], prob = p[i])})

# sampling parameters ------
n_clusters_urban_samp <- round(0.025*n_clusters_urban) #number of urban clusters sampled in each admin1 area
n_clusters_rural_samp <- round(0.02*n_clusters_rural) #round(0.08*n_clusters_rural) #number of rural clusters sampled in each admin1 area
n_hh_samp <- 25 # number of HH sampled from each cluster 

# take sample data set -----
# sample clusters
obs_dat <- NULL
for(area in 1:n_admin1){
  # randomly select clusters from urban strata based on size of cluster
  urban_strata_dat <- all_dat[all_dat$admin1==area & all_dat$U==1,]
  rural_strata_dat <- all_dat[all_dat$admin1==area & all_dat$U==0,]
  obs_dat <- rbind(obs_dat,urban_strata_dat[sample(1:n_clusters_urban[area],n_clusters_urban_samp[area],prob = urban_strata_dat$N_hh),],
                   # randomly select cluster from rural strata  based on size of cluster
                   rural_strata_dat[sample(1:n_clusters_rural[area],n_clusters_rural_samp[area],prob = rural_strata_dat$N_hh),])
}

# sample births
obs_dat$n <- round(n_hh_samp/obs_dat$N_hh*obs_dat$N)
obs_dat$Z <- sapply(1:nrow(obs_dat),function(i){
  sum(sample(c(rep(1,obs_dat$Y[i]),rep(0,(obs_dat$N[i]-obs_dat$Y[i]))), obs_dat$n[i]))
})
obs_dat <- obs_dat[,c('admin1','admin2','U','cluster','Z','n')]

# get admin2 BYM2 estimates ------

overdisp.prior <- list(rho = list(param = c(0, 0.1), initial = 0))
bym2.prior <- list(phi=list(prior="pc", param=c(0.5, 2/3)),
                   prec=list(prior="pc.prec",param=c(1,0.01)))

mod1 <- inla(Z ~ 1 + f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T,
                           hyper=bym2.prior),
             control.family = list(hyper = overdisp.prior),
             family='betabinomial',
             control.compute = list(config=T),
             data=obs_dat, Ntrials=n)
mod1.draws <- get_inla_samples_local(mod1,1000)
mod1.samples <- mod1.draws[,2*n_admin2+1] + mod1.draws[,1:n_admin2]

mod1.res <- data.frame(admin2 = 1:n_admin2, 
                       median = apply(expit(mod1.samples),2,median),
                       sd = apply(expit(mod1.samples),2,sd),
                       lower90 = apply(expit(mod1.samples),2,quantile,probs=0.05),
                       upper90 = apply(expit(mod1.samples),2,quantile,probs=0.95))

# get admin1 + admin2 estimates --------

mod2 <- inla(Z ~ factor(admin1) -1 + f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T,
                                           hyper=bym2.prior),
             control.family = list(hyper = overdisp.prior),
             family='betabinomial',
             control.compute = list(config=T),
             data=dat.new, Ntrials=n)
mod2.draws <- get_inla_samples_local(mod2,1000)

mod2.samples <- mod2.draws[,1:n_admin2] + mod2.draws[,2*n_admin2 + admin.key$admin1]
mod2.res <- data.frame(admin2 = 1:n_admin2, 
                       median = apply(expit(mod2.samples),2,median),
                       sd = apply(expit(mod2.samples),2,sd),
                       lower90 = apply(expit(mod2.samples),2,quantile,probs=0.05),
                       upper90 = apply(expit(mod2.samples),2,quantile,probs=0.95))


# for each K, get Skater model for each set of draws and summarise weighted posterior ------

weighted.post.samples <- skater.res <- n_models <- NULL
# how many models are averaged for each K
K_seq <- c(10,15,20,25)

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
    dat.new <- merge(obs_dat,poly.adm2[,c('admin2','skater')])
    
    # check that all regions have at least 1 death
    death_regions <- (dat.new %>% group_by(skater) %>% summarise(deaths=sum(Z)) %>% summarise(val = sum(deaths>0)))$val
    
    if(death_regions < K){
      message('Too many regions')
    }else{
      skatermod <- inla(Z ~ factor(skater) -1 + 
                          f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T, hyper=bym2.prior),
                        control.family = list(hyper = overdisp.prior),
                        family='betabinomial',
                        control.compute = list(config=T),
                        data=dat.new, Ntrials=n)
      
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

# what do we want to record from each simulation? ------



