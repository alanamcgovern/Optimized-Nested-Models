
library(geoR)
library(sf)
library(terra)
library(tidyverse)
library(tidyterra)
library(spdep)
library(INLA)
source("/Users/alanamcgovern/Desktop/Research/my_helpers.R")

# load population and polygon surfaces (Zambia) ----
data.dir <- "/Users/alanamcgovern/Desktop/Research/UN_Estimates/UN-Subnational-Estimates/Data"

poly.path <- paste0(data.dir,"/shapeFiles_gadm/gadm41_ZMB_shp")
poly.adm2 <- st_read(dsn = poly.path, layer = "gadm41_ZMB_2", options = "ENCODING=UTF-8")
poly.adm2$admin2 <- 1:nrow(poly.adm2)
n_admin2 <- nrow(poly.adm2)
admin2.mat <- nb2mat(poly2nb(poly.adm2), zero.policy = TRUE)
colnames(admin2.mat) <- rownames(admin2.mat) <- poly.adm2$admin2
# get list of neighbors
nb <- poly2nb(poly.adm2)


poly.adm1 <- poly.adm2 %>% group_by(NAME_1) %>% summarise(geometry = st_union(geometry))
poly.adm1$admin1 <- 1:nrow(poly.adm1)
n_admin1 <- nrow(poly.adm1)

poly.adm1.vect <- vect(poly.adm1)
poly.adm1.vect <- project(poly.adm1.vect, "EPSG:32735")
poly.adm2.vect <- vect(poly.adm2)
poly.adm2.vect <- project(poly.adm2.vect, "EPSG:32735")

admin.key <- merge(as.data.frame(poly.adm1[,c('admin1','NAME_1')]),as.data.frame(poly.adm2[,c('admin2','NAME_1')]),by='NAME_1') %>% 
  select(admin1,admin2)

# for 10 year NMR, add up 10 years of U1 surfaces
if(exists('pop_surf')){
  rm(pop_surf)
}
for(year in 2009:2018){
  tmp <- rast(paste0(data.dir,'/Zambia/Population/zmb_u1_',year,'_1km','.tif'))
  if(!exists('pop_surf')){
    pop_surf <- tmp
  }else{
    pop_surf <- pop_surf + tmp
  }
}
# convert to kilometers
pop_surf <- project(pop_surf, "EPSG:32735")

# collapse into lower res
# 6km x 6km gives the approximate right number of clusters
pop_surf <- terra::aggregate(pop_surf,fact=7,fun='sum')
pop_surf <- trim(pop_surf)

names(pop_surf) <- 'N'
values(pop_surf) <- round(values(pop_surf))
pop_surf$pixel <- 1:ncell(pop_surf)
pop_surf$pixel <- ifelse(is.na(values(pop_surf$N)),NA,values(pop_surf$pixel))
pixel_dat <- as.data.frame(pop_surf,xy=T)

# group pixels into  admin1/admin2 areas -------
wp.adm2.list <- lapply(1:nrow(poly.adm2.vect), function(x) {
  list(state_id = x, state_raster = terra::mask(crop(pop_surf,poly.adm2.vect[x,]),
                                                mask = poly.adm2.vect[x,]))
})

pixel.list <- lapply(wp.adm2.list,function(x){
  return(unique(values(x$state_raster)[,'pixel']))
})

pixel.key <- data.frame(admin2 = rep(1:115,times=unlist(lapply(pixel.list,length))),
                          pixel = unlist(pixel.list))
pixel.key <- pixel.key[!is.na(pixel.key$pixel),]

duplicated.pixels <- unique(pixel.key[duplicated(pixel.key$pixel),]$pixel)
#remove all duplicates
pixel.key.tmp <- pixel.key[!(pixel.key$pixel %in% duplicated.pixels),]

#choose area for each duplicate pixel to belong to (randomly)
pixel.alloc <- t(sapply(duplicated.pixels,function(c){
  adm1.areas <-  pixel.key[pixel.key$pixel==c,]$admin2
  return(c(sample(adm1.areas,1),c))
}))
colnames(pixel.alloc) <- c('admin2','pixel')

#combine together
pixel.key <- rbind(pixel.key.tmp,pixel.alloc)
pixel.key <- merge(pixel.key,admin.key)
pixel_dat <- merge(pixel_dat,pixel.key,by='pixel')

# assign clusters to pixels -----------------

# start with the pixels that arent too large
pixel_to_cluster <- pixel_dat[pixel_dat$N <= 500,]
pixel_to_cluster$cluster <- 1:nrow(pixel_to_cluster)

# split cells with high pop and add
high_pop_cells <- pixel_dat[pixel_dat$N>500,]$pixel

for(cell_tmp in high_pop_cells){
  tmp <- pixel_dat[pixel_dat$pixel==cell_tmp,]
  
  # how many clusters to split into?
  n_split <- ceiling(tmp$N/500)
  
  tmp2 <- tmp %>% slice(rep(1, each = n_split))
  tmp2$N <- round(tmp2$N/n_split)
  tmp2$cluster <- max(pixel_to_cluster$cluster) + 1:n_split
  pixel_to_cluster <- rbind(pixel_to_cluster,tmp2)
}

# go through each pixel which is too small and attach to adjacent cluster (preferable in same admin 1)

# find smallest cluster
tmp <- pixel_to_cluster %>% group_by(cluster) %>% summarise(val=sum(N))
smallest_cluster <- tmp[which.min(tmp$val),]$cluster

# about a minute
while(T){
  print(min(tmp$val))
  # move it to adjacent cluster
  tmp2 <- pixel_to_cluster[pixel_to_cluster$cluster==smallest_cluster,]
  # find all adjacent pixels that are not already in the same cluster
  adj_cells <- setdiff(adjacent(pop_surf,tmp2$pixel),tmp2$pixel)
  adm1_cells <- pixel_to_cluster[pixel_to_cluster$admin1 %in% unique(tmp2$admin1),]$pixel
  adj_adm1_cells <- adj_cells[adj_cells %in% adm1_cells]
  
  if(length(adj_adm1_cells)==0){
    cand_cells <- adj_cells[adj_cells %in% pixel_to_cluster$pixel]
  }else{
    cand_cells <- adj_adm1_cells
  }
  
  if(length(cand_cells) == 1){
    buddy <- cand_cells
  }else{
    buddy <- sample(cand_cells,1)
  }
  
  pixel_to_cluster[pixel_to_cluster$cluster==smallest_cluster,]$admin1 <- unique(pixel_to_cluster[pixel_to_cluster$pixel==buddy,]$admin1)
  pixel_to_cluster[pixel_to_cluster$cluster==smallest_cluster,]$admin2 <- unique(pixel_to_cluster[pixel_to_cluster$pixel==buddy,]$admin2)
  pixel_to_cluster[pixel_to_cluster$cluster==smallest_cluster,]$cluster <- unique(pixel_to_cluster[pixel_to_cluster$pixel==buddy,]$cluster)[1]
  
  # repeat until there are no small clusters
  tmp <- pixel_to_cluster %>% group_by(cluster) %>% summarise(val=sum(N))
  smallest_cluster <- tmp[which.min(tmp$val),]$cluster
  if(min(tmp$val)>100)
    break
}

pixel_to_cluster <- pixel_to_cluster[,c('pixel','cluster','admin2','admin1')]

# simulate from random GRF -----
# mem.maxVSize(vsize=50000)

coords <- crds(pop_surf)

system.time({
sim1 <- grf(grid=coords,cov.model = 'matern',
            mean=-3.5, # mean function, either scalar or vector with length equal to number of coordinates
            nugget = 0.0^2,
            kappa=1, # smoothness parameter?
            cov.pars = c(sigmasq = 0.25^2, #partial sill
                         phi = 150e3)) # range parameter (150 km)
})

# add output to population data frame
pixel_dat$r <- expit(sim1$data)

# draw deaths in each pixel (cluster) from simulated risk surface -----

rho <- 0.05
alpha <- pixel_dat$r*(1-rho)/rho
beta <- (1-pixel_dat$r)*(1-rho)/rho
p <- sapply(1:nrow(pixel_dat),function(i){rbeta(1,alpha[i],beta[i])})
pixel_dat$deaths <- sapply(1:nrow(pixel_dat),function(i){rbinom(1, size = pixel_dat$N[i], prob = p[i])})

# risk_surf <- rast(pixel_dat,crs="EPSG:32735")
# 
# plot(log(risk_surf$N))
# plot(risk_surf$r)
# hist(risk_surf$r)

# true admin1 and admin2 level rates
adm1_rates <- pixel_dat %>% group_by(admin1) %>% summarise(true_rate = sum(deaths)/sum(N))
adm2_rates <- pixel_dat %>% group_by(admin2) %>% summarise(true_rate = sum(deaths)/sum(N))

# create cluster level data from pixel level ---------

# start without high pop pixels
pop_dat <- pixel_dat[pixel_dat$N<=500,]
pop_dat <- left_join(pop_dat[,c('pixel','N','deaths')],pixel_to_cluster)

# split cells with high pop and add
high_pop_cells <- pixel_dat[pixel_dat$N>500,]
new_rows <- apply(high_pop_cells,1,function(row){
  pixel_tmp <- as.numeric(row[3])
  new_clusters <- pixel_to_cluster[pixel_to_cluster$pixel==pixel_tmp,]$cluster
  tmp <- data.frame(pixel = pixel_tmp,
                    N = round(row['N']/length(new_clusters)),
                    deaths = round(row['deaths']/length(new_clusters)),
                    cluster = new_clusters,
                    admin2 = row['admin2'],
                    admin1 = row['admin1'],
                    row.names = NULL)
  return(tmp)
})

pop_dat <- rbind(pop_dat,
                 do.call('rbind',new_rows))

pop_dat <- pop_dat %>% group_by(cluster,admin2,admin1) %>% reframe(N=sum(N),deaths=sum(deaths))

# draw some number of clusters from each admin1 ------

#choose which clusters to sample from each admin1, based on their size
sampled.clusters <- unlist(sapply(1:n_admin1,function(area){
  clusters_tmp <- pop_dat[pop_dat$admin1==area,]$cluster
  return(sample(x = clusters_tmp,
                size = round(0.04*length(clusters_tmp)),
                replace = F,
                prob = pop_dat[pop_dat$admin1==area,]$N))
}))

obs_dat <- pop_dat[pop_dat$cluster %in% sampled.clusters,]

# draw some number of births (and deaths) from each cluster -----

# drawing from 25 households, and there are about 2 deaths/HH over 10 years -- average of 50 births sampled
obs_dat$n <- round(pmin(0.75*obs_dat$N,rnorm(nrow(obs_dat),50,10)))
obs_dat$Y <- apply(obs_dat,1,function(x){
  rhyper(1,x['deaths'],x['N'] - x['deaths'], x['n'])
})

# fit basic models ---------
overdisp.prior <- list(rho = list(param = c(0, 0.1), initial = 0))
bym2.prior <- list(phi=list(prior="pc", param=c(0.5, 2/3)),
                   prec=list(prior="pc.prec",param=c(1,0.01)))

fit1 <- inla(Y ~ 1 + f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T,
                                                   hyper=bym2.prior),
                     control.family = list(hyper = overdisp.prior),
                     family='betabinomial',
                     control.compute = list(config=T),
                     data=obs_dat, Ntrials=n)
mod1.draws <- get_inla_samples_local(fit1,1000)
mod1.samples <- mod1.draws[,2*n_admin2+1] + mod1.draws[,1:n_admin2]

adm2_rates$mod1.median <- apply(expit(mod1.samples),2,median)
adm2_rates$mod1.lower <- apply(expit(mod1.samples),2,quantile,probs=0.05)
adm2_rates$mod1.upper <- apply(expit(mod1.samples),2,quantile,probs=0.95)

# compare to true population rate
adm2_rates %>% ggplot() + geom_point(aes(true_rate,mod1.median)) + geom_abline(intercept = 0, slope = 1) +
  xlim(range(adm2_rates$true_rate,adm2_rates$mod1.median)) + ylim(range(adm2_rates$true_rate,adm2_rates$mod1.median))
mean(abs(adm2_rates$true_rate - adm2_rates$mod1.median))
mean(adm2_rates$mod1.upper > adm2_rates$true_rate & adm2_rates$mod1.lower < adm2_rates$true_rate)


fit2 <- inla(Y ~ factor(admin1) - 1 + f(admin2,model='bym2',graph=admin2.mat, scale.model=T, constr=T,
                       hyper=bym2.prior),
             control.family = list(hyper = overdisp.prior),
             family='betabinomial',
             control.compute = list(config=T),
             data=obs_dat, Ntrials=n)
mod2.draws <- get_inla_samples_local(fit2,1000)
mod2.samples <- mod2.draws[,1:n_admin2] + mod2.draws[,2*n_admin2 + admin.key[order(admin.key$admin2),]$admin1]

adm2_rates$mod2.median <- apply(expit(mod2.samples),2,median)
adm2_rates$mod2.lower <- apply(expit(mod2.samples),2,quantile,probs=0.05)
adm2_rates$mod2.upper <- apply(expit(mod2.samples),2,quantile,probs=0.95)

# compare to true population rate
adm2_rates %>% ggplot() + geom_point(aes(true_rate,mod2.median)) + geom_abline(intercept = 0, slope = 1)
mean(abs(adm2_rates$true_rate - adm2_rates$mod2.median))
mean(adm2_rates$mod2.upper > adm2_rates$true_rate & adm2_rates$mod2.lower < adm2_rates$true_rate)


# for each K, get Skater model for each set of draws and summarise weighted posterior ------

weighted.post.samples <- n_models <- NULL
# how many models are averaged for each K
K_seq <- c(5,10,20)

for(K in K_seq){
  cat('K =', K)
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
    death_regions <- (dat.new %>% group_by(skater) %>% summarise(deaths=sum(Y)) %>% summarise(val = sum(deaths>0)))$val
    
    if(death_regions < K){
      message('Cannot fit skater model')
    }else{
      skatermod <- inla(Y ~ factor(skater) -1 + 
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
  
  names(mlik) <- NULL
  eta.samples <- lapply(1:length(eta.stage2.samples),function(m){data.frame(eta.stage2.samples[[m]],y=mlik[m])})
  eta.samples.array <- do.call('rbind',eta.samples)
  
  # record weighted posterior samples (in case we need other summaries later)
  weighted.post.samples[[length(weighted.post.samples)+1]] <- eta.samples.array
  # record how many models were averaged
  n_models[length(n_models)+1] <- length(mlik)
  
}

# get quantiles from weighted posterior
get_weighted_quants <- function(samples,n_areas,q=c(0.05,0.5,0.95)){
  samples$post_weights <- as.numeric(exp(samples$y - (min(samples$y) + log(sum(exp(samples$y-min(samples$y)))))))
  target = q*sum(samples$post_weights)
  
  out <- t(sapply(1:n_areas, function(area){
    tmp_ord <- samples[order(samples[,area]),]
    cdf <- cumsum(tmp_ord$post_weights)
    row <- sapply(target,function(x){tmp_ord[which.min(cdf<x),area]})
    names(row) <- q
    return(row)
  }))
  return(out)
}

# get estimates for each K
quants <- lapply(weighted.post.samples, get_weighted_quants, n_areas = n_admin2)

skater.res <- data.frame(admin2 = rep(1:n_admin2,3), K=rep(K_seq,each=n_admin2))
skater.res[,3:5] <- expit(do.call('rbind',quants))
colnames(skater.res)[3:5] <- c('lower90','median','upper90')

colnames(skater.res) <- c('admin2','K','lower90','median','upper90')
skater.res$ci.width <- skater.res$upper90 - skater.res$lower90

# compare -------------

# Compare ests to true rate
par(mfrow=c(3,2))
plot(adm2_rates$true_rate,adm2_rates$mod1.median)
abline(0,1)
plot(adm2_rates$true_rate,adm2_rates$mod2.median)
abline(0,1)
for(K in K_seq){
  plot(adm2_rates$true_rate,skater.res[skater.res$K==K,]$median)
  abline(0,1)
}

# Compare ci width (to no fixed effects model)
plot(adm2_rates$mod1.upper - adm2_rates$mod1.lower,adm2_rates$mod2.upper - adm2_rates$mod2.lower)
abline(0,1)
for(K in K_seq){
  plot(adm2_rates$mod1.upper - adm2_rates$mod1.lower,skater.res[skater.res$K==K,]$ci.width)
  abline(0,1)
}

# shrinkage factor
sd(adm2_rates$mod1.median)/sd(adm2_rates$true_rate)
sd(adm2_rates$mod2.median)/sd(adm2_rates$true_rate)
for(K in K_seq){
  print(sd(skater.res[skater.res$K==K,]$median)/sd(adm2_rates$true_rate))
}

# MAE 
mean(abs(adm2_rates$true_rate - adm2_rates$mod1.median))
mean(abs(adm2_rates$true_rate - adm2_rates$mod2.median))
for(K in K_seq){
  print(mean(abs(adm2_rates$true_rate - skater.res[skater.res$K==K,]$median)))
}

# Nominal coverage 
mean(adm2_rates$mod1.upper > adm2_rates$true_rate & adm2_rates$mod1.lower < adm2_rates$true_rate)
mean(adm2_rates$mod2.upper > adm2_rates$true_rate & adm2_rates$mod2.lower < adm2_rates$true_rate)
for(K in K_seq){
  print(mean(skater.res[skater.res$K==K,]$upper90 > adm2_rates$true_rate & skater.res[skater.res$K==K,]$lower90 < adm2_rates$true_rate))
}

save(skater.resv1,file='Sim v1.rda')

load('Sim v1.rda')

for(K in K_seq){
  plot(skater.res[skater.res$K==K,]$median,skater.resv1[skater.resv2$K==K,]$median)
}

for(K in K_seq){
  print(mean(abs(adm2_rates$true_rate - skater.resv1[skater.resv1$K==K,]$median)))
}


sqrt(1/(nrow(adm2_rates)-1)*sum((adm2_rates$mod1.median - mean(adm2_rates$mod1.median))^2))

sd(adm2_rates$mod1.median)
