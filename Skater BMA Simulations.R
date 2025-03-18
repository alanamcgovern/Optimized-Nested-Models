
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

