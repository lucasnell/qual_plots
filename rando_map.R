library(maps)
library(sp)
library(maptools)


da0 <- map('county', 'wisconsin,dane', fill = TRUE)
wi0 <- map('state', 'wisconsin', fill = TRUE)
nd0 <- map('state', 'north dakota', fill = TRUE)
wi_nd0 <- map('state', c('wisconsin', 'south dakota', 'north dakota', 'minnesota', 'iowa'), fill = TRUE)

ids <- sapply(strsplit(da0$names, ":"), function(x) x[1])
da <- map2SpatialPolygons(da0, IDs = ids, proj4string = CRS("+proj=longlat +datum=WGS84"))

ids <- sapply(strsplit(wi0$names, ":"), function(x) x[1])
wi <- map2SpatialPolygons(wi0, IDs = ids, proj4string = CRS("+proj=longlat +datum=WGS84"))

ids <- sapply(strsplit(nd0$names, ":"), function(x) x[1])
nd <- map2SpatialPolygons(nd0, IDs = ids, proj4string = CRS("+proj=longlat +datum=WGS84"))

ids <- sapply(strsplit(wi_nd0$names, ":"), function(x) x[1])
wi_nd <- map2SpatialPolygons(wi_nd0, IDs = ids, proj4string = CRS("+proj=longlat +datum=WGS84"))



set.seed(9)
da_aphid_pts <- spsample(da, 100, 'random')
da_wasp_pts <- spsample(da, 10, 'random')

wi_aphid_pts <- spsample(wi, 100, 'random')
wi_wasp_pts <- spsample(wi, 10, 'random')

nd_aphid_pts <- spsample(nd, 50, 'random')
nd_wasp_pts <- spsample(nd, 5, 'random')


plot(da, bg = 'white', lwd = 2)
points(da_aphid_pts, pch = 20, cex = 1.125)
points(da_wasp_pts, pch = 4, cex = 1.5)

plot(wi, bg = 'white', lwd = 1.5)
points(wi_aphid_pts, pch = 20, cex = 0.75)
points(wi_wasp_pts, pch = 4, cex = 1)


plot(wi_nd, bg = 'white', lwd = 1)
points(wi_aphid_pts[1:50,], pch = 20, cex = 0.75)
points(wi_wasp_pts[1:5,], pch = 4, cex = 1)
points(nd_aphid_pts, pch = 20, cex = 0.75)
points(nd_wasp_pts, pch = 4, cex = 1)





# # If you want a plot of Wisconsin with all its counties
# da02 <- map('county', 'wisconsin', fill = TRUE)
# ids2 <- sapply(strsplit(da02$names, ":"), function(x) x[1])
# da2 <- map2SpatialPolygons(da02, IDs = ids2, proj4string = CRS("+proj=longlat +datum=WGS84"))
#
#
# plot(da2, bg = 'transparent')
# points(da_aphid_pts, pch = 20, cex = 0.5)
# points(da_wasp_pts, pch = 4, cex = 0.75)
