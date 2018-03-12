#' # Rainfall spatial reconstruction
#' In this exercise we will try to reconstruct virtual rainfall using IDW
#' algorithm suggested by Goldshtein, (2009).  
#' 
#' We will:
#' - schematize CMLs by points
#' - perform reconstruction
#' - look at results
#'   
#' You can try to play with number of points CML is represented and look at the
#' influence on the results. 

library(raster)
load('spatial_recon.RData')
source('fun_CML.r')

# define an area which will be evaluated
x.lim <- c(-743000, -723000)    
y.lim <- c(-1065000, -1054000) 

head(cml_info)


#' ### Schematize CMLs by points  
#'   
#' Use function cml.to.points() to schematize CMLs by points. Arguments:
#' - table with CML corrdinates (cartesian)
#' - Lmax  - max length of a CML segment represented by one point
#' (Lmax determines how many points will be used to schematize CML)  
#' - use funciton K.dinstance() to calculate distances between points
#' - plot CMLs and points

# represent CMLs by points
K <- cml.to.points(cml_info[, 3:6], Lmax = 2000)
k_dist <- K.distance.mtx(K)

# plot
# define plotting area
options(repr.plot.height = 3, repr.plot.width = 10)

par(mar = c(4, 4, 0.1, 0.1))
plot(0, 0, xlim = x.lim, ylim = y.lim, xlab = 'X [m]', ylab = 'Y [m]')
for (i in 1:nrow(cml_info)){
    lines(data.frame(c(cml_info[i, 3], cml_info[i, 5]),
                     c(cml_info[i, 4], cml_info[i, 6])),
          col = '#00000099', lwd = 1)
}
points(K[, 1:2], col = 2, pch = 20, cex = 1.9)
abline(h = seq(-1070000, by = 5000, length = 15), col = '#00000044', lty = 3)
abline(v = seq(-760000, by = 5000, length = 15), col = '#00000044', lty = 3)



#' ### Iterate to distribute rainfall along the CML segments  
#' - estimate CML accuracy (based on nominal characteristics using ITU
#'   power-law model and quantization info)  
#'   - distribute rainfall along CML, use funciton distribute.cmlR.1D()


# Esitmate CML accuracy
q_var <- get_cml_accuracy(Fr = cml_info$Freq, Pol = cml_info$Pol,
                          length = cml_info$length, std = 0.75)

# Estimate distributed rainfall along CMLs
K2 <- distribute.cmlR.1D(r_cml, K, k_dist, q_var, z = 5, n.iter = 10)

head(K2)

#' ### Extrapolate rainfall to regular grid  
#' - IDW extrapolation
#' - use function Rextrapolate_to_2Dgrid()
#' You can play with parameter rad - radius of influence (decorrelation distance)
#' and z (weight). Playing with other will require changes in the code bellow. 


rec <- Rextrapolate_to_2Dgrid(K2, x.lim, y.lim, g.size = 1000, q.var = q_var,
                              z = 5, rad = 7000)


#' ## Evaluate results  
#' - fill the reconstructed data into raster object
#' - plot the reconstructed and reference rainfall (in the form of raster fields)
#' - do some simple resampling to enable better comparison


# define raster object for reconstructed rainfall 
x.lim <- c(-743000, -723000)    #coordinates of left and right edge of subjected area
y.lim <- c(-1065000, -1054000)

r_rec <- raster(ncol = 20 , nrow = 11 , xmn = x.lim[1], xmx = x.lim[2],
                ymn = y.lim[1], ymx = y.lim[2])
r_rec[] <- rec$R

# Define palettes for plotting  
Rmax <- max(r_rec[], r_ref[], na.rm=T)
Rmax <- ceiling(Rmax)

expo <- 3.3
bks <- ((11:75)*Rmax/75)^expo/Rmax^(expo-1)
palette <- get_palette()

# plotting raster fields 
options(repr.plot.height = 4, repr.plot.width = 10)
par(mfcol=c(1,2))
plot(r_ref, col = palette, breaks = bks, useRaster = F, zlim = c(0, 15),
     xlim = x.lim, ylim = y.lim, legend = F)
for (i in 1:nrow(cml_info)){
    lines(data.frame(c(cml_info[i, 3], cml_info[i, 5]), c(cml_info[i, 4], cml_info[i, 6])),
          col = '#00000099', lwd = 1)
}

plot(r_rec, col = palette, breaks = bks, useRaster = F, zlim = c(0, 15), legend = F)
for (i in 1:nrow(cml_info)){
    lines(data.frame(c(cml_info[i, 3], cml_info[i, 5]), c(cml_info[i, 4], cml_info[i, 6])),
          col='#00000099', lwd = 1)
}



#' ### Do some simple map algebra and plotting  
#' - aaggregate reference rainfall to same resolution as reconstruction to
#' quantitative comparison


# aggregate reference raster
r_agr <- aggregate(r_ref, fact = 10, fun = mean) 

# plot results
options(repr.plot.height = 6, repr.plot.width = 10)
par(mfcol = c(2, 2), mar = c(2.2, 2, 4, .5))
plot(r_agr, col = palette, breaks = bks, useRaster = F, zlim = c(0, 55),
     xlim = x.lim, ylim = y.lim, legend = F, main = 'aggregated reference rainfall')
for (i in 1:nrow(cml_info)){
    lines(data.frame(c(cml_info[i, 3], cml_info[i, 5]),
                     c(cml_info[i, 4], cml_info[i, 6])),
          col = '#00000099', lwd = 1)
}

plot(r_rec, col = palette, breaks=bks, useRaster=F, zlim = c(0,55), legend = F,
     main = 'reconstructed rainfall')
for (i in 1:nrow(cml_info)){
    lines(data.frame(c(cml_info[i, 3], cml_info[i, 5]),
                     c(cml_info[i, 4], cml_info[i, 6])),
          col = '#00000099', lwd = 1)
}

plot(abs(r_rec - r_agr), col = palette, breaks = bks, useRaster = F, zlim = c(0,55), legend = F,
     main = 'absolute error')
for (i in 1:nrow(cml_info)){
    lines(data.frame(c(cml_info[i, 3], cml_info[i, 5]),
                     c(cml_info[i, 4], cml_info[i, 6])),
          col = '#00000099', lwd = 1)
}

hist(r_rec - r_agr, main = 'Error in R [mm/h]', col = 'green', n = 50)
