library(ggplot2)
library(data.table)
library(scales)
library(patchwork)
library(viridis)

allele_area <- function(r, # radius of circle
                        landscape_width = 25
) {
  a=landscape_width/2 # maximum radius before allele exits circle
  if (r <= a) { # circle is completely within the landscape
    area=pi*(r**2)
  }
  else if (r >= (sqrt(2*(a**2)))) { # circle is completely outside the landscape
    area=(a*2)**2
  }
  else { # circle exits the landscape - do some fancy trig
    # https://math.stackexchange.com/questions/1450961/overlapping-area-between-a-circle-and-a-square
    theta=acos(a/r)
    area=(((pi-(4*theta))*(r**2))+(4*a*r*sin(theta)))
  }
  return(area)
}

allele_cone <- function(het_radius, # radius of heterozygotes
                        hom_radius, # radius of homozygotes
                        #hom_area, # area of homozygotes
                        landscape_width = 25){
  if (het_radius == 0){ return(0) }
  rmax = landscape_width/2 # radius the cone starts exiting the circle
  Rmax=sqrt(2*(rmax**2)) # radius @ which circle encompasses landscape
  slope=-1/(het_radius-hom_radius) # slope of hypotenuse
  intercept=-slope * het_radius # intercept of hypotenuse
  fx <- function(x){
    return(x*slope + intercept)
  }
  fy <- function(y){
    return((y - intercept)/slope)
  }
  if (het_radius <= rmax){
    # volume of het cone
    allelevolume <- pi * het_radius^2 * (intercept / 3)
    # cut off top of cone
    alleletop <- pi * hom_radius^2 * ((intercept-1)/3)
    allelevolume <- allelevolume-alleletop
  }
  if ((het_radius > rmax) & (het_radius <= Rmax)) {
    allelevolume <- pi * het_radius^2 * (intercept / 3)
    # cut off top of cone
    alleletop <- pi * hom_radius^2 * ((intercept-1)/3)
    allelevolume <- allelevolume-alleletop
    h = intercept # height of whole cone
    r = het_radius # radius of whole cone
    d = rmax # distance from origin to slice
    # cone slices https://math.stackexchange.com/questions/1660071/how-to-calculate-volume-of-a-right-circular-cones-hyperbola-segment
    slicevol <- (h/3) * 
      ((-2 * d * sqrt(r^2 - d^2)) +
         ((d^3 / r) * log(sqrt(r^2 - d^2) + r)) +
         ((r^2 * acos(d/r)) - ((d^3 / r) * log(d))))
    h1 <- intercept - 1
    r1 <- fy(1)
    if (h1 > 0){
      #cone volume above y=1 was already subtracted, so remove that from the slice
      outside <- h1/3 * 
        ((-2 * d * sqrt(r1^2 - d^2)) +
           ((d^3 / r1) * log(sqrt(r1^2 - d^2) + r1)) +
           (r1^2 * acos(d/r1) - ((d^3 / r1) * log(d))))
      if (is.na(outside)){slicevol = slicevol}
      else {slicevol <- slicevol - outside}
    }
    allelevolume <- allelevolume - 4*slicevol
  }
  if (het_radius > Rmax){
    # full volume of 25x25 box beneath cone
    circleintersectY <- fx(Rmax)
    boxvol <- 25 * 25 * circleintersectY
    # cone on top
    coneheight <- intercept - circleintersectY
    conevol <- pi * Rmax^2 * (coneheight/3)
    conetop <- pi * hom_radius^2 * ((intercept-1)/3)
    conevol <- conevol - conetop
    # slices from the outside
    d <- rmax
    r <- Rmax
    h <- coneheight
    slicevol <- (h/3) * 
      ((-2 * d * sqrt(r^2 - d^2)) +
         ((d^3 / r) * log(sqrt(r^2 - d^2) + r)) +
         ((r^2 * acos(d/r)) - ((d^3 / r) * log(d))))
    h1 <- intercept - 1
    r1 <- fy(1)
    if (h1 > 0){
      # cone volume above y=1 was already subtracted, so remove that from the slice
      outside <- h1/3 * 
        ((-2 * d * sqrt(r1^2 - d^2)) +
           ((d^3 / r1) * log(sqrt(r1^2 - d^2) + r1)) +
           (r1^2 * acos(d/r1) - ((d^3 / r1) * log(d))))
      if (is.na(outside)){slicevol = slicevol}
      else {slicevol <- slicevol - outside}
    }
    allelevolume <- (conevol - 4*slicevol) + boxvol
  }
  return((allelevolume) / (25*25))
}

landscape_width=25
Fisher_wave_H <- function(M, # migration rate
                          S, # selection coefficient
                          landscape_width = 25) {
  landscape_area=landscape_width**2
  rmax=landscape_width/2 # radius @ which circle starts escaping the landscape
  Rmax=sqrt(2*(rmax**2)) # radius @ which circle encompasses landscape
  
  lambd=sqrt(M/S)
  V=M*sqrt(S) # allele velocity
  W=1/(4*sqrt(S)) # width of wavefront?????????????????????????????????????????????????
  hom_radii=seq(0, Rmax, V) # radius of the homozygous allele area circle after each generation
  het_radii=seq(0, Rmax+W, V) # radius of the heterozygous allele area circle after each generation
  hom_radii <- c(rep(0, length(het_radii) - length(hom_radii)), hom_radii) # add het only generations
  het_areas <- sapply(het_radii, allele_area)
  allele_areas <- het_areas
  het_freq <- mapply(allele_cone, het_radius=het_radii, hom_radius=hom_radii)
  frequencies <- het_freq
  return(data.table(allele_areas, frequencies))
}
for (s in c(0.0001, 0.001, 0.01, 0.1)){
  df <- Fisher_wave_H(0.1, s)
  df$s <- s
  df$tick <- seq(nrow(df))
  assign(paste0('S_',s), df)
}

df <- rbind(`S_1e-04`, S_0.001, S_0.01, S_0.1)
fq <- ggplot(df) + theme_bw() + 
      geom_line(aes(x=tick, y=frequencies, group=s, color=s), linewidth=0.75) + 
      scale_color_viridis(trans='log10', option='G', end = 0.9, name='Selection\ncoefficient') +
      scale_x_continuous(trans='log10') +
      xlab('Tick') +
      ylab('Frequency')
jdist <- ggplot(df) + theme_bw() + 
  geom_line(aes(x=frequencies, y=allele_areas, group=s, color=s), linewidth=0.75) + 
  scale_color_viridis(trans='log10', option='G', end = 0.9, name='Selection\ncoefficient') +
  xlab('Frequency') +
  ylab('Area')

fq + jdist + plot_annotation(tag_levels = c('A')) + plot_layout(guides = "collect")
ggsave('figures/supplement-fisherwave.pdf', width=6, height=3, units='in')
