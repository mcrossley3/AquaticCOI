
library(tidyverse)
library(ggmap)
library(sf)
library(raster)
library(rgdal)
library(scales)
library(mgcv)
library(mgcViz)
library(psych)
library(ggplotify)

asintrans <- function(p){# function for arcsine sqrt transformation
  asin(sqrt(p/100))
}
logtrans <- function(p){# function for log transformation
  log(p+.01)
}

zero_scale <- function(x) { # zero-center and scale
  (x-0)/sd(x)
}

center_scale <- function(x) { # mean-center and scale
  scale(x, scale = T)
}

## Function to customize 3D gam plot; Just a modification from MgcViz package
## NOTE: can use 'jet colors' to specify custom color pallette (as done below)
## OR, you can specify one of the 'cookie-cutter' color palettes from MgcViz pkg.

myvis.gam <- function (x, view = NULL, cond = list(), n.grid = 30, too.far = 0, 
                       col = NA, color = "heat", contour.col = NULL, se = -1, type = "link", 
                       plot.type = "persp", zlim = NULL, nCol = 50, ...) 
{
  fac.seq <- function(fac, n.grid) {
    fn <- length(levels(fac))
    gn <- n.grid
    if (fn > gn) 
      mf <- factor(levels(fac))[1:gn]
    else {
      ln <- floor(gn/fn)
      mf <- rep(levels(fac)[fn], gn)
      mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
      mf <- factor(mf, levels = levels(fac))
    }
    mf
  }
  dnm <- names(list(...))
  v.names <- names(x$var.summary)
  if (is.null(view)) {
    k <- 0
    view <- rep("", 2)
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) 
        ok <- FALSE
      else if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]])) <= 1) 
          ok <- FALSE
      }
      else {
        if (length(unique(x$var.summary[[i]])) == 1) 
          ok <- FALSE
      }
      if (ok) {
        k <- k + 1
        view[k] <- v.names[i]
      }
      if (k == 2) 
        break
    }
    if (k < 2) 
      stop("Model does not seem to have enough terms to do anything useful")
  }
  else {
    if (sum(view %in% v.names) != 2) 
      stop(paste(c("view variables must be one of", v.names), 
                 collapse = ", "))
    for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], 
                                 c("numeric", "factor"))) 
      stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }
  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  else {
    if (length(unique(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  if (!ok) 
    stop(paste("View variables must contain more than one value. view = c(", 
               view[1], ",", view[2], ").", sep = ""))
  if (is.factor(x$var.summary[[view[1]]])) 
    m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
  else {
    r1 <- range(x$var.summary[[view[1]]])
    m1 <- seq(r1[1], r1[2], length = n.grid)
  }
  if (is.factor(x$var.summary[[view[2]]])) 
    m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
  else {
    r2 <- range(x$var.summary[[view[2]]])
    m2 <- seq(r2[1], r2[2], length = n.grid)
  }
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  for (i in 1:length(x$var.summary)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) 
        ma <- ma[2]
    }
    if (is.matrix(x$var.summary[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), 
                          byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  if (type == "link") 
    zlab <- paste("linear predictor")
  else if (type == "response") 
    zlab <- paste("")
  else stop("type must be \"link\" or \"response\"")
  fv <- predict.gam(x, newdata = newd, se.fit = TRUE, type = type)
  z <- fv$fit
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], 
                             x$model[, view[2]], dist = too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }
  if (is.factor(m1)) {
    m1 <- as.numeric(m1)
    m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
  }
  if (is.factor(m2)) {
    m2 <- as.numeric(m2)
    m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
  }
  if (se <= 0) {
    old.warn <- options(warn = -1)
    av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, 
                 n.grid - 1)
    options(old.warn)
    max.z <- max(z, na.rm = TRUE)
    z[is.na(z)] <- max.z * 10000
    z <- matrix(z, n.grid, n.grid)
    surf.col <- t(av) %*% z %*% av
    surf.col[surf.col > max.z * 2] <- NA
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      min.z <- min(fv$fit, na.rm = TRUE)
      max.z <- max(fv$fit, na.rm = TRUE)
    }
    surf.col <- surf.col - min.z
    surf.col <- surf.col/(max.z - min.z)
    surf.col <- round(surf.col * nCol)
    con.col <- 1
    if (color == "heat") {
      pal <- heat.colors(nCol)
      con.col <- 3
    }
    else if (color == "topo") {
      pal <- topo.colors(nCol)
      con.col <- 2
    }
    else if (color == "cm") {
      pal <- cm.colors(nCol)
      con.col <- 1
    }
    else if (color == "terrain") {
      pal <- terrain.colors(nCol)
      con.col <- 2
    }
    else if (color == "gray" || color == "bw") {
      pal <- gray(seq(0.1, 0.9, length = nCol))
      con.col <- 1
    }
    ### customized here
    else if (color == 'jet') {
      pal <- jet.colors(nCol)
      con.col = 1
    }
    ####
    else stop("color scheme not recognised")
    if (is.null(contour.col)) 
      contour.col <- con.col
    surf.col[surf.col < 1] <- 1
    surf.col[surf.col > nCol] <- nCol
    if (is.na(col)) 
      col <- pal[as.array(surf.col)]
    z <- matrix(fv$fit, n.grid, n.grid)
    if (plot.type == "contour") {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",main=zlab"), ",...)", 
                    sep = "")
      if (color != "bw") {
        txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
        txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", 
                     ifelse("add" %in% dnm, "", ",add=TRUE"), ",...)", 
                     sep = "")
        eval(parse(text = txt))
      }
      else {
        txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
    else {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",zlab=zlab"), ",...)", 
                    sep = "")
      if (color == "bw") {
        op <- par(bg = "white")
        txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", 
                     stub, sep = "")
        eval(parse(text = txt))
        par(op)
      }
      else {
        txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
  }
  else {
    if (color == "bw" || color == "gray") {
      subs <- paste("grey are +/-", se, "s.e.")
      lo.col <- "gray"
      hi.col <- "gray"
    }
    else {
      subs <- paste("red/green are +/-", se, "s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      z.max <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
      z.min <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
    }
    zlim <- c(z.min, z.max)
    z <- fv$fit - fv$se.fit * se
    z <- matrix(z, n.grid, n.grid)
    if (plot.type == "contour") 
      warning("sorry no option for contouring with errors: try plot.gam")
    stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                  ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
                                                                         dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, 
                                                                                                        "", ",sub=subs"), ",...)", sep = "")
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=lo.col"), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=\"black\""), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit + se * fv$se.fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=hi.col"), stub, sep = "")
    eval(parse(text = txt))
  }
}

######################################################


dat <- read.table("EPT_COI_USA_div.txt", header = TRUE, sep = "\t")

# Cleveland dotplot; no major outliers; lots of zero's
plot(y = 1:nrow(dat),
     x = dat$Nuc_Div,
     xlab = "Nucleotide Diversity",
     ylab = "Order of the data",
     pch = 16, cex = 0.7) 

# Plot response vs. covariates
ggplot(dat, aes(x = ag_pre, y = Nuc_Div))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_grid(~Order)

ggplot(dat, aes(x = ag_post, y = Nuc_Div))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_grid(~Order)

# Pesticide:: major outliers for Ephem's and Trichops
ggplot(dat, aes(x = icide, y = Nuc_Div))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_grid(~Order)
# Pesticide again, log-scale
ggplot(dat, aes(x = log10(icide), y = Nuc_Div))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_grid(~Order)
# Nuc div vs built
ggplot(dat, aes(x = log10(built+1), y = Nuc_Div))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_grid(~Order)

ggplot(dat, aes(x = tmin, y = Nuc_Div))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_grid(~Order)

ggplot(dat, aes(x = precip, y = Nuc_Div))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_grid(~Order)

table(dat$Order)

dat$Medn_dt <- round(dat$Medn_dt, 0)
table(dat$Medn_dt) # years prior to 2008 have very few samples

# Lets see if there is geographic bias in samples prior to 2008
dat <- dat %>%
  mutate(yr = ifelse(Medn_dt < 2008, "Early","Late"))

## Enter your own key below:
register_google(key = "ZZZZZZZZZZZZZZZZZZZZZZ")

mid.lon <- as.numeric(dat %>% summarise(mean(Lon)))
mid.lat <- as.numeric(dat %>% summarise(mean(Lat)))
style1 <- c(feature = "road", element = "all", visibility = "off")

m <- get_googlemap(center = c(lon = (mid.lon-2.5), lat = mid.lat),
                   zoom = 3, scale = 2,
                   maptype = "terrain",
                   style = style1)
ggmap(m, extent = "device")+
  geom_point(aes(x = Lon, y = Lat, color = yr), data = dat)

# Look at sample variation wrt number of pairs
lons <- c(-125, -69)
lats <- c(14.7, 49)
ggmap(m, extent = "panel")+ 
  xlim(lons)+
  ylim(lats)+
  geom_point(aes(x = Lon, y = Lat, colour = Order, size = log(N_pairs+1)),
             data = dat) +
  #scale_size_continuous(range = c(0.5,3))+
  theme_bw()+
  labs(x = "longitude", y = "latitude")
## Plecoptera mainly in the east; no apparent bias in N_pairs; Ephemeroptera sparse in upper midwest

## Calculate weights as log(# sequence pairs) & relative % change in ag land
temp <- dat %>% group_by(Deme) %>% 
  arrange(Species, desc(N_seqs)) %>%
  distinct(Species, .keep_all = T) %>% 
  ungroup %>%
  droplevels %>% 
  mutate('wts' = log(N_seqs),
         'pland_dif' = (ag_post - ag_pre)) %>%
  as.data.frame()


## scale Lat/Lon
maxxy <- max(max(temp$Lat), max(temp$Lon))
temp <- temp %>% mutate(Lat2 = (Lat - min(Lat)/maxxy),
                        Lon2 = (Lon - min(Lon)/maxxy))
temp$Medn_dt <- as.factor(temp$Medn_dt)

## Look at correlations
temp <- temp %>%
  #group_by(Order) %>%
  mutate(log.icide = logtrans(icide),
         log.built = logtrans(built),
         ag_pre.trans = asintrans(ag_pre),
         ag_post.trans = asintrans(ag_post))

# Plot correlations
pairs.panels(temp[,c("ag_pre","ag_post","pland_dif","log.built","tmax","tmin","precip","log.icide")])


temp <- temp %>%
  mutate(ag_pre.scl = center_scale(ag_pre.trans),
         ag_post.scl = center_scale(ag_post.trans),
         ag_dif.scl = zero_scale(pland_dif),
         built.scl = center_scale(log.built),
         icide.scl = center_scale(log.icide),
         tmin.scl = center_scale(tmin),
         precip.scl = center_scale(precip)) %>%
  ungroup()

ggplot(temp, aes(x = ag_dif.scl, y = Nuc_Div))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_grid(~Order, scales = 'free')

png("Correlations.png", width = 7, height = 6, units = "in", res = 300)
pairs.panels(temp[,c("ag_pre.scl","ag_dif.scl","built.scl","tmin.scl","precip.scl","icide.scl")])
dev.off()

## convert columns with scaled covariates to numeric
temp.dat <- temp %>%
  mutate(across(where(is.matrix), as.numeric))


temp.dat$Order <- as.factor(temp.dat$Order)

# Model without pesticides
mod1 <- bam(Nuc_Div ~
              s(Lon, Lat, bs = 'gp', k = 25)+
              s(Mdn_dst, bs = 'tp', k = 4)+
              s(Order, bs = 're', k = 3, m = 1)+
              s(Medn_dt, bs = 're', k = 18, m = 1)+
              s(ag_pre.scl, bs = 'tp', k = 4, m = 1)+
              s(ag_dif.scl, bs = 'tp', k = 4, m = 1)+
              ti(ag_dif.scl, ag_pre.scl, bs = c('tp', 'tp'), k = c(4,4), m = 1)+
              s(tmin.scl, bs = 'tp', k = 4, m = 1)+
              s(precip.scl, bs = 'tp', k = 4, m = 1)+
              s(built.scl, bs = 'tp', k = 4, m = 1),
            data = temp.dat, family = tw, method = 'fREML',
            discrete = TRUE, weights = wts, nthreads = 2, select = TRUE)
summary(mod1)
viz.mod1 <- mgcViz::getViz(mod1)
print(plot(viz.mod1), ask = FALSE)

# Plot 3D on linear predictor scale
vis.gam(mod2, view = c("ag_dif.scl", "ag_pre.scl"), plot.type = "persp", theta = 140, too.far = 100)
# Plot 3D on response scale
vis.gam(mod1, view = c("ag_dif.scl", "ag_pre.scl"), type = 'response', plot.type = "persp", theta = 135,
        too.far = 0.5, color = "terrain")

# Nicer plot with custom colors
## Custom colors for 3D plot
jet.colors <-colorRampPalette(c("#744312","#A6611A","#CCA25C","#E6D3A5",
                                "#F5F5F5","#A7DAD2","#56B5A6","#018571","#004f43"))

as.ggplot(~myvis.gam(mod1, view = c("ag_dif.scl","ag_pre.scl"), plot.type = "persp",
                     theta = 135, type = 'response', color = 'jet', ticktype = "detailed", 
                     xlab = "", ylab = "", cex.axis = 0.8, too.far = .8))

## Run w/ insecticides
mod2 <- bam(Nuc_Div ~
              s(Lon, Lat, bs = 'gp', k = 25)+
              s(Mdn_dst, bs = 'tp', k = 4)+
              s(Order, bs = 're', k = 3, m = 1)+
              s(Medn_dt, bs = 're', k = 18, m = 1)+
              s(ag_pre.scl, bs = 'tp', k = 4, m = 1)+
              s(ag_dif.scl, bs = 'tp', k = 4, m = 1)+
              ti(ag_dif.scl, ag_pre.scl, bs = c('tp', 'tp'), k = c(4,4), m = 1)+
              s(tmin.scl, bs = 'tp', k = 4, m = 1)+
              s(precip.scl, bs = 'tp', k = 4, m = 1)+
              s(icide.scl, bs = 'tp', k = 4, m = 1)+
              s(built.scl, bs = 'tp', k = 4, m = 1),
            data = temp.dat, family = tw, method = 'fREML',
            discrete = TRUE, weights = wts, nthreads = 2, select = T)
summary(mod2)
mod1$aic; mod2$aic # mod2 has more support based on AIC

viz.mod2 <- mgcViz::getViz(mod2)
print(plot(viz.mod2), ask = FALSE)

# Plot 3D on linear predictor scale
vis.gam(mod2, view = c("ag_dif.scl", "ag_pre.scl"), plot.type = "persp", theta = 140, too.far = 0.5)
# Plot 3D on response scale
vis.gam(mod2, view = c("ag_dif.scl", "ag_pre.scl"), plot.type = "persp", theta = 145,
        too.far = 0.5, color = "terrain", type = 'response')

# Nicer plot with custom colors as defined above.
as.ggplot(~myvis.gam(mod2, view = c("ag_dif.scl","ag_pre.scl"), plot.type = "persp",
                     theta = 140, type = 'response', color = 'jet', ticktype = "detailed", 
                     xlab = "", ylab = "", cex.axis = 0.8, too.far = .5))


## Plot individual smoothers
## Ag.pre
sm1 <- plot(sm(viz.mod2, 5))
sm1 <- sm1 + l_ciPoly(level = 0.95)+
  l_fitLine(colour = "black", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.1)+
  xlab("Cropland pre-WWII")+
  scale_y_continuous(
    name = "COI diversity",
    breaks = c(-.9, -.3, 0, .3),
    # transform from predictor scale (adds intercept to log scale response)
    labels = round(exp(c(-5.1559+(-.9),-5.1559+(-.3),(-5.1559),(-5.1559+.3))),3))+
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))
png("Ag_pre.png", width = 7, height = 6, units = "in", res = 300)
sm1
dev.off()

## Ag.diff
sm2 <- plot(sm(viz.mod2, 6))
sm2 <- sm2 + l_ciPoly(level = 0.95)+
  l_fitLine(colour = "black", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.1)+
  xlab("Change in cropland")+
  scale_y_continuous(
    name = "COI diversity",
    breaks = c(-0.4, -0.2, 0, 0.2),
    # transform from predictor scale
    labels = round(exp(c(-5.1559+(-0.4), (-5.1559+(-0.2)), (-5.1559), (-5.1559+ 0.2))),3))+
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

png("Ag_diff.png", width = 7, height = 6, units = "in", res = 300)
sm2
dev.off()

## tmin
sm3 <- plot(sm(viz.mod2, 8))
sm3 <- sm3 + l_ciPoly(level = 0.95)+
  l_fitLine(colour = "black", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.1)+
  xlab("Tmin")+
  scale_y_continuous(
    name = "COI diversity",
    breaks = c(-0.4,0,0.4,0.8),
    # transform from predictor scale
    labels = round(exp(c(-5.1559+(-0.5), (-5.1559+(0)), (-5.1559+.4), (-5.1559+ 0.7))),3))+
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

png("MinTemp.png", width = 7, height = 6, units = "in", res = 300)
sm3
dev.off()

## precip
sm4 <- plot(sm(viz.mod2, 9))
sm4 <- sm4 + l_ciPoly(level = 0.95)+
  l_fitLine(colour = "black", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.1)+
  xlab("Precip")+
  scale_y_continuous(name = "COI diversity",
    breaks = c(-0.4,-.2,0, .2),
    # transform from predictor scale
    labels = round(exp(c(-5.1559+(-0.4), (-5.1559+(-.2)), (-5.1559), (-5.1559+0.2))),3))+
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

png("Precip.png", width = 7, height = 6, units = "in", res = 300)
sm4
dev.off()

## insecticide
sm5 <- plot(sm(viz.mod2, 10))
sm5 <- sm5 + l_ciPoly(level = 0.95)+
  l_fitLine(colour = "black", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.1)+
  xlab("Insecticide application rate")+
  scale_y_continuous(name = "COI diversity",
    breaks = c(-.2, 0, .2),
    # transform from predictor scale
    labels = round(exp(c(-5.1559+(-0.2), (-5.1559), (-5.1559+.2))),3))+
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

png("Insecticides.png", width = 7, height = 6, units = "in", res = 300)
sm5
dev.off()

## built
sm6 <- plot(sm(viz.mod2, 11))
sm6 <- sm6 + l_ciPoly(level = 0.95)+
  l_fitLine(colour = "black", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.1)+
  xlab("Percent impervious")+
  scale_y_continuous(name = "COI diversity",
                     breaks = c(-.2, 0, .2),
                     # transform from predictor scale
                     labels = round(exp(c(-5.1559+(-0.2), (-5.1559), (-5.1559+.2))),3))+
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

png("Built_env.png", width = 7, height = 6, units = "in", res = 300)
sm6
dev.off()


## 3D Interaction
sm.int <- as.ggplot(~myvis.gam(mod2, view = c("ag_dif.scl","ag_pre.scl"), plot.type = "persp",
                                      theta = 135, type = 'response', color = 'jet', ticktype = "detailed", 
                                      xlab = "", ylab = "", cex.axis = 0.5, too.far = .5))

png("Ag_dif_interaction.png", width = 7, height = 6, units = "in", res = 300)
sm.int
dev.off()

## Distance
smX <- plot(sm(viz.mod2, 2))
smX <- smX + l_ciPoly(level = 0.95)+
  l_fitLine(colour = "black", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.1)+
  xlab("Geographic Distance")+
  scale_y_continuous(name = "COI diversity",
                     breaks = c(-0.4, -.2, 0, .2),
                     # transform from predictor scale
                     labels = round(exp(c(-5.1559+(-0.4),-5.1559+(-0.2) , (-5.1559), (-5.1559+.2))),3))+
  theme_classic()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

png("Geographic_dist.png", width = 7, height = 6, units = "in", res = 300)
smX
dev.off()


########################################
# Map predictions 

newD <- list('Lat' = temp$Lat, 'Lon' = temp$Lon,
             'Mdn_dst' = rep(0, nrow(temp)), 'Medn_dt' = rep("2008", nrow(temp)),
             'Order' = temp$Order, 'ag_pre.scl' = rep(0, nrow(temp)),
             'ag_dif.scl' = rep(0, nrow(temp)), 'tmin.scl' = rep(0, nrow(temp)),
             'precip.scl' = rep(0, nrow(temp)), 'icide.scl' = rep(0, nrow(temp)), built.scl = rep(0,nrow(temp)))

temp.pred = temp
temp.pred$fit = predict(mod2, newD, type = 'response', se.fit = F)
max.dif <- quantile(temp.pred$fit, 0.95)
temp.pred$fit[temp.pred$fit > max.dif] <- max.dif
temp.pred$fit.sc <- scales::rescale(temp.pred$fit, to = c(0,100))

temp.pred <- temp.pred %>% mutate('Lon' = (ceiling(Lon)),
                                  'Lat' = (ceiling(Lat))) %>%
  group_by(Lon,Lat, Deme) %>%
  summarise('fit.sc' = mean(fit.sc, na.rm = T),
            'fit.or' = mean(fit), na.rm = T)

temp.pred <- temp.pred %>%
  filter(Lat < 49)

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# Crude map of US
state_map <- ne_states(iso_a2 = "US", returnclass = "sf")
state_map <- state_map %>%
  filter(!name %in% c("Alaska","Hawaii"))

temp.pred2 <- st_as_sf(temp.pred, coords = c("Lon","Lat"),
                       crs = 4326)
tiff("SmoothedCOI_V1.tif", res = 600, width = 6, height = 4, units = "in")
ggplot()+
  geom_sf(data = state_map, fill = "antiquewhite")+
  geom_point(data = temp.pred2, aes(color = fit.sc, geometry = geometry),
             stat = "sf_coordinates", shape = 19, size = 3)+
  scale_color_viridis_c(name = "Smoothed COI diversity",
                        breaks = c(0.1,97), labels = c("0",">0.009"))+
  ylab("")+
  xlab("")+
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(t=-0.6, r=0, b=0, l = 0, unit = "cm"),
        plot.margin = grid::unit(c(0,0,0,0),"mm"))+
  guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5))
dev.off()

temp.pred <- temp.pred %>%
  dplyr::select(-one_of("na.rm"))

write.csv(temp.pred, "COI_predicted_smoothV2.csv")

