#### functions for drawing points in space and polygons behind them

# contents:
# - fn to infer a polygonal shape around a set of points
# - fn to arrange all the xy points from different tissues in a study in a non-overlapping way


#' Function to define a polygon boundary around a set of points in xy space
#' 
#' Fits a convex hull polygon around the points, leaving a bit of margin beyond the points
#' @param x vector of x coords
#' @param y vector of y coords
#' @param marg Amount of extra margin to draw. Default 10%
#' @return a list: x: coords of the polygon. y: y coords on the polygon
#' @example 
#'  x = rnorm(30)
#'  y = rnorm(30)
#'  plot(x, y)
#'  bound = getBoundary(x, y)
#'  polygon(bound, col = rgb(0,0,1,0.5))
getBoundary <- function(x, y, marg = 0.2) {
  # get center of shape: 
  mx = mean(x)
  my = mean(y)
  # expand all points away from center
  x2 = mx + (x - mx) * (1 + marg)
  y2 = my + (y - my) * (1 + marg)
  
  # get convex hull around expanded points:
  ch = chull(x2, y2)
  
  out = list(x = x2[ch],
             y = y2[ch])
  return(out)
}


#' Define non-overlapping x-y coords for plotting multiple tissues in the same frame
#' 
#' Given xy coords for segments from several tissues, shifts each tissue so they don't overlap
#' @param x Vector of x coords
#' @param y Vector of y coords
#' @param tissue Vector of tissue IDs corresponding to x and y
#' @param tissue.order Optional, vector of tissue names giving the order in which to plot them. 
#'  If NULL, then alphabetical order will be used
#' @param rowind Optional, integer vector of row assignments for tissues, aligned to tissue.order.
#'  Must be provided with colind to work.  
#' @param colind Optional, integer vector of column assignments for tissues, aligned to tissue.order
#'  Must be provided with rowind to work.
#' @param expansion A constant scaling factor, for how much to expand the margins between tissues
#' @return A list: x = new x coords. y = new y coords.
#' @examples 
#' # sim data from 5 tissues:
#' set.seed(0)
#' x = rnorm(50)
#' y = rnorm(50)
#' tissue = rep(letters[1:5], each = 10)
#' 
#' # arrange in default layout:
#' newxy = makeTissueGrid(x = x, y = y, tissue = tissue, 
#'                        tissue.order = NULL, rowind = NULL, colind = NULL, 
#'                        nrow = NULL, expansion = 1.2)
#' plot(newxy, pch = 16, col = 0)
#' text(newxy$x, newxy$y, tissue, col = as.numeric(as.factor(tissue)))
#' 
#' # specify a layout:
#' newxy.custom = makeTissueGrid(x = x, y = y, tissue = tissue, 
#'                               rowind = c(1,1,1,2,2), colind = c(1,2,3,1,3), 
#'                               nrow = NULL, expansion = 1.2)
#' plot(newxy.custom, pch = 16, col = 0)
#' text(newxy.custom$x, newxy.custom$y, tissue, col = as.numeric(as.factor(tissue)))
makeTissueGrid <- function(x, y, tissue, 
                           tissue.order = NULL, rowind = NULL, colind = NULL, 
                           nrow = NULL, expansion = 1.2) {
  
  # define the tissue ids and their order:
  if (length(tissue.order) > 0) {
    tissues = tissue.order
  }
  else {
    tissues = sort(unique(tissue))
  }
  
  # get tissue spans and centers:
  xspans = yspans = xcenters = ycenters = c()
  for (tiss in tissue) {
    rx = range(x[tissue == tiss])
    ry = range(y[tissue == tiss])
    xspans[tiss] = diff(rx)
    yspans[tiss] = diff(ry)
    xcenters[tiss] = median(rx)
    ycenters[tiss] = median(ry)
  }
  
  # define how far apart tissues should be in x and y space:
  xmarg = max(xspans) * expansion
  ymarg = max(yspans) * expansion
  
  # if not specified, define a grid of tissue centers:
  if ((length(rowind) == 0) | (length(colind) == 0)) {
    # define number of rows and columns:
    if (length(nrow) == 0) {
      nrow = round(sqrt(length(tissues)))
    }
    ncol = ceiling(length(tissues) / nrow)
    
    # assign rows and columns to each tissue:
    rowind = ceiling((1:length(tissues)) / ncol)
    colind = ceiling(((1:length(tissues))) %% ncol)
    colind = replace(colind, colind == 0, ncol)
    #plot(rowind ~ colind);text(colind, rowind, tissues)  
  }
  rowind = max(rowind) - rowind + 1
  names(rowind) = names(colind) = tissues
  
  # now get xy offsets for each tissue
  xnew = ynew = replace(x, TRUE, NA)
  for (tiss in tissues) {
    tempxoffset = colind[tiss] * xmarg - xcenters[tiss]
    tempyoffset = rowind[tiss] * ymarg - ycenters[tiss]
    xnew[tissue == tiss] = x[tissue == tiss] + tempxoffset
    ynew[tissue == tiss] = y[tissue == tiss] + tempyoffset
  }
  
  out = list(x = xnew, y = ynew)
  return(out)
}


# (this one needs improvements)
# dev note: what's needed here:
# - 0 and tiny values get an empty pch = 1
# - automatically infer and draw polygons within this?
# - and enable polygons to be colored by tissue ID


#' Plot a gene or score in space
#' 
#' Wrapper function for other spatial plotting functions. 
#' Draws plots the x-y coords of tissues, with point size giving variable value.
#' @param x X coords
#' @param y Y coords
#' @param z A non-negative vector of expression levels, pathway/cell scores, etc.
#' @param tissue A vector of tissue IDs
#' @param tissuecols Named vector of colors to use for each tissue's polygon
#' @param use Vector of logicals specifying which elements of the data to use
#' @param rescale Logical for whether to rescale to give z a mean of 1. 
#' @param cex A scalar, controls point size. (Points are further scaled by the values of z.)
#' @param boundaries Optional, a list of x,y vectors defining polygons, in format created by getBoundary()
#' @param ... Arguments passed to plot()
#' @return Draws a plot in xy space. Returns a list: 
#'  x: new x coords
#'  y: new y coords,
#'  outlines: list of polygonal tissue boundaries' x-y coords
#' @example 
#' # sim data from 5 tissues:
#' set.seed(0)
#' x = rnorm(50)
#' y = rnorm(50)
#' z = runif(50)
#' tissue = rep(letters[1:5], each = 10)
#' spaceplot(x, y, z, tissue, tissuecols = NULL, use = TRUE, rescale = FALSE, 
#'           cex = 1, col = "#00008B80",
#'           nrow = NULL, rowind = NULL, colind = NULL, expansion = 1.2) 
spaceplot <- function(x, y, z, tissue, tissuecols = NULL, tissuecols.alpha = 0.2, use = TRUE, rescale = FALSE, 
                      cex = 1, col = "#00008B80",
                      nrow = NULL, rowind = NULL, colind = NULL, expansion = 1.2, ...) {

  # subset everything:
  if (length(col) == length(x)) {
    col = col[use]
  }
  x = x[use]
  y = y[use]
  z = z[use]
  tissue = tissue[use]
  
  # first, get new xy coords:
  newxy.custom = makeTissueGrid(x = x, y = y, tissue = tissue, 
                                rowind = rowind, colind = colind, 
                                nrow = NULL, expansion = expansion)
  x = newxy.custom$x
  y = newxy.custom$y
  
  # also get polygon tissue boundaries:
  boundaries = list()
  for (tiss in unique(tissue)) {
    tempuse = tissue == tiss
    boundaries[[tiss]] = getBoundary(x = x[tempuse], y = y[tempuse], marg = 0.2) 
  }
  
  # transform z:
  if (rescale) {z = z / max(z)}
  
  # plotting:
  plot(x, y, 
       cex = z * cex, 
       col = col,
       pch = 16, 
       xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n", ...)
  
  # tissue boundaries:
  if (length(tissuecols) == 0) {
    tissuecols = rep("grey50", length(unique(tissue)))
    names(tissuecols) = unique(tissue)
  }
  for (tiss in names(boundaries)) {
    polygon(x = boundaries[[tiss]]$x, 
            y = boundaries[[tiss]]$y, 
            col = alpha(tissuecols[[tiss]], tissuecols.alpha), border = NA)
  }
  
  # return coordinates:
  out = list(x = x, y = y, boundaries = boundaries)
  return(out)
}




  