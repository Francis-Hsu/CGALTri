#' Plot a Voronoi Diagram with Edges Restricted to a Rectangle
#' 
#' Computes the Delaunay triangulation of a set of points and prints the Voronoi edges restricted to a given rectangle.
#' @param data Input coordinate matrix, of size \eqn{n} by \eqn{2}.
#' @param bCoord Coordinates of the box boundary.
#' @keywords Voronoi
#' @export
Plot_Cropped_Voronoi_2D = function(data, bCoord = NULL) {
  if (is.null(bCoord)) {
    bCoord = 1.1 * c(min(X[, 1]), min(X[, 2]), max(X[, 1]), max(X[, 2]))
  }
  rX = bCoord[c(1, 3)]
  rY = bCoord[c(2, 4)]
  
  # find data within the box
  dataB = subset(data, data[, 1] >= rX[1] & data[, 1] <= rX[2] & data[, 2] >= rY[1] & data[, 2] <= rY[2])
  xCoord = dataB[, 1]
  yCoord = dataB[, 2]
  
  # some parameters
  gg_color_hue = hcl(h = seq(15, 375, length = 7), l = 65, c = 100)[1:6] # pretty color
  par(mar = rep(0.1, 4)) # remove margin
  
  # plot the Voronoi diagram
  S = Cropped_Voronoi_2D(X, bCoord) # compute Voronoi segements
  plot(xCoord, yCoord, xlim = rX, ylim = rY, pch = 20,
       xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
  rect(bCoord[1], bCoord[2], bCoord[3], bCoord[4], border = gg_color_hue[5], lwd = 2)
  for (i in 1:nrow(S)) {
    segments(S[i, 1], S[i, 2], S[i, 3], S[i, 4], col = gg_color_hue[1], lwd = 1)
  }
}
