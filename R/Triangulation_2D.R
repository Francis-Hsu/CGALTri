#' Plot a Voronoi Diagram with Edges Restricted to a Rectangle
#' 
#' Computes the Delaunay triangulation of a set of points and prints the Voronoi edges restricted to a given rectangle.
#' @param data Input coordinate matrix, of size \eqn{n} by \eqn{2}.
#' @param bCoord Coordinates of the box boundary. In the order of: xmin, ymin, xmax, ymax.
#' @param type Type of triangulation to perform, choose between Delaunay ("Del") and regular ("Reg").
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#' @keywords Voronoi
#' @export
Triangulation_2D = function(data, bCoord = NULL, type = "Del", na.rm = T) {
  # validate type
  type.valid = c("Del", "Reg")
  type = match.arg(type, type.valid)
  
  # validate data
  if (!is.matrix(data) || ncol(data) < 2) {
    stop("Data must be a matrix with ncol >= 2.")
  }
  
  # remove NAs
  if (na.rm) {
    data = data[complete.cases(data), ]
  }
  
  # default box
  if (is.null(bCoord)) {
    bCoord = 1.1 * c(min(data[, 1]), min(data[, 2]), 
                     max(data[, 1]), max(data[, 2]))
  }
  rX = bCoord[c(1, 3)]
  rY = bCoord[c(2, 4)]
  
  # find data within the box
  dataB = subset(data, data[, 1] >= rX[1] & data[, 1] <= rX[2] & data[, 2] >= rY[1] & data[, 2] <= rY[2])
  
  # choose what to plot
  if (type == "Del") {
    S = Delaunay_Tri_2D(dataB[, 1:2], bCoord)
    ret = list(type = "Delaunay", Tri = S$Tri, Vor = S$Vor, Data = dataB, Box = bCoord)
  } else if (type == "Reg") {
    S = Regular_Tri_2D(dataB[, 1:3], bCoord)
    ret = list(type = "Regular", Tri = S$Tri, Lag = S$Lag, Data = dataB, Box = bCoord)
  }
  class(ret) = "CGALTri_2D"
  
  return(ret)
}

plot.CGALTri_2D = function(Obj, type = "Both") {
  # validate type
  type.valid = c("Tri", "Vor", "Lag", "Both")
  type = match.arg(type, type.valid)
  
  # some plotting parameters
  gg_color_hue = hcl(h = seq(15, 375, length = 7), l = 65, c = 100)[1:6] # pretty color
  par(mar = rep(0.1, 4)) # remove margin
  
  # triangulation type
  tri.id = switch(Obj$type, Delaunay = 1, Regular = 2)
  
  # extract data
  box = Obj$Box
  x = Obj$Data[, 1]
  y = Obj$Data[, 2]
  if (tri.id == 2) {
    w = Obj$Data[, 3]
  }
  rX = box[c(1, 3)]
  rY = box[c(2, 4)]
  
  # plot triangulation
  if (type == "Tri" || type == "Both") {
    S = Obj$Tri
    # plot the data
    plot(x, y, xlim = rX, ylim = rY, pch = 20,
         xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
    rect(box[1], box[2], box[3], box[4], border = gg_color_hue[5], lwd = 2)
    for (i in 1:nrow(S)) {
      segments(S[i, 1], S[i, 2], S[i, 3], S[i, 4], col = gg_color_hue[1], lwd = 1)
    }
  }
  
  # plot Laguerre-Voronoi diagram
  if (type == "Vor" || type == "Pow" || type == "Both") {
    if (tri.id == 1) {
      S = Obj$Vor
    } else if (tri.id == 2) {
      S = Obj$Lag
    }
    plot(x, y, xlim = rX, ylim = rY, pch = 20,
         xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
    rect(box[1], box[2], box[3], box[4], border = gg_color_hue[5], lwd = 2)
    for (i in 1:nrow(S)) {
      segments(S[i, 1], S[i, 2], S[i, 3], S[i, 4], col = gg_color_hue[1], lwd = 1)
    }
    if (tri.id == 2) {
      circ.col = rgb(t(col2rgb(gg_color_hue[3])), alpha = 100, maxColorValue = 255)
      symbols(x, y, sqrt(w), fg = NULL, bg = circ.col, inches = F, add = T)
    }
  }
}
