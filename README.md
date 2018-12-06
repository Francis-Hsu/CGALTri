# CGALTri
## Description
An experimental R package that computes triangulation using CGAL.

## Installation
In R console, run `devtools::install_github("Francis-Hsu/CGALTri")`.

## Requirement
CGAL-4.13

## Example
```R
library(CGALTri)

# example given in the document
X = matrix(c(0, 0, 1, 1, 0, 1), 3, 2, byrow = T)
B = c(-1, -1, 2, 2)
Cropped_Voronoi_2D(X, B)
Plot_Cropped_Voronoi_2D(X, B)
Plot_Cropped_Voronoi_2D(X, B, Vor = F)

# generate some random data
X = 2 * mvrnorm(100, rep(0, 2), diag(2))

# plot an uncropped Voronoi diagram
Plot_Cropped_Voronoi_2D(X)
Plot_Cropped_Voronoi_2D(X, F)

# now crop it
B = c(-1, -1, 1, 1)
Plot_Cropped_Voronoi_2D(X, B)
Plot_Cropped_Voronoi_2D(X, B, F)
```

## Reference
CGAL, [*Computational Geometry Algorithms Library*](https://www.cgal.org)
