# CGALTri
## Description
An experimental R package that performs triangulation & 2D optimal transportation for uniform measure.

## Installation
In R console, run `devtools::install_github("Francis-Hsu/CGALTri")`.

For Windows user, `Makevars.win` must be modified to provide correct paths for Boost and CGAL libraries. [`vcpkg`](https://github.com/Microsoft/vcpkg) is highly recommended.

## Requirement
CGAL 4.13

## Example
```R
library(MASS)
library(CGALTri)

# example given in CGAL document
X = matrix(c(0, 0, 1, 1, 0, 1), 3, 2, byrow = T)
B = c(-1, -1, 2, 2)
TX = Triangulation_2D(X, B)
plot(TX)

# generate some random data
n = 25
X = 2 * mvrnorm(n, rep(0, 2), diag(2))

# uncropped triangulation
TX = Triangulation_2D(X)
plot(TX)

# now crop it
B = c(-2, -2, 2, 2)
TX = Triangulation_2D(X, B)
plot(TX)

# regular triangulation
w = sample(1:3 / 10, n, replace = T)
wX = cbind(X, w)
TwX = Triangulation_2D(wX, type = "Reg")
plot(TwX)
```

## Reference
CGAL, [*Computational Geometry Algorithms Library*](https://www.cgal.org).

Bruno Lévy, Erica L. Schwindt, *Notions of optimal transport theory and how to implement them on a computer*, Computers & Graphics,
Volume 72, 2018, Pages 135-148, ISSN 0097-8493, https://doi.org/10.1016/j.cag.2018.01.009.
