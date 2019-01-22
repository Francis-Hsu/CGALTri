#include <Rcpp.h>
using namespace Rcpp;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_utils_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT Weight;
typedef K::Point_3 Point_3;
typedef K::Weighted_point_3 Weighted_point_3;
typedef CGAL::Regular_triangulation_3<K> Regular_triangulation;

// [[Rcpp::export]]
NumericMatrix Regular_Tri_3D(const NumericMatrix &data) {
  // load data
  std::vector<Weighted_point_3> wpoints;
  for (int i = 0; i < data.nrow(); i++) {
    Point_3 p(data(i, 0), data(i, 1), data(i, 2));
    Weighted_point_3 wp(p, data(i, 3));
    wpoints.push_back(wp);
  }
  Regular_triangulation rt(wpoints.begin(), wpoints.end());
  
  // iteration indices
  int itdx = 0;
  
  // matrix storing edges of the triangulation
  Weighted_point_3 s, t;
  NumericMatrix tEdges(rt.number_of_finite_edges(), 6);
  for(auto it = rt.finite_edges_begin(); it != rt.finite_edges_end(); ++it) {
    s = it->first->vertex(it->second)->point();
    t = it->first->vertex(it->third)->point();
    tEdges(itdx, 0) = s.x();
    tEdges(itdx, 1) = s.y();
    tEdges(itdx, 2) = s.z();
    tEdges(itdx, 3) = t.x();
    tEdges(itdx, 4) = t.y();
    tEdges(itdx, 5) = t.z();
    itdx++;
  }
  
  return tEdges;
}