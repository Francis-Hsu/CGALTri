#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <iterator>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Weighted_point_2 Weighted_point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation_2;
typedef CGAL::Regular_triangulation_2<K> Regular_triangulation;

// A class to recover dual diagram from stream.
// Rays, lines and segments are cropped to a rectangle
// so that only segments are stored
// adapted from print_cropped_voronoi.cpp
struct Crop_Dual {
  std::list<Segment_2> dual_seg;
  Iso_rectangle_2 Rec;
  
  Crop_Dual(const Iso_rectangle_2 &rec):Rec(rec) {}
  
  template <class RSL>
  void Crop_to_Segment(const RSL &rsl) {
    CGAL::Object obj = CGAL::intersection(rsl, Rec);
    const Segment_2 *s = CGAL::object_cast<Segment_2>(&obj);
    if (s) {
      dual_seg.push_back(*s);
    }
  }
  
  void operator<<(const Ray_2 &ray) { 
    Crop_to_Segment(ray);
  }
  void operator<<(const Line_2 &line) { 
    Crop_to_Segment(line);
  }
  void operator<<(const Segment_2 &seg) { 
    Crop_to_Segment(seg);
  }
};

// [[Rcpp::export]]
List Delaunay_Tri_2D(const NumericMatrix &data, const NumericVector &bCoord) {
  // load data
  std::vector<Point_2> points;
  for (int i = 0; i < data.nrow(); i++) {
    Point_2 p(data(i, 0), data(i, 1));
    points.push_back(p);
  }
  
  // insert points into the triangulation
  Delaunay_triangulation_2 dt2;
  dt2.insert(points.begin(), points.end());
  
  // construct a rectangle
  Iso_rectangle_2 bbox(bCoord(0), bCoord(1), bCoord(2), bCoord(3));
  Crop_Dual vor(bbox);
  
  // extract the cropped Voronoi diagram
  dt2.draw_dual(vor);
  
  // iteration indices
  int itdx = 0;
  
  // matrix storing edges of the triangulation
  Point_2 v;
  NumericMatrix tEdges(3 * dt2.number_of_faces(), 4);
  for(auto it = dt2.finite_faces_begin(); it != dt2.finite_faces_end(); ++it) {
    for(int j = itdx; j <= itdx + 2; j++) {
      v = it->vertex(j - itdx)->point();
      tEdges(j, 0) = v.x();
      tEdges(j, 1) = v.y();
      tEdges(itdx + (j - itdx + 2) % 3, 2) = v.x();
      tEdges(itdx + (j - itdx + 2) % 3, 3) = v.y();
    }
    itdx = itdx + 3;
  }
  
  // matrix storing segments of the Voronoi diagram
  itdx = 0;
  Point_2 s, t;
  NumericMatrix vEdges(vor.dual_seg.size(), 4);
  for(auto it = vor.dual_seg.begin(); it != vor.dual_seg.end(); ++it) {
    s = it->source();
    t = it->target();
    vEdges(itdx, 0) = s.x();
    vEdges(itdx, 1) = s.y();
    vEdges(itdx, 2) = t.x();
    vEdges(itdx, 3) = t.y();
    itdx++;
  }
  
  List ret;
  ret["Tri"] = tEdges;
  ret["Vor"] = vEdges;
    
  return ret;
}

// [[Rcpp::export]]
List Regular_Tri_2D(const NumericMatrix &data, const NumericVector &bCoord) {
  // load data
  std::vector<Weighted_point_2> wpoints;
  for (int i = 0; i < data.nrow(); i++) {
    Point_2 p(data(i, 0), data(i, 1));
    Weighted_point_2 wp(p, data(i, 2));
    wpoints.push_back(wp);
  }
  Regular_triangulation rt(wpoints.begin(), wpoints.end());
  
  // construct a rectangle
  Iso_rectangle_2 bbox(bCoord(0), bCoord(1), bCoord(2), bCoord(3));
  Crop_Dual pow(bbox);
  
  // extract the cropped power diagram
  rt.draw_dual(pow);
  
  // iteration indices
  int itdx = 0;
  
  // matrix storing edges of the triangulation
  Weighted_point_2 v;
  NumericMatrix tEdges(3 * rt.number_of_faces(), 4);
  for(auto it = rt.finite_faces_begin(); it != rt.finite_faces_end(); ++it) {
    for(int j = itdx; j <= itdx + 2; j++) {
      v = it->vertex(j - itdx)->point();
      tEdges(j, 0) = v.x();
      tEdges(j, 1) = v.y();
      tEdges(itdx + (j - itdx + 2) % 3, 2) = v.x();
      tEdges(itdx + (j - itdx + 2) % 3, 3) = v.y();
    }
    itdx = itdx + 3;
  }
  
  // matrix storing segments of the power diagram
  itdx = 0;
  Point_2 s, t;
  NumericMatrix pEdges(pow.dual_seg.size(), 4);
  for(auto it = pow.dual_seg.begin(); it != pow.dual_seg.end(); ++it) {
    s = it->source();
    t = it->target();
    pEdges(itdx, 0) = s.x();
    pEdges(itdx, 1) = s.y();
    pEdges(itdx, 2) = t.x();
    pEdges(itdx, 3) = t.y();
    itdx++;
  }
  
  List ret;
  ret["Tri"] = tEdges;
  ret["Lag"] = pEdges;
  
  return ret;
}