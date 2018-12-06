#include <Rcpp.h>
using namespace Rcpp;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <iterator>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;
typedef CGAL::Delaunay_triangulation_2<K>  Delaunay_triangulation_2;

// A class to recover Voronoi diagram from stream.
// Rays, lines and segments are cropped to a rectangle
// so that only segments are stored
struct Cropped_voronoi_from_delaunay{
  std::list<Segment_2> m_cropped_vd;
  Iso_rectangle_2 m_bbox;
  
  Cropped_voronoi_from_delaunay(const Iso_rectangle_2& bbox):m_bbox(bbox){}
  
  template <class RSL>
  void crop_and_extract_segment(const RSL& rsl) {
    CGAL::Object obj = CGAL::intersection(rsl, m_bbox);
    const Segment_2* s = CGAL::object_cast<Segment_2>(&obj);
    if (s) {
      m_cropped_vd.push_back(*s);
    }
  }
  
  void operator<<(const Ray_2& ray)    { crop_and_extract_segment(ray); }
  void operator<<(const Line_2& line)  { crop_and_extract_segment(line); }
  void operator<<(const Segment_2& seg){ crop_and_extract_segment(seg); }
};

// [[Rcpp::export]]
List Cropped_Voronoi_2D(const NumericMatrix &data, const NumericVector &bCoord){
  // load some points
  std::vector<Point_2> points;
  for (int i = 0; i < data.nrow(); i++) {
    double x = data(i, 0);
    double y = data(i, 1);
    Point_2 p(x, y);
    points.push_back(p);
  }
  
  // insert points into the triangulation
  Delaunay_triangulation_2 dt2;
  dt2.insert(points.begin(), points.end());
  
  // construct a rectangle
  Iso_rectangle_2 bbox(bCoord(0), bCoord(1), bCoord(2), bCoord(3));
  Cropped_voronoi_from_delaunay vor(bbox);
  
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
  itdx = 0;
  
  // matrix storing segments of the Voronoi diagram
  Point_2 s, t;
  NumericMatrix vEdges(vor.m_cropped_vd.size(), 4);
  for(auto it = vor.m_cropped_vd.begin(); it != vor.m_cropped_vd.end(); ++it) {
    s = it->source();
    t = it->target();
    vEdges(itdx, 0) = s.x();
    vEdges(itdx, 1) = s.y();
    vEdges(itdx, 2) = t.x();
    vEdges(itdx, 3) = t.y();
    itdx++;
  }
  
  List ret;
  ret["tEdges"] = tEdges;
  ret["vEdges"] = vEdges;
    
  return ret;
}

