#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_adaptation_traits_2.h>
#include <CGAL/Regular_triangulation_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <iterator>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Weighted_point_2 Weighted_point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;
typedef CGAL::Regular_triangulation_2<K> RT2;
typedef CGAL::Regular_triangulation_adaptation_traits_2<RT2> AT;
typedef CGAL::Regular_triangulation_caching_degeneracy_removal_policy_2<RT2> DRP;
typedef CGAL::Voronoi_diagram_2<RT2, AT, DRP> VD;
typedef VD::Face_handle Face_handle;
typedef VD::Ccb_halfedge_circulator Ccb_halfedge_circulator;

// shoelace formula for polygon area
double Area_Shoelace(const std::vector<Point_2> &V) {
  if (V.size() <= 2) {
    return 0.0;
  }
  
  std::size_t n = V.size() - 1;
  double area = 0.0;
  for (std::size_t i = 1; i < n; i++) {
    area += V[i].x() * (V[i + 1].y() - V[i - 1].y());
  }
  area += V[n].x() * (V[1].y() - V[n - 1].y());
  
  return std::abs(area) / 2.0;
}

// compute the area of a Voronoi face cropped to a rectangle
double Box_Crop_Area(Face_handle f, const VD &vd, const Iso_rectangle_2 &rec) {
  // useful containers
  std::vector<Point_2> vorVert;
  CGAL::Object obj;
  Segment_2 seg;
  Ray_2 ray;
  Point_2 s, t, sB, tB;
  
  // circulators, (--) gives clockwise rotation
  Ccb_halfedge_circulator ec_start = f->ccb();
  Ccb_halfedge_circulator ec = ec_start;
  Ccb_halfedge_circulator ec_prev;
  do {
    if (!ec->is_valid()) {
      continue;
    }
    
    // recover the halfedge object from circulator
    // the result is of clockwise orientation
    const CGAL::Object heDual = vd.dual().dual(ec->dual()); 
    if (ec->is_segment()) {
      // convert the halfedge to a segment
      seg = CGAL::object_cast<Segment_2>(heDual);
      s = seg.source();
      t = seg.target();
      
      // check for intersection
      obj = CGAL::intersection(seg, rec);
      if (obj) {
        seg = CGAL::object_cast<Segment_2>(obj);
        sB = seg.source();
        tB = seg.target();
        
        // if intersection happens, then:
        // 1) only the source is in the box, save the source and the boundary point
        // 2) only the target is in the box, save the boundary point
        // 3) both source and target are in the box, save the source
        // 4) both source and target are outside of the box, save both boundaries
        if (rec.has_on_bounded_side(s) && rec.has_on_unbounded_side(t)) {
          vorVert.push_back(s);
          if (rec.has_on_boundary(sB))  {
            vorVert.push_back(sB);
          } else {
            vorVert.push_back(tB);
          }
        } else if (rec.has_on_bounded_side(t) && rec.has_on_unbounded_side(s)) {
          if (rec.has_on_boundary(sB))  {
            vorVert.push_back(sB);
          } else {
            vorVert.push_back(tB);
          }
        } else if (rec.has_on_bounded_side(s) && rec.has_on_bounded_side(t)){
          vorVert.push_back(s);
        } else {
          vorVert.push_back(sB);
          vorVert.push_back(tB);
        }
        // Rcout << s.x() << "," << s.y() << "-I-" << t.x() << "," << t.y() << std::endl;
      } else {
        continue; // no intersection, check the next halfedge
      }
    } else if (ec->is_ray()) {
      ray = CGAL::object_cast<Ray_2>(heDual);
      s = ray.source();
      
      // check for intersection
      obj = CGAL::intersection(ray, rec);
      if (obj) {
        seg = CGAL::object_cast<Segment_2>(obj);
        sB = seg.source();
        tB = seg.target();
        
        // if intersection happens, then:
        // 1) if the previous halfedge is not a ray, save the source and the boundary point
        // 2) if the previous halfedge is a ray, save the boundary point(s)
        ec_prev = ec;
        ++ec_prev;
        if (ec_prev->is_segment()) {
          vorVert.push_back(s);
        }
        
        if (rec.has_on_boundary(sB))  {
          vorVert.push_back(sB);
        } else {
          vorVert.push_back(tB);
        }
        
      } else {
        continue; // no intersection, check next halfedge
      }
    }
  } while (--ec != ec_start);
  
  // repeat the first vertex to close the circulation
  if (vorVert.size() > 2) {
    vorVert.push_back(vorVert[0]);
  }
  
  // insert rectangle vertices as needed
  
  return Area_Shoelace(vorVert);
}

// compute the area for each cell of the power diagram of Y
// [[Rcpp::export]]
NumericVector Voronoi_Area_2D(const NumericMatrix &Y) {
  // load data
  std::vector<Weighted_point_2> wpoints;
  for (int i = 0; i < Y.nrow(); i++) {
    Point_2 p(Y(i, 0), Y(i, 1));
    Weighted_point_2 wp(p, Y(i, 2));
    wpoints.push_back(wp);
  }
  RT2 rt(wpoints.begin(), wpoints.end());
  VD vd(rt);
  
  // construct a uniform box
  Iso_rectangle_2 ubox(0.0, 0.0, 1.0, 1.0);
  
  List Test_Crop;
  
  // iterate through vertices of rt, extract corresponding faces from vd
  int faceIndex;
  Ccb_halfedge_circulator ec_start;
  Ccb_halfedge_circulator ec;
  NumericVector VorArea(Y.nrow());
  std::vector<Point_2> VorVert;
  for (auto f = vd.faces_begin(); f != vd.faces_end(); ++f) {
    faceIndex = std::distance(wpoints.begin(), std::find(wpoints.begin(), wpoints.end(), f->dual()->point()));
    VorArea[faceIndex] = Box_Crop_Area(f, vd, ubox);
  }
  
  return VorArea;
}