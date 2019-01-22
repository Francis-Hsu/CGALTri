#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_adaptation_traits_2.h>
#include <CGAL/Regular_triangulation_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/number_utils.h>
#include <iterator>
#include <utility>
#include <tuple>
#include <algorithm>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
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

// which side is point p on? 
// index starts from the bottom side in counterclockwise order.
int which_side(const Point_2 &p, const Iso_rectangle_2 &b, bool is_exit) {
  int side;
  if (!b.has_on_boundary(p)) {
    side = -1; // p is not on the boundary
  } else if (p.y() == b.ymin() && p != b[0]) {
    side = 0; // bottom
  } else if (p.x() == b.xmax() && p != b[1]) {
    side = 1; // right
  } else if (p.y() == b.ymax() && p != b[2]) {
    side = 2; // top
  } else {
    side = 3; // left
  }
  
  if (p == b[0] || p == b[1] || p == b[2] || p == b[3]) {
    side += !is_exit;
    side %= 4;
  }
  
  return side;
}

// shoelace area formula for polygon
double Area_Shoelace(const std::vector<Point_2> &V) {
  if (V.size() <= 2) {
    return 0.0;
  }
  
  std::size_t n = V.size() - 1;
  CGAL::Lazy_exact_nt<CGAL::Gmpq> area = 0.0;
  for (std::size_t i = 1; i < n; i++) {
    area += V[i].x() * (V[i + 1].y() - V[i - 1].y());
  }
  area += V[n].x() * (V[1].y() - V[n - 1].y());
  
  return std::abs(CGAL::to_double(area)) / 2.0;
}

// compute the area of a Voronoi face cropped to a rectangle
std::vector<Point_2> Face_Rect_Crop(Face_handle f, const VD &vd, const Iso_rectangle_2 &rec) {
  // useful containers
  std::vector<Point_2> cropped_face_vertices;
  CGAL::Object obj;
  Segment_2 seg;
  Ray_2 ray;
  bool is_outgoing;
  Point_2 s, t, sB, tB;
  
  // containers for helping inserting vertices from the rectangle
  std::vector<std::tuple<size_t, int, bool>> intersect_points; // format: (index, side = 0-3, is_exit)
  std::vector<std::pair<size_t, std::vector<Point_2>>> vert_insertion;
  
  // circulators, (--) gives clockwise rotation
  Ccb_halfedge_circulator ec_start = f->ccb();
  Ccb_halfedge_circulator ec = ec_start;
  
  // extract vertices of the cropped face in clockwise order
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
        // if intersection happens, then:
        seg = CGAL::object_cast<Segment_2>(obj);
        sB = seg.source();
        tB = seg.target();
        
        if (rec.has_on_bounded_side(s) && rec.has_on_unbounded_side(t)) {
          // 1. only the source is in the box, save the source and the boundary point
          cropped_face_vertices.push_back(s);
          if (rec.has_on_boundary(sB))  {
            cropped_face_vertices.push_back(sB);
            intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(sB, rec, true), true));
          } else {
            cropped_face_vertices.push_back(tB);
            intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(tB, rec, true), true));
          }
        } else if (rec.has_on_bounded_side(t) && rec.has_on_unbounded_side(s)) {
          // 2. only the target is in the box, save the boundary point
          if (rec.has_on_boundary(sB))  {
            cropped_face_vertices.push_back(sB);
            intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(sB, rec, false), false));
          } else {
            cropped_face_vertices.push_back(tB);
            intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(tB, rec, false), false));
          }
        } else if (rec.has_on_bounded_side(s) && rec.has_on_bounded_side(t)){
          // 3. both source and target are in the box, save the source
          cropped_face_vertices.push_back(s);
        } else {
          // 4. both source and target are outside of the box, save both boundaries
          cropped_face_vertices.push_back(sB);
          intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(sB, rec, false), false));
          
          cropped_face_vertices.push_back(tB);
          intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(tB, rec, true), true));
        }
      } else {
        continue; // no intersection, check the next halfedge
      }
    } else if (ec->is_ray()) {
      is_outgoing = ec->has_target(); // a halfedge ray is outgoing if it has target
      ray = CGAL::object_cast<Ray_2>(heDual);
      s = ray.source();
      
      // check for intersection
      obj = CGAL::intersection(ray, rec);
      if (obj) {
        // if intersection happens, then:
        seg = CGAL::object_cast<Segment_2>(obj);
        sB = seg.source();
        tB = seg.target();
        
        if (rec.has_on_bounded_side(s)) {
          // 1. the source of the converted ray is in the box
          if (is_outgoing) {
            // a. the halfedge is outgoing, save the source and the boundary point.
            cropped_face_vertices.push_back(s);
            if (rec.has_on_boundary(sB)) {
              cropped_face_vertices.push_back(sB);
              intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(sB, rec, true), true));
            } else {
              cropped_face_vertices.push_back(tB);
              intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(tB, rec, true), true));
            }
          } else {
            // b. the halfedge is incoming, save the boundary point(s).
            if (rec.has_on_boundary(sB)) {
              cropped_face_vertices.push_back(sB);
              intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(sB, rec, false), false));
            }
            
            if (rec.has_on_boundary(tB)) {
              cropped_face_vertices.push_back(tB);
              intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(tB, rec, false), false));
            }
          }
        } else {
          // 2. the source of the converted ray is not in the box, save the boundary point(s).
          if (is_outgoing) {
            // if the ray is outgoing, the source of intersection is the entry point
            if (rec.has_on_boundary(sB)) {
              cropped_face_vertices.push_back(sB);
              intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(sB, rec, !is_outgoing), !is_outgoing));
            }
            
            if (rec.has_on_boundary(tB)) {
              cropped_face_vertices.push_back(tB);
              intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(tB, rec, is_outgoing), is_outgoing));
            }
          } else {
            // if the ray is incoming, the target of intersection is the entry point
            if (rec.has_on_boundary(tB)) {
              cropped_face_vertices.push_back(tB);
              intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(tB, rec, is_outgoing), is_outgoing));
            }
            if (rec.has_on_boundary(sB)) {
              cropped_face_vertices.push_back(sB);
              intersect_points.push_back(std::make_tuple(cropped_face_vertices.size() - 1, which_side(sB, rec, !is_outgoing), !is_outgoing));
            }
          }
        }
      } else {
        continue; // no intersection, check the next halfedge
      }
    }
  } while (--ec != ec_start);
  
  // insert rectangle vertices as needed
  int exit_side, entry_side;
  std::vector<Point_2> to_insert;
  if (!intersect_points.empty()) {
    // rotate to the right if the first intersection point is not an exit
    if (!std::get<2>(intersect_points[0])) {
      std::rotate(intersect_points.rbegin(), intersect_points.rbegin() + 1, intersect_points.rend());
    }
    
    // find out vertices to insert for each exit point
    for (size_t i = 0; i < intersect_points.size() - 1; i += 2) {
      exit_side = std::get<1>(intersect_points[i]);
      entry_side = std::get<1>(intersect_points[i + 1]);
      if (exit_side != entry_side && entry_side != -1 && exit_side != -1) {
        for (int i = exit_side; i != entry_side; i = ((i - 1) % 4 + 4) % 4) {
          to_insert.push_back(rec.vertex(i));
        }
        vert_insertion.push_back(std::make_pair(std::get<0>(intersect_points[i]), to_insert));
        to_insert.clear();
      }
    }
    
    if (!vert_insertion.empty()) {
      // sort indices into descending order for convenient insertion
      std::sort(vert_insertion.rbegin(), vert_insertion.rend());
      for (auto it = vert_insertion.begin(); it != vert_insertion.end(); ++it) {
        cropped_face_vertices.insert(cropped_face_vertices.begin() + std::get<0>(*it) + 1, (std::get<1>(*it)).begin(), (std::get<1>(*it)).end());
      }
    }
  }
  
  // repeat the first vertex to close the circulation
  if (cropped_face_vertices.size() > 2) {
    cropped_face_vertices.push_back(cropped_face_vertices[0]);
  }
  
  return cropped_face_vertices;
}

// compute the area for each cell of the power diagram of Y
NumericVector Voronoi_Area_2D(const NumericMatrix &X, const NumericVector &w, const NumericVector &b) {
  // load data
  std::vector<Weighted_point_2> wpoints;
  for (int i = 0; i < X.nrow(); i++) {
    Point_2 p(X(i, 0), X(i, 1));
    Weighted_point_2 wp(p, w(i));
    wpoints.push_back(wp);
  }
  RT2 rt(wpoints.begin(), wpoints.end());
  VD vd(rt);
  
  // construct a rectangle
  Iso_rectangle_2 B(b(0), b(1), b(2), b(3));
  
  // iterate through vertices of rt, extract corresponding faces from vd
  int face_index;
  NumericVector voronoi_area(X.nrow());
  std::vector<Point_2> cropped_face_vertices;
  for (auto f = vd.faces_begin(); f != vd.faces_end(); ++f) {
    face_index = std::distance(wpoints.begin(), std::find(wpoints.begin(), wpoints.end(), f->dual()->point()));
    voronoi_area[face_index] = Area_Shoelace(Face_Rect_Crop(f, vd, B));
  }
  
  return voronoi_area;
}

// simple gradient descent for now
// [[Rcpp::export]]
NumericVector Transport_2D(const NumericMatrix &X) {
  NumericVector b = NumericVector::create(0.0, 0.0, 1.0, 1.0); // 2D uniform measure
  std::size_t n = X.nrow();
  double nr = 1 / static_cast<double>(n);
  
  // initial weight
  NumericVector v(n, nr);
  NumericVector w(n);
  NumericVector grad;
  
  // optmize
  for (int i = 0; i < 500; i++) {
    grad = v - Voronoi_Area_2D(X, w, b);
    w = w + 0.02 * grad;
    // Rcout << max(abs(grad)) / min(v) << std::endl;
  }

  return w;
}