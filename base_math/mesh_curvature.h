// core functions from Discrete Differential Geometry used to compute curvature on meshes. 
#include <iostream>
#include <string>

#include "base_definitions.h"

#include "vector.h"
#include "geometry.h"
#include "algebra.h"

#include "mesh.h"

namespace base_math {
    //  This function computes the Voronoi area weight of a vertex in a triangle of a half-edge mesh. 
    //  It’s a standard piece used in discrete differential geometry (e.g., mean curvature, Laplace-Beltrami).
    //  -  mesh: the half-edge mesh.
    //  -  v_id: index of the vertex in the mesh whose Voronoi area contribution (from this face) we are computing.
    //  -  f: the triangular face that contains the vertex v_id.
    //  -  Return value: a scalar T which is the area weight assigned to vertex v_id from triangle f.
    template <typename T>
    T voronoi_based_weighting(half_edge_mesh<T>& mesh, size_t v_id, face<T>& f)
    {
        // Initialize area weight to zero
        T w_area = 0.0;

        // Get the vertex object from the mesh
        vertex<T>& v = *(mesh.vertices[v_id]);

        // Get the coordinates of vertex v
        basevector<T, 3>& P = v.coords;
        // Get the other two vertices of the triangle face f
        std::pair<size_t, size_t> other = f.get_other_vertices(v_id);
        basevector<T, 3>& Q = mesh.vertices[other.first]->coords;
        basevector<T, 3>& R = mesh.vertices[other.second]->coords;

        // Edge vectors
        basevector<T, 3> PR = R - P;
        basevector<T, 3> PQ = Q - P;
        basevector<T, 3> QR = R - Q;

        // Total area of this triangle
        T area_T = f.area;

        // Check if triangle is obtuse
        // is_obtuse checks whether the triangle has an obtuse angle (angle > 90°) at any vertex
        bool is_obtuse = PQ.dot(PR) < 0 || (-PQ).dot(QR) < 0 || (-PR).dot(-QR) < 0;

        // For an angle at vertex P, cosine is proportional to PQ* PR.
        // If PQ.dot(PR) < 0, angle at P is obtuse.
        // 
        // For angle at Q :
        // Two edges from Q are - PQ(from Q to P) and QR(from Q to R).
        // If(-PQ).dot(QR) < 0, angle at Q is obtuse.
        // 
        // For angle at R :
        // Two edges from R are - PR(from R to P) and -QR(from R to Q).
        // If(-PR).dot(-QR) < 0, angle at R is obtuse.
        if (!is_obtuse) {
            // 1. Voronoi Area using cotangent weights
            // In a non-obtuse triangle, the Voronoi region of each vertex can be nicely defined inside the triangle. 
            // The code is computing the Voronoi area around vertex P using cotangent weights.
            // Key points:
            // -  PQ.dot(-QR) / PQ.cross(-QR).length():
            // -    This is cot(angle at Q):
            // -      dot(a, b) = |a||b| cos(θ)
            // -      |a * b| = |a||b| sin(θ)
            // -      So dot / |cross| = cos(θ) / sin(θ) = cot(θ).
            // -  PR.dot(QR) / PR.cross(QR).length():
            // -    This is cot(angle at R).
            // 
            // -  PQ.dot(PQ) is |PQ|^2 (squared edge length).
            // -  PR.dot(PR) is |PR|^2.
            // -  The formula:
            //      [ A_Voronoi(P) = (|PQ|^2 * cot(angle R) + |PR|^2 * cot(angle Q))/8 ]
            // This is a standard mixed Voronoi area formula used in discrete differential geometry 
            // (e.g., from Meyer et al., “Discrete Differential-Geometry Operators for Triangulated 2-Manifolds”).
            // This gives the portion of the triangle’s area that is assigned to vertex P when the triangle is acute.
            T cot_Q = PQ.dot(-QR) / PQ.cross(-QR).length();
            T cot_R = PR.dot(QR) / PR.cross(QR).length();
            w_area = (PQ.dot(PQ) * cot_R + PR.dot(PR) * cot_Q) / T(8.0);
        }
        else {
            // 2. Obtuse case
            // For obtuse triangles the standard Voronoi partition (cutting by perpendicular bisectors) would extend outside the triangle, 
            // so the usual formula is not appropriate. The typical remedy:
            // -  If the obtuse angle is at the vertex P:
            // -    PQ.dot(PR) < 0 means angle at P is obtuse.
            // -    Then the whole triangle area is split so that the obtuse vertex gets half the triangle area:
            // -      w_area = area_T / 2.0;
            // -  Otherwise (obtuse at Q or R, but not at P):
            // -    The vertex P gets one quarter of the triangle area:
            // -      w_area = area_T / 4.0;
            // This matches standard “mixed area” definitions used in many mesh processing algorithms.
            if (PQ.dot(PR) < 0) {
                // Angle at P is obtuse
                w_area = area_T / T(2.0);
            }
            else {
                // Angle at Q or R is obtuse
                w_area = area_T / T(4.0);
            }
        }
        return w_area;
    }
}
