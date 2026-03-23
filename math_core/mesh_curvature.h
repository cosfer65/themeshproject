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
    T voronoi_based_weighting(mesh<T>& fmesh, size_t v_id, meshFace<T>& f)
    {
        // Initialize area weight to zero
        T w_area = 0.0;

        // Get the vertex object from the mesh
        meshVertex<T>& v = *(fmesh.vertices[v_id]);

        // Get the coordinates of vertex v
        basevector<T, 3>& P = v.position;
        // Get the other two vertices of the triangle face f
        std::pair<meshVertex<T>*, meshVertex<T>*> other = f.getIncidentVertices(v_id);
        if (!other.first || !other.second) {
            // If the vertex is not part of the face, return zero area contribution
            return w_area;
        }
        basevector<T, 3>& Q = other.first->position;
        basevector<T, 3>& R = other.second->position;

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

    /*
    Purpose
        worldCoordsToFaceCoords converts the 3D coordinates of a triangle's vertices (in world space) into 2D coordinates in the face's
        local tangent plane, defined by basis vectors u_f and v_f.
        This is a standard step in curvature computation: you want to work in a 2D parameterization of the triangle (u,v) instead of
        raw 3D coordinates.
    Parameters
        -const dvec3& u_f, const dvec3& v_f
        Orthonormal tangent basis for the face. Together with the face normal, they form a local 3D frame, but only the tangent
        u_f, v_f are used here.
        -const dvec3* vertices
        Pointer to an array of three 3D vertex positions of the face: vertices[0], vertices[1], vertices[2].
        -const dvec3& centroid
        The 3D centroid of the face. Used as origin of the local coordinate system.
        -dvec2* results
    Output array of three 2D coordinates (u,v) for each vertex in the face local frame:
        -results[0] -> 2D coords of vertices[0]
        -results[1] -> 2D coords of vertices[1]
        -results[2] -> 2D coords of vertices[2]
    */
    template <typename T>
    void worldCoordsToFaceCoords(const basevector<T, 3>& u_f, const basevector<T, 3>& v_f, const
        basevector<T, 3>* vertices, const basevector<T, 3>& centroid, basevector<T, 2>* results)
    {
        // Project vertices onto plane defined by centroid
        basevector<T, 3> proj_n1 = vertices[0] - centroid;
        basevector<T, 3> proj_n2 = vertices[1] - centroid;
        basevector<T, 3> proj_n3 = vertices[2] - centroid;

        // Now proj_n1, proj_n2, proj_n3 are vectors from the centroid to each vertex.
        // So the centroid effectively becomes the(0, 0) origin in the local coordinates.

        // u,v in a flat array for better performance. Conceptually this is a 2x3 matrix:
        T uvf[] = { u_f.x(), u_f.y(), u_f.z(), v_f.x(), v_f.y(), v_f.z() };

        // Unrolled matrix-vector multiplication for better performance on fixed 2x3 operation
        results[0][0] = uvf[0] * proj_n1.x() + uvf[1] * proj_n1.y() + uvf[2] * proj_n1.z();
        results[0][1] = uvf[3] * proj_n1.x() + uvf[4] * proj_n1.y() + uvf[5] * proj_n1.z();

        results[1][0] = uvf[0] * proj_n2.x() + uvf[1] * proj_n2.y() + uvf[2] * proj_n2.z();
        results[1][1] = uvf[3] * proj_n2.x() + uvf[4] * proj_n2.y() + uvf[5] * proj_n2.z();

        results[2][0] = uvf[0] * proj_n3.x() + uvf[1] * proj_n3.y() + uvf[2] * proj_n3.z();
        results[2][1] = uvf[3] * proj_n3.x() + uvf[4] * proj_n3.y() + uvf[5] * proj_n3.z();
    }

    /*
    Purpose
        convert3dTo2dCoords prepares data for computing the second fundamental form of a triangular mesh face.
        It transforms 3D normal differences into a 2D coordinate system and packages them with 2D position
        differences for later least-squares fitting.

        This function is a key step in discrete curvature computation: it converts the 3D geometric changes
        (how normals vary across the triangle) into 2D measurements that can be used to fit the second
        fundamental form coefficients [ee, ff, gg].
    Parameters
        -const dmat2x3& uv_map
        A 2x3 transformation matrix that projects 3D vectors into the face's local 2D coordinate system.
        Row-major layout: [u_x, u_y, u_z, v_x, v_y, v_z]
        where (u, v) define the orthonormal tangent basis of the face plane.

        -const dvec2(&nodes_local)[3]
        Array of three 2D vertex positions in the face's local coordinate system.
        These represent the triangle's vertices projected onto the tangent plane:
        nodes_local[0]: vertex 0 in 2D face coordinates
        nodes_local[1]: vertex 1 in 2D face coordinates
        nodes_local[2]: vertex 2 in 2D face coordinates

        -const dvec3(&normals)[3]
        Array of three 3D vertex normals at the triangle's vertices:
        normals[0]: normal at vertex 0
        normals[1]: normal at vertex 1
        normals[2]: normal at vertex 2

        -double* result
        Output array of 12 doubles organized as follows:
        [0-5]:   Projected normal differences (6 floats: 3 edges × 2 components)
                 [0-1]: diff(n2-n1) in 2D
                 [2-3]: diff(n0-n2) in 2D
                 [4-5]: diff(n1-n0) in 2D
        [6-11]:  Edge position differences (6 floats: 3 edges × 2 components)
                 [6-7]:   diff(pos2-pos1) in 2D
                 [8-9]:   diff(pos0-pos2) in 2D
                 [10-11]: diff(pos1-pos0) in 2D

    Mathematical Background
        The output feeds into leastSquares(), which solves the linear system:
            E * [ee, ff, gg]^T = N
        where:
        - N (result[0-5]):  How normals change across edges (measures surface bending)
        - E (result[6-11]): How positions change across edges (edge vectors)
        - [ee, ff, gg]:     Second fundamental form coefficients describing surface curvature

        The second fundamental form captures how the surface normal changes, which is the geometric
        definition of curvature. The least-squares fit finds the best quadratic approximation to
        this normal variation.
    */
    template <typename T>
    void convert3dTo2dCoords(const basematrix<T, 2, 3>& uv_map, const basevector<T, 2>(&nodes_local)[3], const basevector<T, 3>(&normals)[3], T* result)
    {
        // Calculates how vertex normals change across the triangle's three edges. This captures how the surface is bending.
        basevector<T, 3> diff_n2n1(normals[2] - normals[1]);
        basevector<T, 3> diff_n0n2(normals[0] - normals[2]);
        basevector<T, 3> diff_n1n0(normals[1] - normals[0]);

        // Project these 3D normal differences into the 2D local coordinate system defined by uv_map
        basevector<T, 3> diff_n2_n1_3d(diff_n2n1.x(), diff_n2n1.y(), diff_n2n1.z());
        basevector<T, 3> diff_n0_n2_3d(diff_n0n2.x(), diff_n0n2.y(), diff_n0n2.z());
        basevector<T, 3> diff_n1_n0_3d(diff_n1n0.x(), diff_n1n0.y(), diff_n1n0.z());

        basevector<T, 2> diff_n2_n1_2d(0);
        basevector<T, 2> diff_n0_n2_2d(0);
        basevector<T, 2> diff_n1_n0_2d(0);
        // Unrolled matrix-vector multiplication for 2x3 * 3x1 operations
        // uv_map is 2x3 row-major: [u_x, u_y, u_z, v_x, v_y, v_z]
        diff_n2_n1_2d[0] = uv_map[0] * diff_n2_n1_3d.x() + uv_map[1] * diff_n2_n1_3d.y() + uv_map[2] * diff_n2_n1_3d.z();
        diff_n2_n1_2d[1] = uv_map[3] * diff_n2_n1_3d.x() + uv_map[4] * diff_n2_n1_3d.y() + uv_map[5] * diff_n2_n1_3d.z();

        diff_n0_n2_2d[0] = uv_map[0] * diff_n0_n2_3d.x() + uv_map[1] * diff_n0_n2_3d.y() + uv_map[2] * diff_n0_n2_3d.z();
        diff_n0_n2_2d[1] = uv_map[3] * diff_n0_n2_3d.x() + uv_map[4] * diff_n0_n2_3d.y() + uv_map[5] * diff_n0_n2_3d.z();

        diff_n1_n0_2d[0] = uv_map[0] * diff_n1_n0_3d.x() + uv_map[1] * diff_n1_n0_3d.y() + uv_map[2] * diff_n1_n0_3d.z();
        diff_n1_n0_2d[1] = uv_map[3] * diff_n1_n0_3d.x() + uv_map[4] * diff_n1_n0_3d.y() + uv_map[5] * diff_n1_n0_3d.z();

        // Also compute edge differences in local 2D coordinates
        basevector<T, 2> diff_nodes_n2n1 = nodes_local[2] - nodes_local[1];
        basevector<T, 2> diff_nodes_n0n2 = nodes_local[0] - nodes_local[2];
        basevector<T, 2> diff_nodes_n1n0 = nodes_local[1] - nodes_local[0];

        // Store results in the output array
        // Layout:
        // [0-5]: Normal difference (6 floats: 3 edges x 2 components)
        // [6-11]: Edge coordinate differences (6 floats: 3 edges x 2 components)
        // The output feeds into leastSquares(), which solves
        // E * [ee, ff, gg] ^ T = N
        // N(result[0 - 5]) : How normals change
        // E(result[6 - 11]) : How positions change
        // [ee, ff, gg] : Second fundamental form coefficients describing surface curvature
        result[0] = diff_n2_n1_2d[0];
        result[1] = diff_n2_n1_2d[1];
        result[2] = diff_n0_n2_2d[0];
        result[3] = diff_n0_n2_2d[1];
        result[4] = diff_n1_n0_2d[0];
        result[5] = diff_n1_n0_2d[1];
        result[6] = diff_nodes_n2n1[0];
        result[7] = diff_nodes_n2n1[1];
        result[8] = diff_nodes_n0n2[0];
        result[9] = diff_nodes_n0n2[1];
        result[10] = diff_nodes_n1n0[0];
        result[11] = diff_nodes_n1n0[1];
    }


}
