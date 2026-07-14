// Curvature.h
#pragma once

#include <vector>
#include <fstream>
#include "vector.h"
#include "eigen.h"
#include "mesh_explicit.h"

namespace curvature
{

    /**
     * @brief Private mesh curvature calculator class for computing discrete curvature properties.
     *
     * This class provides methods to compute various curvature properties for vertices in an
     * explicit mesh representation, including Gaussian curvature, mean curvature, and principal
     * curvature directions. It uses discrete differential geometry techniques such as cotangent
     * weighting, the discrete Laplace-Beltrami operator, angle defects, and least-squares fitting.
     *
     * @note This class operates on a reference to an existing MeshExplicit<double> object.
     */
    template <typename T>
    class CurvatureCalculator {
        btm::MeshExplicit<T>& mesh;

        /**
         * @brief Computes the cotangent weight for an edge between two vertices.
         *
         * This function calculates the cotangent weight used in discrete differential geometry
         * for mesh processing operations such as mean curvature computation or the discrete
         * Laplace-Beltrami operator. The weight is computed as half the sum of cotangents of
         * the angles opposite to the edge in the triangles sharing that edge.
         *
         * For an edge (v1, v2), the cotangent weight is:
         *   w_ij = 0.5 * (cot(α) + cot(β))
         * where α and β are the angles opposite to the edge in the two incident triangles.
         *
         * @param v1_index Index of the first vertex of the edge
         * @param v2_index Index of the second vertex of the edge
         * @return double The computed cotangent weight for the edge
         *
         * @note The function iterates through faces incident to v1 to find triangles
         *       that contain both vertices.
         * @note If the edge is on the boundary (only one incident triangle), only one
         *       cotangent term will be included in the sum.
         * @note The cotangent is computed as cos(θ)/sin(θ) where the cosine is calculated
         *       using the dot product and sine using the cross product magnitude.
         */
        T computeCotangentWeight(size_t v1_index, size_t v2_index) {
            // Find the two triangles that share edge (i,j)
            T cot_sum = 0.0;

            VertexExplicit<T>* v1 = &mesh.vertices[v1_index];
            VertexExplicit<T>* v2 = &mesh.vertices[v2_index];

            // for (meshFace<double>* face : v1->adjacentFaces) {
            for (auto face_index : mesh.vertex_adjacency[v1_index].incident_faces) {
                FaceExplicit<T>* face = &mesh.faces[face_index];
                // Check if this triangle contains both vertices i and j
                bool has_i = false, has_j = false;
                size_t other_vertex = (size_t)-1;

                for (int k = 0; k < 3; ++k) {
                    if (face->v(k) == v1_index) has_i = true;
                    else if (face->v(k) == v2_index) has_j = true;
                    else other_vertex = face->v(k);
                }

                if (has_i && has_j && other_vertex != (size_t)-1) {
                    basevec3<T> vi = v1->position;
                    basevec3<T> vj = v2->position;
                    basevec3<T> vk = mesh.vertices[other_vertex].position;

                    basevec3<T> e1 = vk - vi;
                    basevec3<T> e2 = vk - vj;
                    T cos_angle = e1.normalize().dot(e2.normalize());
                    T sin_angle = e1.cross(e2).length() / (e1.length() * e2.length());
                    cot_sum += cos_angle / sin_angle;
                }
            }

            return cot_sum / 2.0;
        }

        /**
         * @brief Computes principal curvature directions and values at a given vertex using least-squares fitting.
         *
         * This method calculates the principal curvatures (k1, k2) and their corresponding principal directions
         * at a specific vertex by fitting a quadratic surface to the local neighborhood using least-squares approximation.
         * The method constructs a tangent coordinate system at the vertex, builds a normal equations system from
         * neighbor vertex information, and solves for the shape operator coefficients via eigendecomposition.
         *
         * @param v_index The index of the vertex in the mesh for which to compute principal directions.
         * @param curvatures Reference to a vector storing VertexCurvature data for all mesh vertices.
         *                   The curvature information for v_index will be updated.
         *
         * @details
         * Algorithm steps:
         * 1. **Tangent Frame Construction**: Builds an orthonormal tangent frame (t1, t2, n) at the vertex,
         *    where n is the vertex normal and t1, t2 span the tangent plane.
         *
         * 2. **Least-Squares System Setup**: For each neighboring vertex, computes:
         *    - Position offset projected into tangent plane: (u, v)
         *    - Normal difference projected into tangent plane: (b1, b2)
         *    Constructs two linear constraints per neighbor:
         *    - [u, v, 0] * x = b1
         *    - [0, u, v] * x = b2
         *    where x = [s11, s12, s22]^T are the unknown shape operator coefficients.
         *    Accumulates A^T*A and A^T*b matrices directly.
         *
         * 3. **Linear System Solution**: Solves the 3x3 system ATA * x = ATb using Gaussian elimination
         *    with partial pivoting to obtain the shape operator coefficients.
         *
         * 4. **Eigendecomposition**: Reconstructs the symmetric 2x2 shape operator matrix S = [s11 s12; s12 s22]
         *    and computes its eigenvalues (principal curvatures) and eigenvectors (principal directions)
         *    via eigenSymmetric2x2().
         *
         * 5. **Result Storage**: Maps 2D eigenvectors back to 3D tangent directions, normalizes them,
         *    and stores results in the curvature structure. Ensures k1 <= k2 by swapping if necessary.
         *
         * @note Special cases:
         * - If the vertex has fewer than 3 neighbors, returns zero curvature with basis directions set to t1, t2.
         * - If the linear system is degenerate (pivot < 1e-12), returns zero curvature with basis directions.
         *
         * @see eigenSymmetric2x2() for 2x2 symmetric matrix eigendecomposition
         * @see VertexCurvature for the curvature data structure
         */
        void computePrincipalDirections(size_t v_index, btm::fast_vec<VertexCurvature<T>>& curvatures)    // input vertex for which to compute principal directions
        {
            VertexExplicit<T>* v = &mesh.vertices[v_index];
            VertexCurvature<T>& vc = curvatures[v_index];

            // k1, k2 and dir1, dir2 are references to vc.
            T& k1 = vc.k1;
            T& k2 = vc.k2;
            basevec3<T>& dir1 = vc.k1_dir;
            basevec3<T>& dir2 = vc.k2_dir;

            // Build tangent frame (t1, t2, n)
            // use vertex normal ni and a safe reference vector ref.
            basevec3<T>& pi = v->position;
            basevec3<T>& ni = v->normal;

            basevec3<T> ref = (std::fabs(ni.x()) > std::fabs(ni.z()))
                ? basevec3<T>(-ni.y(), ni.x(), 0.0)
                : basevec3<T>(0.0, -ni.z(), ni.y());

            // reate orthonormal tangent vectors t1, t2 so local coordinates can be measured in the tangent plane.
            basevec3<T> t1(ref); t1.normalize();
            basevec3<T> t2(cross(ni, t1)); t2.normalize();

            // build least - squares normal equations

            // Accumulate least-squares system A^T A x = A^T b
            // Unknowns are shape-operator coefficients: x = [s11, s12, s22]^T.
            T ATA[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
            T ATb[3] = { 0,0,0 };

            VertexAdjacency& va = mesh.vertex_adjacency[v_index];
            fast_vec<std::uint32_t>& nbrs = va.neighbor_vertices;

            // handle underconstrained case
            // If fewer than 3 neighbors, return zero curvature and basis directions.
            if (nbrs.size() < 3) {
                // Not enough neighbors to fit a tensor
                k1 = k2 = 0.0;
                dir1 = t1;
                dir2 = t2;
                return;
            }

            // For each neighbor :
            //  Position offset projected into tangent plane : (u, v).
            // Normal difference projected into tangent plane : (b1, b2).
            // Adds two linear constraints :
            // [u, v, 0] * x = b1
            // [0, u, v] * x = b2
            //  Accumulates A^T A (ATA) and A^T b (ATb) directly.
            for (std::uint32_t nbr_index : nbrs) {
                const basevec3<T>& pj = mesh.vertices[nbr_index].position;
                const basevec3<T>& nj = mesh.vertices[nbr_index].normal;

                basevec3<T> dij = pj - pi;
                T u = dot(dij, t1);
                T v = dot(dij, t2);

                basevec3<T> dn = nj - ni;
                T b1 = dot(dn, t1);
                T b2 = dot(dn, t2);

                // Row for b1: [u, v, 0] * x = b1
                {
                    T a0 = u;
                    T a1 = v;
                    T a2 = 0.0;

                    ATA[0][0] += a0 * a0;
                    ATA[0][1] += a0 * a1;
                    ATA[0][2] += a0 * a2;
                    ATA[1][1] += a1 * a1;
                    ATA[1][2] += a1 * a2;
                    ATA[2][2] += a2 * a2;

                    ATb[0] += a0 * b1;
                    ATb[1] += a1 * b1;
                    ATb[2] += a2 * b1;
                }

                // Row for b2: [0, u, v] * x = b2
                {
                    T a0 = 0.0;
                    T a1 = u;
                    T a2 = v;

                    ATA[0][0] += a0 * a0;
                    ATA[0][1] += a0 * a1;
                    ATA[0][2] += a0 * a2;
                    ATA[1][1] += a1 * a1;
                    ATA[1][2] += a1 * a2;
                    ATA[2][2] += a2 * a2;

                    ATb[0] += a0 * b2;
                    ATb[1] += a1 * b2;
                    ATb[2] += a2 * b2;
                }
            } // (for (std::uint32_t nbr_index : nbrs))

            // Symmetrize ATA
            ATA[1][0] = ATA[0][1];
            ATA[2][0] = ATA[0][2];
            ATA[2][1] = ATA[1][2];

            // Solve 3x3 system ATA * x = ATb (simple Gaussian elimination)
            // Solves 3x3 system
            // Solves ATA * x = ATb via Gaussian elimination with partial pivoting.
            T M[3][4] = {
                { ATA[0][0], ATA[0][1], ATA[0][2], ATb[0] },
                { ATA[1][0], ATA[1][1], ATA[1][2], ATb[1] },
                { ATA[2][0], ATA[2][1], ATA[2][2], ATb[2] }
            };

            // Gaussian elimination
            for (int r = 0; r < 3; ++r) {
                // pivot
                int pivot = r;
                T maxA = std::fabs(M[r][r]);
                for (int rr = r + 1; rr < 3; ++rr) {
                    T val = std::fabs(M[rr][r]);
                    if (val > maxA) { maxA = val; pivot = rr; }
                }
                //  If pivot is too small(< 1e-12), treats system as degenerate and returns zeros.
                if (maxA < 1e-12) {
                    // Degenerate system
                    k1 = k2 = 0.0;
                    dir1 = t1;
                    dir2 = t2;
                    return;
                }
                if (pivot != r) {
                    for (int c = 0; c < 4; ++c)
                        std::swap(M[r][c], M[pivot][c]);
                }
                // normalize row r
                T diag = M[r][r];
                for (int c = r; c < 4; ++c) M[r][c] /= diag;
                // eliminate below
                for (int rr = r + 1; rr < 3; ++rr) {
                    T factor = M[rr][r];
                    for (int c = r; c < 4; ++c)
                        M[rr][c] -= factor * M[r][c];
                }
            }
            // back substitution
            for (int r = 2; r >= 0; --r) {
                for (int rr = 0; rr < r; ++rr) {
                    T factor = M[rr][r];
                    M[rr][r] = 0.0;
                    M[rr][3] -= factor * M[r][3];
                }
            }

            {
                // Extract principal curvatures**
                //   -Reconstruct symmetric 2x2 matrix `S = [s11 s12; s12 s22]`.
                //   - Call `eigenSymmetric2x2(...)`:
                //   - eigenvalues → principal curvatures(`k1`, `k2` in this code’s assignment order),
                //   -eigenvectors → 2D principal directions.
                T s11 = M[0][3];
                T s12 = M[1][3];
                T s22 = M[2][3];

                // Eigen decomposition of S = [s11 s12; s12 s22]
                T l1, l2, v1x, v1y, v2x, v2y;
                eigenSymmetric2x2<T>(s11, s12, s22, l1, l2, v1x, v1y, v2x, v2y);

                // Map eigenvectors back to 3D
                //  - `d1 = v1x * t1 + v1y * t2`, `d2 = v2x * t1 + v2y * t2`.
                basevec3<T> d1 = v1x * t1 + v1y * t2;
                basevec3<T> d2 = v2x * t1 + v2y * t2;

                // Store
                // - Normalize and store into `k_1_dir`, `k_2_dir`.
                k1 = l1;
                k2 = l2;
                dir1 = d1.normalize();
                dir2 = d2.normalize();
                // k_1 is assigned the smaller curvature, k_2 the larger, so swap if necessary.
                if (fabs(k1) > fabs(k2)) {
                    // swap to ensure k1 <= k2
                    std::swap(k1, k2);
                    std::swap(dir1, dir2);
                }
            }
        }
        std::fstream logFile_gaus; // Optional log file for debugging or output
        bool logFileOpen_gaus = false; // Flag to indicate if the log file is open
    public:
        /**
         * @brief Constructs a CurvatureCalculator for the specified mesh.
         *
         * @param _mesh Reference to the MeshExplicit<T> object for which curvatures will be computed.
         */
        CurvatureCalculator(btm::MeshExplicit<T>& _mesh)
            : mesh(_mesh) {}

        /**
         * @brief Calculates the curvature properties of a vertex.
         *
         * This function computes both the Gaussian and Mean curvature for a specified vertex
         * using its geometric properties and those of its neighbors. The calculation
         * is based on the discrete Laplace-Beltrami operator for mean curvature and the
         * angle defect for Gaussian curvature.
         *
         * Gaussian curvature is computed using the angle defect method:
         * - For interior vertices: K = (2π - Σθ) / A_mixed
         * - For boundary vertices: K = (π - Σθ) / A_mixed
         *
         * Mean curvature is computed using the discrete Laplace-Beltrami operator:
         * H = (1/2A) * Σ cot(α_ij) * (p_j - p_i)
         *
         * @param vertex_index The index of the vertex for which to calculate curvature.
         * @return A VertexCurvature struct containing the calculated curvature values.
         */
        VertexCurvature<T> calculate_vertex_curvature(size_t vertex_index) {
            VertexCurvature<T> vc;

            VertexExplicit<T>& v = mesh.vertices[vertex_index];
            // fetch angle defect and mixed area
            T area_mixed = v.voronoi_area;
            T angle_defect = 2.0 * PI<T> -v.angle_sum;
            bool is_boundary = v.is_boundary;
            if (v.is_boundary) {
                // For boundary vertices, the angle defect is different
                angle_defect = PI<T> -v.angle_sum;
            }
            if (logFileOpen_gaus) {
                logFile_gaus << "Vertex " << vertex_index << ": is_boundary = " << is_boundary
                    << ", angle_sum = " << v.angle_sum
                    // << ", angle_defect = " << angle_defect
                    << ", area_mixed = " << area_mixed << std::endl;
            }

            vc.gaussian = angle_defect / area_mixed;
            vc.absGaussCurvature = fabs(vc.gaussian);

            if (vc.absGaussCurvature < 0.02) {
                vc.signGauss = 0;
            }
            else {
                vc.signGauss = (vc.gaussian > 0) ? 1 : ((vc.gaussian < 0) ? -1 : 0);
            }

            // Compute mean curvature using discrete Laplace-Beltrami operator
            dvec3 laplace(0, 0, 0);

            for (auto neighbor_index : mesh.vertex_adjacency[vertex_index].neighbor_vertices) {
                T cot_weight = computeCotangentWeight(vertex_index, neighbor_index);
                laplace += cot_weight * (mesh.vertices[neighbor_index].position - v.position);
            }

            // Mean curvature = 0.5 * ||Laplace-Beltrami||
            vc.mean = 0.5 * laplace.dot(v.normal) / area_mixed;
            vc.meanCurvatureDir = v.normal * vc.mean * 2;

            return vc;
        }

        /**
         * @brief Computes curvature properties for all vertices in the mesh.
         *
         * This method iterates through all vertices in the mesh and computes their curvature
         * properties including Gaussian curvature, mean curvature, and principal curvature
         * directions. Results are stored in the provided vector.
         *
         * @param curvatures Output vector to be filled with VertexCurvature data for each vertex.
         *                   The vector will be populated with one entry per mesh vertex, where
         *                   curvatures[i] corresponds to mesh.vertices[i].
         *
         * @note The curvatures vector is filled sequentially, so it should be empty or properly
         *       sized before calling this method.
         * @note This method calls both calculate_vertex_curvature() and computePrincipalDirections()
         *       for each vertex to compute complete curvature information.
         */
        void compute_vertex_curvatures() {
            btm::fast_vec<VertexCurvature<T>>& curvatures = mesh.vertex_curvatures;
            curvatures.clear();
            // logFile_gaus.open("c:\\temp\\curvature_log.txt", std::ios::out);
            logFileOpen_gaus = logFile_gaus.is_open();
            // Compute curvatures and fill the curvatures vector
            size_t num_vertices = mesh.vertices.size();
            for (size_t i = 0; i < num_vertices; ++i) {
                // Compute curvature for vertex i and fill vc
                VertexCurvature<T> vc = calculate_vertex_curvature(i);
                curvatures.push_back(vc);
                computePrincipalDirections(i, curvatures);
            }
        }
    };

    // Compute per-vertex curvature for a triangle mesh.
    // Assumes Mesh provides positions, faces, and adjacency.
    template <typename T>
    int compute_vertex_curvatures(btm::MeshExplicit<T>& mesh) {
        CurvatureCalculator<T> calculator(mesh);
        calculator.compute_vertex_curvatures();

        return 1; // Return 1 for success, 0 for failure
    }
}
