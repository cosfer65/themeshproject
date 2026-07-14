
#if 0

#include "curvature.h"
#include "vector.h"
#include "eigen.h"
#include <Windows.h>

using namespace btm;
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
    class CurvatureCalculator {
        btm::MeshExplicit<double>& mesh;


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
        double computeCotangentWeight(size_t v1_index, size_t v2_index) {
            // Find the two triangles that share edge (i,j)
            double cot_sum = 0.0;

            VertexExplicit<double>* v1 = &mesh.vertices[v1_index];
            VertexExplicit<double>* v2 = &mesh.vertices[v2_index];

            // for (meshFace<double>* face : v1->adjacentFaces) {
            for (auto face_index : mesh.vertex_adjacency[v1_index].incident_faces) {
                FaceExplicit<double>* face = &mesh.faces[face_index];
                // Check if this triangle contains both vertices i and j
                bool has_i = false, has_j = false;
                size_t other_vertex = (size_t)-1;

                for (int k = 0; k < 3; ++k) {
                    if (face->v(k) == v1_index) has_i = true;
                    else if (face->v(k) == v2_index) has_j = true;
                    else other_vertex = face->v(k);
                }

                if (has_i && has_j && other_vertex != (size_t)-1) {
                    dvec3 vi = v1->position;
                    dvec3 vj = v2->position;
                    dvec3 vk = mesh.vertices[other_vertex].position;

                    dvec3 e1 = vk - vi;
                    dvec3 e2 = vk - vj;
                    double cos_angle = e1.normalize().dot(e2.normalize());
                    double sin_angle = e1.cross(e2).length() / (e1.length() * e2.length());
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
        void computePrincipalDirections(size_t v_index, std::vector<VertexCurvature>& curvatures)    // input vertex for which to compute principal directions
        {
            VertexExplicit<double>* v = &mesh.vertices[v_index];
            VertexCurvature& vc = curvatures[v_index]; 

            // k1, k2 and dir1, dir2 are references to vc.
            double& k1 = vc.k1;
            double& k2 = vc.k2;
            basevector<double, 3>& dir1 = vc.k1_dir;
            basevector<double, 3>& dir2 = vc.k2_dir;

            // Build tangent frame (t1, t2, n)
            // use vertex normal ni and a safe reference vector ref.
            basevector<double, 3>& pi = v->position;
            basevector<double, 3>& ni = v->normal;

            basevector<double, 3> ref = (std::fabs(ni.x()) > std::fabs(ni.z()))
                ? basevector<double, 3>(-ni.y(), ni.x(), 0.0)
                : basevector<double, 3>(0.0, -ni.z(), ni.y());

            // reate orthonormal tangent vectors t1, t2 so local coordinates can be measured in the tangent plane.
            basevector<double, 3> t1(ref); t1.normalize();
            basevector<double, 3> t2(cross(ni, t1)); t2.normalize();

            // build least - squares normal equations

            // Accumulate least-squares system A^T A x = A^T b
            // Unknowns are shape-operator coefficients: x = [s11, s12, s22]^T.
            double ATA[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
            double ATb[3] = { 0,0,0 };

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
                const basevector<double, 3>& pj = mesh.vertices[nbr_index].position;
                const basevector<double, 3>& nj = mesh.vertices[nbr_index].normal;

                basevector<double, 3> dij = pj - pi;
                double u = dot(dij, t1);
                double v = dot(dij, t2);

                basevector<double, 3> dn = nj - ni;
                double b1 = dot(dn, t1);
                double b2 = dot(dn, t2);

                // Row for b1: [u, v, 0] * x = b1
                {
                    double a0 = u;
                    double a1 = v;
                    double a2 = 0.0;

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
                    double a0 = 0.0;
                    double a1 = u;
                    double a2 = v;

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
            double M[3][4] = {
                { ATA[0][0], ATA[0][1], ATA[0][2], ATb[0] },
                { ATA[1][0], ATA[1][1], ATA[1][2], ATb[1] },
                { ATA[2][0], ATA[2][1], ATA[2][2], ATb[2] }
            };

            // Gaussian elimination
            for (int r = 0; r < 3; ++r) {
                // pivot
                int pivot = r;
                double maxA = std::fabs(M[r][r]);
                for (int rr = r + 1; rr < 3; ++rr) {
                    double val = std::fabs(M[rr][r]);
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
                double diag = M[r][r];
                for (int c = r; c < 4; ++c) M[r][c] /= diag;
                // eliminate below
                for (int rr = r + 1; rr < 3; ++rr) {
                    double factor = M[rr][r];
                    for (int c = r; c < 4; ++c)
                        M[rr][c] -= factor * M[r][c];
                }
            }
            // back substitution
            for (int r = 2; r >= 0; --r) {
                for (int rr = 0; rr < r; ++rr) {
                    double factor = M[rr][r];
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
                double s11 = M[0][3];
                double s12 = M[1][3];
                double s22 = M[2][3];

                // Eigen decomposition of S = [s11 s12; s12 s22]
                double l1, l2, v1x, v1y, v2x, v2y;
                eigenSymmetric2x2<double>(s11, s12, s22, l1, l2, v1x, v1y, v2x, v2y);

                // Map eigenvectors back to 3D
                //  - `d1 = v1x * t1 + v1y * t2`, `d2 = v2x * t1 + v2y * t2`.
                basevector<double, 3> d1 = v1x * t1 + v1y * t2;
                basevector<double, 3> d2 = v2x * t1 + v2y * t2;

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

    public:
        /**
         * @brief Constructs a CurvatureCalculator for the specified mesh.
         *
         * @param _mesh Reference to the MeshExplicit<double> object for which curvatures will be computed.
         */
        CurvatureCalculator(btm::MeshExplicit<double>& _mesh)
            : mesh(_mesh) {}

        /**
         * @brief Calculates the curvature properties of a vertex.
         *
         * This function computes both the Gaussian and Mean curvature for a specified vertex
         * using its geometric properties and those of its neighbors. The calculation
         * is based on the discrete Laplace-Beltrami operator for mean curvature and the
         * angle defect for Gaussian curvature.
         *
         * @param vertex_index The index of the vertex for which to calculate curvature.
         * @return A VertexCurvature struct containing the calculated curvature values.
         */
        VertexCurvature calculate_vertex_curvature(size_t vertex_index) {
            VertexCurvature vc;

            VertexExplicit<double>& v = mesh.vertices[vertex_index];
            // fetch angle defect and mixed area
            double area_mixed = v.voronoi_area;
            double angle_defect = 2.0 * PI<double> -v.angle_sum;
            if (v.is_boundary) {
                // For boundary vertices, the angle defect is different
                angle_defect = PI<double> -v.angle_sum;
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

            // for (meshVertex<double>* neighbor : v->neighbourVertices) {
            for (auto neighbor_index : mesh.vertex_adjacency[vertex_index].neighbor_vertices) {
                double cot_weight = computeCotangentWeight(vertex_index, neighbor_index);
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
        void compute_vertex_curvatures(std::vector<VertexCurvature>& curvatures) {
            // Compute curvatures and fill the curvatures vector
            size_t num_vertices = mesh.vertices.size();
            for (size_t i = 0; i < num_vertices; ++i) {
                // Compute curvature for vertex i and fill vc
                VertexCurvature vc = calculate_vertex_curvature(i);
                curvatures.push_back(vc);
                computePrincipalDirections(i, curvatures);
            }
        }
    };

    /**
     * @brief Computes curvature properties for all vertices in a mesh.
     *
     * This is the main entry point function for computing vertex curvatures in a mesh.
     * It creates a CurvatureCalculator instance and uses it to compute Gaussian curvature,
     * mean curvature, and principal curvature directions for all vertices.
     *
     * @param mesh Reference to the MeshExplicit<double> object containing the mesh geometry.
     *             The mesh must have valid vertex positions, normals, face connectivity,
     *             vertex adjacency information, Voronoi areas, and angle sums pre-computed.
     * @param curvatures Output vector that will be filled with VertexCurvature data.
     *                   Upon return, curvatures[i] contains the curvature information for
     *                   mesh.vertices[i].
     *
     * @return int Returns 1 for success, 0 for failure.
     *
     * @note This function requires that mesh preprocessing has been completed, including:
     *       - Vertex normals computed
     *       - Voronoi areas computed
     *       - Angle sums computed
     *       - Vertex adjacency information populated
     */
    int compute_vertex_curvatures(btm::MeshExplicit<double>& mesh, std::vector<VertexCurvature>& curvatures)
    {
        CurvatureCalculator calculator(mesh);
        calculator.compute_vertex_curvatures(curvatures);

        return 1; // Return 1 for success, 0 for failure
    }
}
#endif
