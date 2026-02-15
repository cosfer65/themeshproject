#include "base_definitions.h"

#include "vector.h"

#include "mesh.h"
#include "mesh_curvature.h"

using namespace base_math;

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
static void worldCoordsToFaceCoords(const dvec3& u_f, const dvec3& v_f, const dvec3* vertices, const dvec3& centroid, dvec2* results)
{
    // Project vertices onto plane defined by centroid
    dvec3 proj_n1 = vertices[0] - centroid;
    dvec3 proj_n2 = vertices[1] - centroid;
    dvec3 proj_n3 = vertices[2] - centroid;
    // Now proj_n1, proj_n2, proj_n3 are vectors from the centroid to each vertex.
    // So the centroid effectively becomes the(0, 0) origin in the local coordinates.

    // u,v in a flat array for better performance. Conceptually this is a 2x3 matrix:
    double uvf[] = { u_f.x(), u_f.y(), u_f.z(), v_f.x(), v_f.y(), v_f.z() };

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
    [0-5]:   Projected normal differences (6 floats: 3 edges √ó 2 components)
             [0-1]: diff(n2-n1) in 2D
             [2-3]: diff(n0-n2) in 2D
             [4-5]: diff(n1-n0) in 2D
    [6-11]:  Edge position differences (6 floats: 3 edges √ó 2 components)
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
static void convert3dTo2dCoords(const dmat2x3& uv_map, const dvec2(&nodes_local)[3], const dvec3(&normals)[3], double* result)
{
    // Calculates how vertex normals change across the triangle's three edges. This captures how the surface is bending.
    dvec3 diff_n2n1(normals[2] - normals[1]);
    dvec3 diff_n0n2(normals[0] - normals[2]);
    dvec3 diff_n1n0(normals[1] - normals[0]);

    // Project these 3D normal differences into the 2D local coordinate system defined by uv_map
    dvec3 diff_n2_n1_3d(diff_n2n1.x(), diff_n2n1.y(), diff_n2n1.z());
    dvec3 diff_n0_n2_3d(diff_n0n2.x(), diff_n0n2.y(), diff_n0n2.z());
    dvec3 diff_n1_n0_3d(diff_n1n0.x(), diff_n1n0.y(), diff_n1n0.z());

    dvec2 diff_n2_n1_2d(0);
    dvec2 diff_n0_n2_2d(0);
    dvec2 diff_n1_n0_2d(0);

    // Unrolled matrix-vector multiplication for 2x3 * 3x1 operations
    // uv_map is 2x3 row-major: [u_x, u_y, u_z, v_x, v_y, v_z]
    diff_n2_n1_2d[0] = uv_map[0] * diff_n2_n1_3d.x() + uv_map[1] * diff_n2_n1_3d.y() + uv_map[2] * diff_n2_n1_3d.z();
    diff_n2_n1_2d[1] = uv_map[3] * diff_n2_n1_3d.x() + uv_map[4] * diff_n2_n1_3d.y() + uv_map[5] * diff_n2_n1_3d.z();

    diff_n0_n2_2d[0] = uv_map[0] * diff_n0_n2_3d.x() + uv_map[1] * diff_n0_n2_3d.y() + uv_map[2] * diff_n0_n2_3d.z();
    diff_n0_n2_2d[1] = uv_map[3] * diff_n0_n2_3d.x() + uv_map[4] * diff_n0_n2_3d.y() + uv_map[5] * diff_n0_n2_3d.z();

    diff_n1_n0_2d[0] = uv_map[0] * diff_n1_n0_3d.x() + uv_map[1] * diff_n1_n0_3d.y() + uv_map[2] * diff_n1_n0_3d.z();
    diff_n1_n0_2d[1] = uv_map[3] * diff_n1_n0_3d.x() + uv_map[4] * diff_n1_n0_3d.y() + uv_map[5] * diff_n1_n0_3d.z();

    // Also compute edge differences in local 2D coordinates
    dvec2 diff_nodes_n2n1 = nodes_local[2] - nodes_local[1];
    dvec2 diff_nodes_n0n2 = nodes_local[0] - nodes_local[2];
    dvec2 diff_nodes_n1n0 = nodes_local[1] - nodes_local[0];

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

/*
Purpose
    Assemble the per-face data needed to estimate the second fundamental form on a single triangle.

High-level Description
    For the given triangle face `cur_face`, this function:
    1. Collects the 3D vertex positions and transforms them into 2D coordinates in the local
       face tangent frame spanned by `(u_f, v_f)` using `worldCoordsToFaceCoords`.
    2. Gathers the three vertex normals of the face.
    3. Builds the 2x3 matrix `uv_map` that maps 3D vectors into the 2D face coordinate system.
    4. Calls `convert3dTo2dCoords` to:
       - Project normal differences into 2D.
       - Compute edge vectors in 2D.
       - Pack everything into the `results` array, which will later be used by `leastSquares`
         to recover the second fundamental form coefficients `[e, f, g]` for this face.

Parameters
    - mesh
        The half-edge mesh that owns the face. Used to access vertex positions and normals
        via the face's vertex indices (`cur_face.v1`, `cur_face.v2`, `cur_face.v3`).
    - cur_face
        The current triangle face for which second fundamental form data is being prepared.
        Its `center` is used as the origin for the local 2D parameterization.
    - u_f, v_f
        Orthonormal tangent basis vectors of the face in world space.
        They define the 2D local (u, v) coordinate system of the face plane.
    - results
        Pointer to an array of 12 doubles.
        On output, it contains:
          [0-5]   Projected normal differences in 2D (three edges √ó two components).
          [6-11]  Edge position differences in 2D (three edges √ó two components).
        This layout matches the expectations of `convert3dTo2dCoords` and `leastSquares`,
        which solve for the second fundamental form coefficients of the face.

Notes
    - `elemCentroid` is used as the 3D origin for the local 2D parameterization of the face.
    - `world_coords` holds the 3D positions of the triangle's three vertices.
    - `local_coords` holds the corresponding 2D coordinates after projection into the
      face's tangent plane.
    - `uv_map` encodes the same basis (`u_f`, `v_f`) as a 2x3 matrix in row-major order:
      [u_f.x, u_f.y, u_f.z, v_f.x, v_f.y, v_f.z].
    - The normals passed to `convert3dTo2dCoords` are the three vertex normals of `cur_face`,
      which are used to measure how the normal field bends across the face.
*/
static void prepareSecondFundamentalCoefficients(half_edge_mesh<double>& mesh, const face<double>& cur_face, const dvec3& u_f, const dvec3& v_f, double* results)
{
    dvec3 elemCentroid = cur_face.center;
    dvec3 world_coords[3] = { mesh.vertices[cur_face.v1]->coords, mesh.vertices[cur_face.v2]->coords,mesh.vertices[cur_face.v3]->coords };
    dvec2 local_coords[3];
    worldCoordsToFaceCoords(u_f, v_f, world_coords, elemCentroid, local_coords);
    const dvec3& node_normal_1 = mesh.vertices[cur_face.v1]->normal;
    const dvec3& node_normal_2 = mesh.vertices[cur_face.v2]->normal;
    const dvec3& node_normal_3 = mesh.vertices[cur_face.v3]->normal;
    dmat2x3 uv_map({ u_f.x(), u_f.y(), u_f.z(), v_f.x(), v_f.y(), v_f.z() });
    convert3dTo2dCoords(uv_map, local_coords, { node_normal_1, node_normal_2, node_normal_3 }, results);
}

/*
project_uv_to_face_2d
This function performs a coordinate transformation by projecting two 3D basis vectors (u_p and v_p)
from vertex space into a 2D element coordinate system.

High-Level Purpose
The function transforms vertex plane basis vectors into the 2D parametric space of a triangular face element.
This is essential for computing curvature, as it allows you to express geometric quantities consistently
between different coordinate frames.

How It Works
Parameters
-   Map: A 2x3 transformation matrix (stored as a flat double array of 6 elements)
-   u_p, v_p: 3D basis vectors in vertex space
-   res1, res2: Output 2x1 vectors (2 floats each) representing u_p and v_p in 2D element space

The Math
The function performs two matrix-vector multiplications:
res1 = map_mat x u_p    (2x3 matrix times 3x1 vector = 2x1 result)
res2 = map_mat x v_p
*/
static void project_uv_to_face_2d(const basematrix<double, 2, 3>& map_mat, const dvec3& u_p, const dvec3& v_p,
    basematrix<double, 2, 1>& res1, basematrix<double, 2, 1>& res2)
{
    // map_mat is 2x3 row-major: [u_x, u_y, u_z, v_x, v_y, v_z]
    basematrix<double, 2, 1> res_mat1(map_mat * u_p);
    basematrix<double, 2, 1> res_mat2(map_mat * v_p);
    for (int r = 0; r < 2; ++r) {
        res1[r] = res_mat1[r];
        res2[r] = res_mat2[r];
    }
}

/*
Purpose
    Rotate the face-local tangent basis `(u_f, v_f)` so that its normal aligns with the vertex normal,
    and export the rotated basis as a `2x3` matrix.

High-level Description
    The face has its own tangent basis `(u_f, v_f)` and corresponding normal `face_normal`.
    The vertex has a (possibly different) normal `vertex_normal`, defining the vertex tangent plane.
    When face and vertex normals are not coplanar, curvature quantities expressed in the face frame
    need to be related to the vertex frame. This function computes a rotation that takes the face
    plane into the vertex plane and applies it to the face basis vectors.

    Geometrically:
    - Find the rotation axis as the cross product between `vertex_normal` and `face_normal`.
    - Find the rotation angle from the dot product between the same normals.
    - Build the 3x3 rotation matrix using Rodrigues' formula.
    - Rotate `u_f` and `v_f` by this matrix.
    - Normalize and store the rotated basis in `result` as a 2x3 row-major matrix:
          row 0: rotated `u_f` (x, y, z)
          row 1: rotated `v_f` (x, y, z)

Parameters
    - u_f, v_f
        Input face tangent basis vectors (3D) before rotation; they span the original face plane.
        They are also conceptually the basis of the 2D face parameterization (u, v).

    - vertex_normal
        The normal vector at the vertex. Defines the target plane orientation; after rotation,
        the rotated face basis should be aligned with this normal (vertex tangent plane).

    - face_normal
        The original normal of the face. Together with `vertex_normal`, it determines the rotation
        axis and angle needed to align the two planes.

    - result
        Output `basematrix<double, 2, 3>` that stores the rotated basis vectors in row-major order:
            result[0..2] = rotated_u_f (x, y, z)
            result[3..5] = rotated_v_f (x, y, z)

Notes
    - The rotation axis is `vertex_normal x face_normal`, normalized to unit length.
    - The rotation angle Ë is computed from the arccos of the normalized dot product between
      `vertex_normal` and `face_normal`, with a sign convention given by the code (`-acos(...)`).
    - Rodriguesí rotation formula is used to explicitly build the 3x3 rotation matrix R:
          R = I * cos(Ë) + (1 - cos(Ë)) * (k k^T) + [k]_x * sin(Ë)
      where k is the unit rotation axis.
    - After applying R to `u_f` and `v_f`, the resulting vectors are normalized to ensure they
      remain a clean orthonormal tangent basis in the rotated (vertex-aligned) plane.
    - This rotated basis is later used to express the vertex tangent directions in a face-aligned
      (or face-rotated) 2D coordinate system when converting second fundamental form coefficients
      between frames.
*/
static void rotateElemPlane(dvec3& u_f, dvec3& v_f, dvec3& vertex_normal, dvec3& face_normal, basematrix<double, 2, 3>& result)
{
    // Compute rotation axis
    dvec3 rotation_axis = vertex_normal.cross(face_normal).normalize();
    dvec3 rotated_u_f;
    dvec3 rotated_v_f;

    // compute rotation angle and store sin/cos and rotation matrix components
    double norm_1 = vertex_normal.length();
    double norm_2 = face_normal.length();
    double theta = -acos((vertex_normal.dot(face_normal)) / (norm_1 * norm_2));
    double cth = cos(theta);
    double sth = sin(theta);
    // Rotation matrix components using Rodrigues' rotation formula
    // This is the standard Rodrigues rotation matrix for rotation around unit axis k = rotation_axis by angle Œ∏ :
    // R = I * cos(theta) + (1 - ?cos(theta)) * kk ^ T + [k]_x sin(theta)
    double r11 = cth + rotation_axis.x() * rotation_axis.x() * (1 - cth);
    double r12 = rotation_axis.x() * rotation_axis.y() * (1 - cth) - rotation_axis.z() * sth;
    double r13 = rotation_axis.x() * rotation_axis.z() * (1 - cth) + rotation_axis.y() * sth;
    double r21 = rotation_axis.x() * rotation_axis.y() * (1 - cth) + rotation_axis.z() * sth;
    double r22 = cth + rotation_axis.y() * rotation_axis.y() * (1 - cth);
    double r23 = rotation_axis.z() * rotation_axis.y() * (1 - cth) - rotation_axis.x() * sth;
    double r31 = rotation_axis.z() * rotation_axis.x() * (1 - cth) - rotation_axis.y() * sth;
    double r32 = rotation_axis.z() * rotation_axis.y() * (1 - cth) + rotation_axis.x() * sth;
    double r33 = cth + rotation_axis.z() * rotation_axis.z() * (1 - cth);
    // so R is
    // [r11 r12 r13]
    // [r21 r22 r23]
    // [r31 r32 r33]
    
    // Rotate the basis vectors u_f and v_f using the rotation matrix
    // rotated_u_f = R * u_f
    // rotated_v_f = R * v_f
    rotated_u_f[0] = u_f.x() * r11 + u_f.y() * r12 + u_f.z() * r13;
    rotated_u_f[1] = u_f.x() * r21 + u_f.y() * r22 + u_f.z() * r23;
    rotated_u_f[2] = u_f.x() * r31 + u_f.y() * r32 + u_f.z() * r33;

    rotated_v_f[0] = v_f.x() * r11 + v_f.y() * r12 + v_f.z() * r13;
    rotated_v_f[1] = v_f.x() * r21 + v_f.y() * r22 + v_f.z() * r23;
    rotated_v_f[2] = v_f.x() * r31 + v_f.y() * r32 + v_f.z() * r33;

    rotated_u_f.normalize();
    rotated_v_f.normalize();

    // Store rotated basis in result
    // row 0 : rotated_u_f(x, y, z)
    // row 1 : rotated_v_f(x, y, z)
    result[0] = rotated_u_f.x();
    result[1] = rotated_u_f.y();
    result[2] = rotated_u_f.z();
    result[3] = rotated_v_f.x();
    result[4] = rotated_v_f.y();
    result[5] = rotated_v_f.z();
}

/*
Purpose
    Compute the second fundamental form coefficients [e, f, g] for a single triangular face
    in the face-local tangent frame spanned by (u_f, v_f). These coefficients quantify
    how the surface normal bends across the face and are later converted/accumulated into
    per-vertex curvature data.

High-level behavior
    1. Assemble a compact 12-entry buffer `form_II` via `prepareSecondFundamentalCoefficients`:
       - It projects the triangle vertices into the 2D face-local parameter space (u_f, v_f).
       - It measures how the vertex normals change along the triangle edges.
       - It packs both normal differences and 2D edge vectors into `form_II`.
    2. Pass `form_II` to `leastSquares`, which solves a small linear system (in normal-equation
       form) to fit the second fundamental form coefficients [e, f, g] in a least-squares sense.

Parameters
    - mesh
        Pointer to the owning half-edge mesh. Used indirectly by
        `prepareSecondFundamentalCoefficients` to access per-vertex positions and normals
        for the current face `n_face`.

    - n_face
        The triangle face for which the second fundamental form is being estimated.
        Its vertex indices and centroid are used to build the local 2D parameterization
        and to look up the associated vertex normals.

    - currNode
        The vertex currently being processed at a higher level in the curvature pipeline.
        It is not used inside this function, but is kept for interface symmetry with
        surrounding code and for potential future extensions (e.g., vertex-specific weighting).

    - u_f, v_f
        Orthonormal 3D tangent basis vectors associated with `n_face`. They span the
        face-local tangent plane and define the 2D coordinate system in which [e, f, g]
        are computed (u-direction, v-direction and mixed term).

    - e_f_g
        Output pointer to a 3-element array of doubles, filled as:
            e_f_g[0] = e  (curvature in the u_f direction)
            e_f_g[1] = f  (mixed curvature term between u_f and v_f)
            e_f_g[2] = g  (curvature in the v_f direction)

Internal data layout (form_II)
    `form_II` is a 12-element array populated by `prepareSecondFundamentalCoefficients` as:
        - form_II[0..5]   : normal differences projected into 2D (three edges ? two components)
        - form_II[6..11]  : 2D edge vectors of the triangle (three edges ? two components)
    These entries form the right-hand side (normal variations) and the design matrix
    (edge directions) of the least-squares problem used to recover [e, f, g].

Mathematical background
    The second fundamental form in the face tangent frame is:
        II(u, v) = e du^2 + 2 f du dv + g dv^2
    For each triangle edge, the change in the normal vector is approximately II applied
    to the edge direction in (u, v) coordinates. `leastSquares` solves the overdetermined
    linear system that relates edge directions to observed normal changes and returns
    the best-fit (e, f, g) in the least-squares sense.

Notes
    - This function is called per face, typically from a per-vertex accumulation loop
      (see `calculate_vertex_curvatures`), and its output is later transformed into
      the vertex tangent frame before being area-weighted and summed.
    - Any changes to the layout of `form_II` must be mirrored in `leastSquares` and in
      `prepareSecondFundamentalCoefficients`, as they form a tightly coupled data pipeline.
*/
static void computeSecondFundamentalFromComponents(half_edge_mesh<double>* mesh,
    face<double>& n_face,            // The specific triangle face to analyze
    vertex<double>& currNode,        // Current vertex (note: appears unused!)
    dvec3& u_f, dvec3& v_f,          // Local tangent basis vectors for the face
    double* e_f_g)                // Output: [ee, ff, gg] coefficients
{
    // Create temporary storage (fdifferences[12])
    // This 12-element array will hold intermediate data structured as:
    // Elements [0-5]: Normal differences between triangle vertices
    // Elements [6-11]: Edge coordinate differences in 2D local space
    double fdifferences[12];
    // Evaluate the second fundamental form
    // Projects the triangle into 2D local coordinates using basis(u_f, v_f)
    // Measures how vertex normals change across the triangle
    // Stores edge differences in both normal space and position space
    prepareSecondFundamentalCoefficients(*mesh, n_face, u_f, v_f, fdifferences);
    // Solve for ee,ff,gg using least squares
    // Treats the problem as a linear system : E * [ee, ff, gg] ^ T = N
    // Uses Cholesky decomposition to solve the symmetric normal equations
    // Outputs the three second fundamental form coefficients
    // The second fundamental form coefficients describe :
    // ee: Curvature in the u - direction
    // ff : Twist / mixed curvature term
    // gg : Curvature in the v - direction
    // These are later used to compute :
    // Gaussian curvature(K = ee * gg - ff ^ 2)
    // Mean curvature(H = 0.5 * (ee + gg))
    // Principal curvatures and directions
    leastSquares(fdifferences, e_f_g);
}


/*
convertSecondFundamentalFormComponentsToVertexFrame takes the second fundamental form coefficients
of a face (triangle element) expressed in the face's local tangent basis (u_f, v_f) and re-expresses
them in the vertex tangent basis (u_p, v_p).
- Input e_f_g is [e, f, g] for the triangle face.
- Output e_f_g_V is [e_v, f_v, g_v] for the vertex tangent frame.
This is needed because we:
1.	Compute curvature per face (element).
2.	Want curvature per vertex in a consistent local coordinate system of that vertex.
*/
/*
* How this fits into the overall curvature computation pipeline:
1.	For each face, computeSecondFundamentalFromComponents computes [e, f, g] in that face's local tangent frame.
2.	For each vertex, and each incident face:
-  convertSecondFundamentalFormComponentsToVertexFrame converts the face's [e, f, g] into the vertex's frame [e_v, f_v, g_v].
-  These are then area-weighted and accumulated to get per-vertex second fundamental form.
3.	From the per-vertex second fundamental form, the code later computes:
-  Gaussian & mean curvature.
-  Principal curvatures and directions.
*/
static void convertSecondFundamentalFormComponentsToVertexFrame(
    double* iso_map,          // 2x3 matrix (flattened) mapping 3D -> face 2D (u_f, v_f)
    dvec3& u_p, dvec3& v_p,  // vertex tangent basis (3D)
    dvec3& u_f, dvec3& v_f,  // face tangent basis (3D)
    dvec3& node_normal,      // vertex normal
    dvec3& elem_normal,      // face normal
    double* e_f_g,         // [e, f, g] in face basis
    double* e_f_g_V        // [e_v, f_v, g_v] in vertex basis (output)
)
{
    // Transform vertex basis(u_p, v_p) into face 2D coordinates
    basematrix<double, 2, 1> u_p_in_face_coords(0);
    basematrix<double, 2, 1> v_p_in_face_coords(0);

    // Check if face and vertex planes are coplanar
    if (isCoplanar<double>(u_f, v_f, u_p, v_p)) {
        // Directly express vertex basis in face 2D coordinates
        project_uv_to_face_2d(iso_map, u_p, v_p, u_p_in_face_coords, v_p_in_face_coords);
    }
    else {
        // Rotate face basis into vertex plane
        basematrix<double, 2, 3> iso_map_new(0);
        // Adjust face basis(u_f, v_f) to align with vertex normal
        rotateElemPlane(u_f, v_f, node_normal, elem_normal, iso_map_new);
        // Now express vertex basis in the rotated face 2D coordinates
        project_uv_to_face_2d(iso_map_new, u_p, v_p, u_p_in_face_coords, v_p_in_face_coords);
    }

    // Normalize the resulting 2D vectors
    if (!u_p_in_face_coords[0] && !u_p_in_face_coords[1])
    {
        u_p_in_face_coords[0] = TOLLERANCE<double>;
        u_p_in_face_coords[1] = TOLLERANCE<double>;
    }
    if (!v_p_in_face_coords[0] && !v_p_in_face_coords[1])
    {
        v_p_in_face_coords[0] = TOLLERANCE<double>;
        v_p_in_face_coords[1] = TOLLERANCE<double>;
    }
    dvec3 u_p_f = { u_p_in_face_coords[0], u_p_in_face_coords[1], 0 };
    dvec3 v_p_f = { v_p_in_face_coords[0], v_p_in_face_coords[1], 0 };
    u_p_f.normalize();
    v_p_f.normalize();

    // Re-express second fundamental form coefficients in vertex frame
    double ee = e_f_g[0];
    double ff = e_f_g[1];
    double gg = e_f_g[2];
    /*
    Mathematically, the second fundamental form coefficients transform under a change of basis in the tangent plane.
    - Original form in face basis:
    [ II(u, v) = e, du^2 + 2f, du,dv + g, dv^2 ]
    - When changing to new basis (u_p_f, v_p_f), the new coefficients [e_v, f_v, g_v] are obtained by the standard quadratic form change of coordinates:
    - For direction u_p_f = (a, b) in the face (u, v) coordinates: [ e_v = II(u_p_f, u_p_f) = a^2 e + 2ab f + b^2 g ]
    - For mixed term f_v for basis (u_p_f, v_p_f): [ f_v = II(u_p_f, v_p_f) ] which expands to exactly the combination in the code.
    - For direction v_p_f = (c, d): [ g_v = II(v_p_f, v_p_f) = c^2 e + 2cd f + d^2 g ]
    The code is implementing this basis transformation explicitly:
    - e_f_g_V[0] -> e_v (curvature along vertex u-direction).
    - e_f_g_V[1] -> f_v (mixed curvature in vertex basis).
    - e_f_g_V[2] -> g_v (curvature along vertex v-direction).
    */
    e_f_g_V[0] = u_p_f.x() * u_p_f.x() * ee + 2 * u_p_f.x() * u_p_f.y() * ff + u_p_f.y() * u_p_f.y() * gg;
    e_f_g_V[1] = v_p_f.x() * u_p_f.x() * ee
        + (u_p_f.y() * v_p_f.x() + u_p_f.x() * v_p_f.y()) * ff
        + (u_p_f.y() * v_p_f.y()) * gg;
    e_f_g_V[2] = v_p_f.x() * v_p_f.x() * ee + 2 * v_p_f.x() * v_p_f.y() * ff + v_p_f.y() * v_p_f.y() * gg;
}

/*
Inputs and outputs
- efg is a pointer to 3 floats:
- efg[0] = e (a.k.a. E in some texts)
- efg[1] = f (mixed term)
- efg[2] = g (a.k.a. G in some texts)
These are the vertex-accumulated second fundamental form coefficients in the vertex tangent frame
(they come from the per-element e_f_g_V accumulations we computed earlier).
- principal_curvs[0] = k_min
- principal_curvs[1] = k_max
after the swap step, index 0 is guaranteed to be the smaller value.
- absKmin, absKmax (outputs): magnitudes of the principal curvatures.
- signGauss, signMean (outputs): sign indicators (?1, 0, +1-ish via tolerance) for Gaussian and ‚Äúmean‚Äù curvature.
*/
void transformComponentsToCurvatureData(const double* efg, VertexCurvatureData<double>& curvature_data)
{
    // The Weingarten (shape operator) matrix in this vertex frame is effectively:
    // W = | e  f |
    //     | f  g |
    // -  For a symmetric 2x2 matrix [a b; b c]:
    // -  Trace œ = a + c = k1 + k2 = 2H
    // -  Determinant D = ac - b^2 = K
    // So :
    // -  Gauss_curvature = e * g - f^2 = K
    // -  Mean_curvature = 0.5 * (e + g) = H
    // These are the standard continuous formulas.
    // 
    // Compute Gaussian and Mean curvature from second fundamental form coefficients
    double Gauss_curvature = efg[0] * efg[2] - efg[1] * efg[1];
    curvature_data.gaussCurvature = Gauss_curvature;
    curvature_data.absGaussCurvature = double(fabs(Gauss_curvature));
    // Mean curvature H = 0.5 * (e + g)
    double Mean_curvature = double(0.5 * (efg[0] + efg[2]));
    curvature_data.meanCurvature = Mean_curvature;

    // Solve the characteristic polynomial of the Weingarten matrix to get principal curvatures
    // They solve the characteristic polynomial of the 2x2 matrix :
    // -  For principal curvatures k1, k2 :
    // -  k1 + k2 = 2H
    // -  k1 * k2 = K
    // The eigenvalues k satisfy :
    // -  k^2 - (k1 + k2)k + k1 * k2 = 0
    // -  With k1 + k2 = 2H and k1 * k2 = K :
    // So solve_quadratic(1, -2H, K, ...) returns the two principal curvatures.
    // -  k^2 - 2Hk + k = 0
    std::tuple<double, double> res;
    solve_quadratic<double>(1, -2 * Mean_curvature, Gauss_curvature, res);
    curvature_data.principal_curvatures[0] = std::get<0>(res);
    curvature_data.principal_curvatures[1] = std::get<1>(res);
    if (curvature_data.principal_curvatures[0] > curvature_data.principal_curvatures[1]) {
        std::swap(curvature_data.principal_curvatures[0], curvature_data.principal_curvatures[1]);
    }
    // Determine signs of Gaussian and Mean curvature, and magnitudes of principal curvatures
    double mean = double(0.5 * (curvature_data.principal_curvatures[0] + curvature_data.principal_curvatures[1]));
    curvature_data.signGauss = 0;
    if (mean > TOLLERANCE<double>)
        curvature_data.signMean = 1;
    else if (mean < -TOLLERANCE<double>)
        curvature_data.signMean = -1;
    if (curvature_data.principal_curvatures[0] * curvature_data.principal_curvatures[1] > TOLLERANCE<double>)
        curvature_data.signGauss = 1;
    else if (curvature_data.principal_curvatures[0] * curvature_data.principal_curvatures[1] < -TOLLERANCE<double>)
        curvature_data.signGauss = -1;
    // Store absolute values of principal curvatures
    if (fabs(curvature_data.principal_curvatures[0]) < fabs(curvature_data.principal_curvatures[1])) {
        curvature_data.absKmin = double(fabs(curvature_data.principal_curvatures[0]));
        curvature_data.absKmax = double(fabs(curvature_data.principal_curvatures[1]));
    }
    else {
        curvature_data.absKmax = double(fabs(curvature_data.principal_curvatures[0]));
        curvature_data.absKmin = double(fabs(curvature_data.principal_curvatures[1]));
    }
}

/*
Convert the 2D principal curvature directions (in the local tangent plane basis) into 3D vectors in world space.
-  dmin_2d: 2D direction corresponding to minimum curvature, written as coordinates in the (u_p, v_p) tangent basis.
-  dmax_2d: 2D direction corresponding to maximum curvature, also in (u_p, v_p) coordinates.
-  u_p, v_p: two orthonormal 3D tangent vectors at the vertex (the vertex's local tangent frame).
-  dmin, dmax: resulting 3D principal direction vectors (outputs).

calculatePrincipalDirections computes 2D eigenvectors of the Weingarten matrix W (principal directions in 2D).
principalDirections2Dto3D converts those 2D eigenvectors into 3D using the vertex tangent basis (u_p, v_p).
*/
void principalDirections2Dto3D(double* dmin_2d, double* dmax_2d, dvec3& u_p, dvec3& v_p, dvec3& dmin, dvec3& dmax)
{
    dmin = dvec3(
        dmin_2d[0] * u_p.x() + dmin_2d[1] * v_p.x(),
        dmin_2d[0] * u_p.y() + dmin_2d[1] * v_p.y(),
        dmin_2d[0] * u_p.z() + dmin_2d[1] * v_p.z()
    );
    dmax = dvec3(
        dmax_2d[0] * u_p.x() + dmax_2d[1] * v_p.x(),
        dmax_2d[0] * u_p.y() + dmax_2d[1] * v_p.y(),
        dmax_2d[0] * u_p.z() + dmax_2d[1] * v_p.z()
    );
    dmin.normalize();
    dmax.normalize();
}

/*
It computes principal curvature directions at a vertex:
-  Input W is the 2x2 Weingarten (shape) matrix in the vertex tangent frame:
-  W[0] = e (top-left entry)
-  W[1] = f (off-diagonal entry; matrix is symmetric so bottom-left is also f)
-  W[2] = g (bottom-right entry)
So the matrix is: [ W = [ e  f ][ f  g ] ]
-  Input princ[2] contains the principal curvatures (the eigenvalues of W):
-  princ[0], princ[1]
-  u_p, v_p are orthonormal 3D tangent basis vectors at the vertex (local tangent frame).
-  Output:
-  dmin, dmax are 3D unit vectors giving the minimum and maximum curvature directions at the vertex.
The function:
-  Computes the 2D eigenvectors of W corresponding to the eigenvalues in princ.
-  Converts those 2D vectors into 3D directions using the vertex tangent basis (u_p, v_p) via principalDirections2Dto3D.
*/
void calculatePrincipalDirections(double* W, double* princ, dvec3& u_p, dvec3& v_p, dvec3& dmin, dvec3& dmax)
{
    // Compute eigenvector components with numerical stability checks
    /*
    For a symmetric 2x2 matrix
    [ W = [ a  b ][ b  c ] ]
    and eigenvalue E, an eigenvector [x, y]^T satisfies:
    [ (W - E I)[x  y] = 0 ]
    One row is:
    [ (a - E) x + b y = 0 ]
    They choose to set y = 1, and solve for x:
    [ x = -b/(a - E) ]
    In the code:
    -  a = W[0]
    -  b = W[1]
    -  E is one of princ[0] or princ[1]
    -  denom_1 = a - princ[0], denom_2 = a - princ[1]
    So these are the denominators for the formula x = -b / (a - E).
    */
    double denom_1 = W[0] - princ[0];
    double denom_2 = W[0] - princ[1];

    // Handle near-zero denominators to prevent division by zero
    /*
    If a is very close to E, then a - E is ~0 and the division could explode numerically.
    They clamp denom_1 and denom_2 away from zero by enforcing a minimum magnitude of TOLLERANCE<double>.
    This avoids division by zero / huge eigenvector components.
    */
    if (fabs(denom_1) < TOLLERANCE<double>) denom_1 = TOLLERANCE<double>;
    if (fabs(denom_2) < TOLLERANCE<double>) denom_2 = TOLLERANCE<double>;

    // Compute the 2D eigenvector components (in UV coordinates)
    /*
    For each eigenvalue E:
    -  They define the eigenvector in 2D as [x, 1] where:
    -  x = -b / (a - E) which is implemented as -W[1] / denom.
    So:
    -  For eigenvalue princ[0]: 2D eigenvector ~ [eigvec_1, 1]
    -  For eigenvalue princ[1]: 2D eigenvector ~ [eigvec_2, 1]
    These are still in the 2D tangent coordinate system (u_p, v_p basis).
    */
    double eigvec_1 = -W[1] / denom_1;
    double eigvec_2 = -W[1] / denom_2;

    // Determine which eigenvector corresponds to min/max curvature and convert to 3D
    if (fabs(princ[0]) <= fabs(princ[1]))
    {
        double princ_dir_min[] = { eigvec_1, 1 };
        double princ_dir_max[] = { eigvec_2, 1 };
        principalDirections2Dto3D(princ_dir_min, princ_dir_max, u_p, v_p, dmin, dmax);
    }
    else
    {
        double princ_dir_max[] = { eigvec_1, 1 };
        double princ_dir_min[] = { eigvec_2, 1 };
        principalDirections2Dto3D(princ_dir_min, princ_dir_max, u_p, v_p, dmin, dmax);
    }
}

// calculate_vertex_curvatures computes curvature information at a single vertex v of a triangle mesh:
// -  It aggregates the second fundamental form coefficients from all faces incident to v.
// -  It converts those per-face coefficients into the vertex's tangent frame.
// -  It averages them with area weights.
// -  From that, it later derives:
// -    Principal curvatures (min/max curvature values)
// -    Principal directions (3D directions of min/max curvature)
// -  It finally stores the principal direction vectors on v.
// This function is called from compute_vertex_curvatures for each vertex in the mesh.
//  
// Inputs and Outputs
// Inputs:
// -  half_edge_mesh<double>* mesh
// Full mesh structure. Needed to:
// -  Access incident faces of the vertex (via incident_faces + mesh->faces).
// -  Compute per-face weights via voronoi_based_weighting.
// -  vertex<double>& v
// The specific vertex for which we want curvature:
// -  normal is used to define the vertex tangent frame.
// -  incident_faces is the list of face indices touching this vertex.
// -  id is used by voronoi_based_weighting.
// Outputs:
// -  principal_dir_min - 3D direction of minimum principal curvature at the vertex.
// -  principal_dir_max - 3D direction of maximum principal curvature at the vertex.
static void calculate_vertex_curvatures(half_edge_mesh<double>* mesh, vertex<double>& v) {
    // -  sum_ee, sum_ff, sum_gg -> running area - weighted sums of second fundamental form components in the vertex frame :
    // -  ee -> curvature in local u direction
    // -  ff -> mixed term
    // -  gg -> curvature in local v direction
    // -  sum_of_V_areas -> sum of per - face vertex areas(weights) for normalization.
    double sum_ee = 0;
    double sum_ff = 0;
    double sum_gg = 0;
    double sum_of_V_areas = 0;

    // -  Given the vertex normal normal, create_uv_reference_plane constructs an orthonormal tangent basis:
    // -  u_vertex, v_vertex are tangent to the surface at the vertex and orthogonal to normal.
    // -  This (u_vertex, v_vertex) is the vertex tangent frame in which we ultimately want curvatures and directions.
    dvec3 u_vertex, v_vertex;
    create_uv_reference_plane(v.normal, u_vertex, v_vertex);

    // Iterate over all faces incident to the vertex
    for (auto& i : v.incident_faces) {
        // Get the face and its normal
        face<double>* n_face = mesh->faces[i];
        dvec3& face_normal = n_face->normal;

        // Create local tangent basis for the face
        // (u_face, v_face) is a 3D tangent basis aligned with the face
        dvec3 u_face, v_face;
        create_uv_reference_plane(face_normal, u_face, v_face);
        // -  ref_map encodes a 2x3 matrix mapping 3D vectors into the face's 2D (u_f, v_f) coordinate system:
        // -  Row 0: components of u_face
        // -  Row 1: components of v_face
        // So ref_map * [x y z]^T gives 2D coordinates of a 3D vector in the face's tangent plane.
        double ref_map[] = { u_face.x(), u_face.y(), u_face.z(), v_face.x(), v_face.y(), v_face.z() };

        // Compute per face second fundamental form (face frame)
        // -  This function:
        // -    Projects the triangle into the face's 2D coordinates.
        // -    Looks at normal differences between vertices on that face.
        // -    Uses a least squares fit to get [e, f, g] for that face, in the face frame (u_face, v_face).
        // After this call:
        // -  e_f_g[0]: curvature along u_face.
        // -  e_f_g[1]: mixed curvature term in face frame.
        // -  e_f_g[2]: curvature along v_face.
        // Note: they are still in the face tangent frame, not yet in the vertex tangent frame.
        double e_f_g[3];
        computeSecondFundamentalFromComponents(mesh, *n_face, v, u_face, v_face, e_f_g);

        // Transform face-frame [e, f, g] into vertex-frame [e_v, f_v, g_v]
        // -  convertSecondFundamentalFormComponentsToVertexFrame:
        // -    Expresses the vertex tangent basis u_vertex, v_vertex in the face's coordinates (or a rotated version if face and vertex planes are not coplanar).
        // -    Applies the quadratic form basis change to re-express [e, f, g] in the vertex frame.
        // -  On output:
        // -    e_f_g_v[0] = e_v (curvature along u_vertex)
        // -    e_f_g_v[1] = f_v (mixed term in vertex frame)
        // -    e_f_g_v[2] = g_v (curvature along v_vertex)
        // So these are per-face contributions already expressed in the vertex's tangent frame.
        double e_f_g_v[3];
        convertSecondFundamentalFormComponentsToVertexFrame(ref_map, u_vertex, v_vertex, u_face, v_face, v.normal, face_normal, e_f_g, e_f_g_v);

        // Area weighting and accumulation

        // -  v_area is an area-like weight for the pair (vertex v, face n_face):
        // -    Typically a Voronoi area or similar local area assigned to v from that face.
        // -    Encodes how much of the surface "near" v is contributed by this face.
        double v_area = voronoi_based_weighting<double>(*mesh, v.id, *n_face);

        // -  Each coefficient is summed with this area weight.
        // -  Conceptually: you're computing an area-weighted average of the second fundamental form over the 1-ring of faces around v.        sum_ee += v_area * e_f_g_v[0];
        sum_ff += v_area * e_f_g_v[1];
        sum_gg += v_area * e_f_g_v[2];

        sum_of_V_areas += v_area;
    }
    // Prevent division by zero in normalization
    if (!sum_of_V_areas)
        sum_of_V_areas = 1;

    // Normalize to get averaged second fundamental form at vertex
    // If a sum is exactly zero, it's replaced by a small TOLLERANCE_SQ<double> 
    // to avoid degenerate zero matrices in later computations.
    double normalized_ee = (!sum_ee) ? TOLLERANCE_SQ<double> : (sum_ee / sum_of_V_areas);
    double normalized_ff = (!sum_ff) ? TOLLERANCE_SQ<double> : (sum_ff / sum_of_V_areas);
    double normalized_gg = (!sum_gg) ? TOLLERANCE_SQ<double> : (sum_gg / sum_of_V_areas);

    // Transform second fundamental form components into curvature data
    // W represents the entries of the 2x2 Weingarten (shape) matrix in the vertex tangent frame:
    // W = [ e  f ]
    //     [ f  g ]
    // 
    // where:
    //   e = W[0]
    //   f = W[1]
    //   g = W[2]
    double W[] = { normalized_ee,normalized_ff,normalized_gg };

    // Derive principal curvature values and curvature signs
    VertexCurvatureData<double>& curvature_data = v.curvature_data;

    // -  transformComponentsToCurvatureData:
    // -    Computes Gaussian curvature K = e*g - f^2.
    // -    Computes Mean curvature H = 0.5*(e + g).
    // -    Solves the eigenvalue problem of W to get the 2 principal curvatures (principal_curvs[0], principal_curvs[1]).
    // -    Determines sign information and absolute magnitudes.
    // At this point, we have the eigenvalues (principal curvature magnitudes) at the vertex.
    transformComponentsToCurvatureData(W, curvature_data);

    //Get position of Node
    //dvec3 Node_3d = v.coords;

    // Compute principal directions in 3D
    // -  calculatePrincipalDirections:
    // -  For the 2x2 symmetric matrix W, uses each eigenvalue k1, k2 to compute its 2D eigenvector in (u_vertex, v_vertex) coordinates.
    // -  Chooses which eigenvector corresponds to min and max principal curvature.
    // -  Transforms those 2D directions into 3D using u_vertex and v_vertex, and normalizes them.
    // So:
    // -  curvature_data.principal_directions[0]: 3D direction of minimum curvature.
    // -  curvature_data.principal_directions[1]: 3D direction of maximum curvature.
    calculatePrincipalDirections(W, curvature_data.principal_curvatures, u_vertex, v_vertex,
        curvature_data.principal_directions[0], curvature_data.principal_directions[1]);
}

// compute_vertex_curvatures computes curvature information for all vertices in the mesh.
void compute_vertex_curvatures(half_edge_mesh<double>* mesh) {
    // Precompute face properties and vertex normals
    mesh->compute_face_properties();
    mesh->compute_vertex_normals();

    // Compute curvature for each vertex in the mesh
    for (auto& v : mesh->vertices) {
        calculate_vertex_curvatures(mesh, *(v.second));
    }
    mesh->curvatures_computed() = true;
}
