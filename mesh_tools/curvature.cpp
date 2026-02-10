#include "base_definitions.h"

#include "vector.h"

#include "mesh.h"
#include "mesh_curvature.h"

using namespace base_math;

/*
Purpose
    WorldCoordsToFaceCoords converts the 3D coordinates of a triangle's vertices (in world space) into 2D coordinates in the face's local tangent plane, defined by basis vectors u_f and v_f.
    This is a standard step in curvature computation: you want to work in a 2D parameterization of the triangle (u,v) instead of raw 3D coordinates.
Parameters
    -const fvec3& u_f, const fvec3& v_f
    Orthonormal tangent basis for the face. Together with the face normal, they form a local 3D frame, but only the tangent u_f, v_f are used here.
    -const fvec3* vertices
    Pointer to an array of three 3D vertex positions of the face: vertices[0], vertices[1], vertices[2].
    -const fvec3& centroid
    The 3D centroid of the face. Used as origin of the local coordinate system.
    -fvec2* results
Output array of three 2D coordinates (u,v) for each vertex in the face local frame:
    -results[0] -> 2D coords of vertices[0]
    -results[1] -> 2D coords of vertices[1]
    -results[2] -> 2D coords of vertices[2]
*/
static void WorldCoordsToFaceCoords(const fvec3& u_f, const fvec3& v_f, const fvec3* vertices, const fvec3& centroid, fvec2* results)
{
    // Project vertices onto plane defined by centroid
    fvec3 proj_n1 = vertices[0] - centroid;
    fvec3 proj_n2 = vertices[1] - centroid;
    fvec3 proj_n3 = vertices[2] - centroid;
    // Now proj_n1, proj_n2, proj_n3 are vectors from the centroid to each vertex.
    // So the centroid effectively becomes the(0, 0) origin in the local coordinates.

    // u,v in a flat array for better performance. Conceptually this is a 2x3 matrix:
    float uvf[] = { u_f.x(), u_f.y(), u_f.z(), v_f.x(), v_f.y(), v_f.z() };

    // Unrolled matrix-vector multiplication for better performance on fixed 2x3 operation
    results[0][0] = uvf[0] * proj_n1.x() + uvf[1] * proj_n1.y() + uvf[2] * proj_n1.z();
    results[0][1] = uvf[3] * proj_n1.x() + uvf[4] * proj_n1.y() + uvf[5] * proj_n1.z();

    results[1][0] = uvf[0] * proj_n2.x() + uvf[1] * proj_n2.y() + uvf[2] * proj_n2.z();
    results[1][1] = uvf[3] * proj_n2.x() + uvf[4] * proj_n2.y() + uvf[5] * proj_n2.z();

    results[2][0] = uvf[0] * proj_n3.x() + uvf[1] * proj_n3.y() + uvf[2] * proj_n3.z();
    results[2][1] = uvf[3] * proj_n3.x() + uvf[4] * proj_n3.y() + uvf[5] * proj_n3.z();
}

// This function prepares data for computing the second fundamental form of a triangular mesh face.
// It transforms 3D normal differences into a 2D coordinate system and packages them with 2D position
// differences for later least - squares fitting.
static void convert3dTo2dCoords(const fmat2x3& uv_map, const fvec2(&nodes_local)[3], const fvec3(&normals)[3], float* result)
{
    // Calculates how vertex normals change across the triangle's three edges. This captures how the surface is bending.
    fvec3 diff_n2n1(normals[2] - normals[1]);
    fvec3 diff_n0n2(normals[0] - normals[2]);
    fvec3 diff_n1n0(normals[1] - normals[0]);

    // Project these 3D normal differences into the 2D local coordinate system defined by uv_map
    fvec3 diff_n2_n1_3d(diff_n2n1.x(), diff_n2n1.y(), diff_n2n1.z());
    fvec3 diff_n0_n2_3d(diff_n0n2.x(), diff_n0n2.y(), diff_n0n2.z());
    fvec3 diff_n1_n0_3d(diff_n1n0.x(), diff_n1n0.y(), diff_n1n0.z());

    fvec2 diff_n2_n1_2d(0);
    fvec2 diff_n0_n2_2d(0);
    fvec2 diff_n1_n0_2d(0);

    // Unrolled matrix-vector multiplication for 2x3 * 3x1 operations
    // uv_map is 2x3 row-major: [u_x, u_y, u_z, v_x, v_y, v_z]
    diff_n2_n1_2d[0] = uv_map[0] * diff_n2_n1_3d.x() + uv_map[1] * diff_n2_n1_3d.y() + uv_map[2] * diff_n2_n1_3d.z();
    diff_n2_n1_2d[1] = uv_map[3] * diff_n2_n1_3d.x() + uv_map[4] * diff_n2_n1_3d.y() + uv_map[5] * diff_n2_n1_3d.z();

    diff_n0_n2_2d[0] = uv_map[0] * diff_n0_n2_3d.x() + uv_map[1] * diff_n0_n2_3d.y() + uv_map[2] * diff_n0_n2_3d.z();
    diff_n0_n2_2d[1] = uv_map[3] * diff_n0_n2_3d.x() + uv_map[4] * diff_n0_n2_3d.y() + uv_map[5] * diff_n0_n2_3d.z();

    diff_n1_n0_2d[0] = uv_map[0] * diff_n1_n0_3d.x() + uv_map[1] * diff_n1_n0_3d.y() + uv_map[2] * diff_n1_n0_3d.z();
    diff_n1_n0_2d[1] = uv_map[3] * diff_n1_n0_3d.x() + uv_map[4] * diff_n1_n0_3d.y() + uv_map[5] * diff_n1_n0_3d.z();

    // Also compute edge differences in local 2D coordinates
    fvec2 diff_nodes_n2n1 = nodes_local[2] - nodes_local[1];
    fvec2 diff_nodes_n0n2 = nodes_local[0] - nodes_local[2];
    fvec2 diff_nodes_n1n0 = nodes_local[1] - nodes_local[0];

    // Store results in the output array
    // Layout:
    // [0-5]: Normal difference (6 floats: 3 edges x 2 components)
    // [6-11]: Edge coordinate differences (6 floats: 3 edges x 2 components)
    // The output feeds into leastSquares(), which solves
    // E * [ee, ff, gg]^T = N
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

static void prepareSeconFundamentalCoefficients(half_edge_mesh<float>& mesh, const face<float>& cur_face, const fvec3& u_f, const fvec3& v_f, float* results)
{
    fvec3 elemCentroid = cur_face.center;
    fvec3 world_coords[3] = { mesh.vertices[cur_face.v1]->coords, mesh.vertices[cur_face.v2]->coords,mesh.vertices[cur_face.v3]->coords };
    fvec2 local_coords[3];
    WorldCoordsToFaceCoords(u_f, v_f, world_coords, elemCentroid, local_coords);
    const fvec3& node_normal_1 = mesh.vertices[cur_face.v1]->normal;
    const fvec3& node_normal_2 = mesh.vertices[cur_face.v2]->normal;
    const fvec3& node_normal_3 = mesh.vertices[cur_face.v3]->normal;
    fmat2x3 uv_map({ u_f.x(), u_f.y(), u_f.z(), v_f.x(), v_f.y(), v_f.z() });
    convert3dTo2dCoords(uv_map, local_coords, { node_normal_1, node_normal_2, node_normal_3 }, results);
}

/*
 * The leastSquares helper extracts the per-triangle second fundamental form components (ee, ff, gg)
 * from a previously assembled 12-float array
**/
static void leastSquares(const float* second_fundamental, float* ee_ff_gg)
{
    // N is the 6x1 RHS vector storing measured normal differences second_fundamental[0..5]
    basematrix<float, 6, 1> N({ second_fundamental[0],second_fundamental[1],second_fundamental[2],
                    second_fundamental[3],second_fundamental[4],second_fundamental[5] });
    // E is the 6x3 matrix storing the local edge differences second_fundamental[6..17]
    // arranged as rows: [du1 dv1 0 0 du1 dv1; du2 dv2 0 0 du2 dv2; du3 dv3 0 0 du3 dv3]
    basematrix<float, 6, 3> E({ second_fundamental[6], second_fundamental[7],0,0,second_fundamental[6],second_fundamental[7] ,
                   second_fundamental[8], second_fundamental[9],0,0,second_fundamental[8],second_fundamental[9] ,
                  second_fundamental[10], second_fundamental[11],0,0,second_fundamental[10],second_fundamental[11] });

    // Transpose of E (3x6)
    basematrix<float, 3, 6> Trans_E = E.transpose();
    // Compute normal equations components
    basematrix<float, 3, 3> tr_x_e = Trans_E * E;           // yields the symmetric 3x3 normal matrix -> E^T * E (3x3)
    basematrix<float, 3, 1> rhs = Trans_E * N;              // accumulates  E^T * N (3x1)
    basematrix<float, 3, 1> sol;                            // will hold [ee ff gg]^T

    // Solving tr_x_e * sol = rhs recovers[ee, ff, gg]^T.
    // Cholesky is used because E^T*E is symmetric positive semi-definite.
    if (!cholesky_solve_3x3<float>(tr_x_e, rhs, sol)) {
        // If the first solve fails(e.g., because E^T*E is near-singular), 
        // the diagonal gets a tiny TOLLERANCE<float> bump, 
        // effectively applying Tikhonov regularization, and the solve is retried.
        basematrix<float, 3, 3> reg = tr_x_e;
        for (int i = 0; i < 3; ++i) reg[i * 3 + i] += TOLLERANCE<float>;
        cholesky_solve_3x3<float>(reg, rhs, sol);
    }

    ee_ff_gg[0] = sol[0];
    ee_ff_gg[1] = sol[1];
    ee_ff_gg[2] = sol[2];
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
-   Map: A 2x3 transformation matrix (stored as a flat float array of 6 elements)
-   u_p, v_p: 3D basis vectors in vertex space
-   res1, res2: Output 2x1 vectors (2 floats each) representing u_p and v_p in 2D element space

The Math
The function performs two matrix-vector multiplications:
res1 = map_mat x u_p    (2x3 matrix times 3x1 vector = 2x1 result)
res2 = map_mat x v_p
*/
static void project_uv_to_face_2d(const basematrix<float, 2, 3>& map_mat, const fvec3& u_p, const fvec3& v_p,
    basematrix<float, 2, 1>& res1, basematrix<float, 2, 1>& res2)
{
    // map_mat is 2x3 row-major: [u_x, u_y, u_z, v_x, v_y, v_z]
    basematrix<float, 2, 1> res_mat1(map_mat * u_p);
    basematrix<float, 2, 1> res_mat2(map_mat * v_p);
    for (int r = 0; r < 2; ++r) {
        res1[r] = res_mat1[r];
        res2[r] = res_mat2[r];
    }
}

/*
High-level purpose
rotateElemPlane takes a local tangent basis of a face (u_f, v_f) and rotates it so that the face’s normal
aligns with the vertex normal (vertex_normal).
This is used when the face plane and vertex plane are not coplanar: you want to “tilt” the face tangent
basis into the vertex’s tangent plane so that curvature quantities can be consistently expressed in the vertex frame.
*/
static void rotateElemPlane(fvec3& u_f, fvec3& v_f, fvec3& vertex_normal, fvec3& face_normal, basematrix<float, 2, 3>& result)
{
    // Compute rotation axis
    fvec3 rotation_axis = (vertex_normal * face_normal).normalize();
    fvec3 rotated_u_f;
    fvec3 rotated_v_f;

    // compute rotation angle and store sin/cos and rotation matrix components
    float norm_1 = vertex_normal.length();
    float norm_2 = face_normal.length();
    float theta = -acos((vertex_normal.dot(face_normal)) / (norm_1 * norm_2));
    float cth = cos(theta);
    float sth = sin(theta);
    // Rotation matrix components using Rodrigues' rotation formula
    // This is the standard Rodrigues rotation matrix for rotation around unit axis k = rotation_axis by angle è :
    // R = I * cosè + (1?cosè) * k, k ^ T + [k]_x sinè
    float r11 = cth + rotation_axis.x() * rotation_axis.x() * (1 - cth);
    float r12 = rotation_axis.x() * rotation_axis.y() * (1 - cth) - rotation_axis.z() * sth;
    float r13 = rotation_axis.x() * rotation_axis.z() * (1 - cth) + rotation_axis.y() * sth;
    float r21 = rotation_axis.x() * rotation_axis.y() * (1 - cth) + rotation_axis.z() * sth;
    float r22 = cth + rotation_axis.y() * rotation_axis.y() * (1 - cth);
    float r23 = rotation_axis.z() * rotation_axis.y() * (1 - cth) - rotation_axis.x() * sth;
    float r31 = rotation_axis.z() * rotation_axis.x() * (1 - cth) - rotation_axis.y() * sth;
    float r32 = rotation_axis.z() * rotation_axis.y() * (1 - cth) + rotation_axis.x() * sth;
    float r33 = cth + rotation_axis.z() * rotation_axis.z() * (1 - cth);

    // Rotate u_f and v_f using the rotation matrix
    rotated_u_f[0] = u_f.x() * r11 + u_f.y() * r12 + u_f.z() * r13;
    rotated_u_f[1] = u_f.x() * r21 + u_f.y() * r22 + u_f.z() * r23;
    rotated_u_f[2] = u_f.x() * r31 + u_f.y() * r32 + u_f.z() * r33;
    rotated_v_f[0] = v_f.x() * r11 + v_f.y() * r12 + v_f.z() * r13;
    rotated_v_f[1] = v_f.x() * r21 + v_f.y() * r22 + v_f.z() * r23;
    rotated_v_f[2] = v_f.x() * r31 + v_f.y() * r32 + v_f.z() * r33;
    rotated_u_f.normalize();
    rotated_v_f.normalize();

    // Store rotated basis in result
    result[0] = rotated_u_f.x();
    result[1] = rotated_u_f.y();
    result[2] = rotated_u_f.z();
    result[3] = rotated_v_f.x();
    result[4] = rotated_v_f.y();
    result[5] = rotated_v_f.z();
}

static void computeSecondFundamentalFormComponentsPerElem(half_edge_mesh<float>* mesh,
    face<float>& n_face,            // The specific triangle face to analyze
    vertex<float>& currNode,        // Current vertex (note: appears unused!)
    fvec3& u_f, fvec3& v_f,         // Local tangent basis vectors for the face
    float* ee_ff_gg)                // Output: [ee, ff, gg] coefficients
{
    // Create temporary storage (form_II[12])
    // This 12-element array will hold intermediate data structured as:
    // Elements [0-5]: Normal differences between triangle vertices
    // Elements [6-11]: Edge coordinate differences in 2D local space
    float form_II[12];
    // Evaluate the second fundamental form
    // Projects the triangle into 2D local coordinates using basis(u_f, v_f)
    // Measures how vertex normals change across the triangle
    // Stores edge differences in both normal space and position space
    prepareSeconFundamentalCoefficients(*mesh, n_face, u_f, v_f, form_II);
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
    leastSquares(form_II, ee_ff_gg);
}


/*
convertSecondFundamentalFormComponentsToVertexFrame takes the second fundamental form coefficients
of a face (triangle element) expressed in the face’s local tangent basis (u_f, v_f) and re-expresses 
them in the vertex tangent basis (u_p, v_p).
- Input ee_ff_gg is [e, f, g] for the triangle face.
- Output ee_ff_gg_V is [e_v, f_v, g_v] for the vertex tangent frame.
This is needed because we:
1.	Compute curvature per face (element).
2.	Want curvature per vertex in a consistent local coordinate system of that vertex.
*/
/*
* How this fits into the overall curvature computation pipeline:
1.	For each face, computeSecondFundamentalFormComponentsPerElem computes [e, f, g] in that face's local tangent frame.
2.	For each vertex, and each incident face:
-  convertSecondFundamentalFormComponentsToVertexFrame converts the face's [e, f, g] into the vertex's frame [e_v, f_v, g_v].
-  These are then area-weighted and accumulated to get per-vertex second fundamental form.
3.	From the per-vertex second fundamental form, the code later computes:
-  Gaussian & mean curvature.
-  Principal curvatures and directions.
*/
static void convertSecondFundamentalFormComponentsToVertexFrame(
    float* iso_map,          // 2x3 matrix (flattened) mapping 3D -> face 2D (u_f, v_f)
    fvec3& u_p, fvec3& v_p,  // vertex tangent basis (3D)
    fvec3& u_f, fvec3& v_f,  // face tangent basis (3D)
    fvec3& node_normal,      // vertex normal
    fvec3& elem_normal,      // face normal
    float* ee_ff_gg,         // [e, f, g] in face basis
    float* ee_ff_gg_V        // [e_v, f_v, g_v] in vertex basis (output)
)
{
    // Transform vertex basis(u_p, v_p) into face 2D coordinates
    basematrix<float, 2, 1> u_p_in_face_coords(0);
    basematrix<float, 2, 1> v_p_in_face_coords(0);

    // Check if face and vertex planes are coplanar
    if (isCoplanar<float>(u_f, v_f, u_p, v_p)) {
        // Directly express vertex basis in face 2D coordinates
        project_uv_to_face_2d(iso_map, u_p, v_p, u_p_in_face_coords, v_p_in_face_coords);
    }
    else {
        // Rotate face basis into vertex plane
        basematrix<float, 2, 3> iso_map_new(0);
        // Adjust face basis(u_f, v_f) to align with vertex normal
        rotateElemPlane(u_f, v_f, node_normal, elem_normal, iso_map_new);
        // Now express vertex basis in the rotated face 2D coordinates
        project_uv_to_face_2d(iso_map_new, u_p, v_p, u_p_in_face_coords, v_p_in_face_coords);
    }

    // Normalize the resulting 2D vectors
    if (!u_p_in_face_coords[0] && !u_p_in_face_coords[1])
    {
        u_p_in_face_coords[0] = TOLLERANCE<float>;
        u_p_in_face_coords[1] = TOLLERANCE<float>;
    }
    if (!v_p_in_face_coords[0] && !v_p_in_face_coords[1])
    {
        v_p_in_face_coords[0] = TOLLERANCE<float>;
        v_p_in_face_coords[1] = TOLLERANCE<float>;
    }
    fvec3 u_p_f = { u_p_in_face_coords[0], u_p_in_face_coords[1], 0 };
    fvec3 v_p_f = { v_p_in_face_coords[0], v_p_in_face_coords[1], 0 };
    u_p_f.normalize();
    v_p_f.normalize();

    // Re-express second fundamental form coefficients in vertex frame
    float ee = ee_ff_gg[0];
    float ff = ee_ff_gg[1];
    float gg = ee_ff_gg[2];
    /*
    Mathematically, the second fundamental form coefficients transform under a change of basis in the tangent plane.
    - Original form in face basis:
    [ II(u, v) = e, du^2 + 2f, du,dv + g, dv^2 ]
    - When changing to new basis (u_p_f, v_p_f), the new coefficients [e_v, f_v, g_v] are obtained by the standard quadratic form change of coordinates:
    - For direction u_p_f = (a, b) in the face (u, v) coordinates: [ e_v = II(u_p_f, u_p_f) = a^2 e + 2ab f + b^2 g ]
    - For mixed term f_v for basis (u_p_f, v_p_f): [ f_v = II(u_p_f, v_p_f) ] which expands to exactly the combination in the code.
    - For direction v_p_f = (c, d): [ g_v = II(v_p_f, v_p_f) = c^2 e + 2cd f + d^2 g ]
    The code is implementing this basis transformation explicitly:
    - ee_ff_gg_V[0] -> e_v (curvature along vertex u-direction).
    - ee_ff_gg_V[1] -> f_v (mixed curvature in vertex basis).
    - ee_ff_gg_V[2] -> g_v (curvature along vertex v-direction).
    */
    ee_ff_gg_V[0] = u_p_f.x() * u_p_f.x() * ee + 2 * u_p_f.x() * u_p_f.y() * ff + u_p_f.y() * u_p_f.y() * gg;
    ee_ff_gg_V[1] = v_p_f.x() * u_p_f.x() * ee
        + (u_p_f.y() * v_p_f.x() + u_p_f.x() * v_p_f.y()) * ff
        + (u_p_f.y() * v_p_f.y()) * gg;
    ee_ff_gg_V[2] = v_p_f.x() * v_p_f.x() * ee + 2 * v_p_f.x() * v_p_f.y() * ff + v_p_f.y() * v_p_f.y() * gg;
}

/*
Inputs and outputs
•	accum_II is a pointer to 3 floats:
•	accum_II[0] = e (a.k.a. E in some texts)
•	accum_II[1] = f (mixed term)
•	accum_II[2] = g (a.k.a. G in some texts)
These are the vertex-accumulated second fundamental form coefficients in the vertex tangent frame
(they come from the per-element ee_ff_gg_V accumulations we computed earlier).
•	principal_curvs[0] = k_min
•	principal_curvs[1] = k_max
after the swap step, index 0 is guaranteed to be the smaller value.
•	absKmin, absKmax (outputs): magnitudes of the principal curvatures.
•	signGauss, signMean (outputs): sign indicators (?1, 0, +1-ish via tolerance) for Gaussian and “mean” curvature.
*/
void transformComponentsToCurvatureData(const float* accum_II, VertexCurvatureData<float>& curvature_data)
    // float* principal_curvs,
    // float& absKmin, float& absKmax,
    // int& signGauss, int& signMean)
{
    // -  e = accum_II[0]
    // -  f = accum_II[1]
    // -  g = accum_II[2]
    // The Weingarten (shape operator) matrix in this vertex frame is effectively:
    // W = | e  f |
    //     | f  g |
    // -  For a symmetric 2?2 matrix [a b; b c]:
    // -  Trace ô = a + c = k1 + k2 = 2H
    // -  Determinant Ä = ac - b^2 = K
    // So :
    // -  Gauss_curvature = e * g - f^2 = K
    // -  Mean_curvature = 0.5 * (e + g) = H
    // These are the standard continuous formulas.
    // 
    // Compute Gaussian and Mean curvature from second fundamental form coefficients
    float Gauss_curvature = accum_II[0] * accum_II[2] - accum_II[1] * accum_II[1];
    curvature_data.gaussCurvature = Gauss_curvature;
    curvature_data.absGaussCurvature = float(fabs(Gauss_curvature));
    // Mean curvature H = 0.5 * (e + g)
    float Mean_curvature = float(0.5 * (accum_II[0] + accum_II[2]));
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
    std::tuple<float, float> res;
    solve_quadratic<float>(1, -2 * Mean_curvature, Gauss_curvature, res);
    curvature_data.principal_curvatures[0] = std::get<0>(res);
    curvature_data.principal_curvatures[1] = std::get<1>(res);
    if (curvature_data.principal_curvatures[0] > curvature_data.principal_curvatures[1]) {
        std::swap(curvature_data.principal_curvatures[0], curvature_data.principal_curvatures[1]);
    }
    // Determine signs of Gaussian and Mean curvature, and magnitudes of principal curvatures
    float mean = float(0.5 * (curvature_data.principal_curvatures[0] + curvature_data.principal_curvatures[1]));
    curvature_data.signGauss = 0;
    if (mean > TOLLERANCE<float>)
        curvature_data.signMean = 1;
    else if (mean < -TOLLERANCE<float>)
        curvature_data.signMean = -1;
    if (curvature_data.principal_curvatures[0] * curvature_data.principal_curvatures[1] > TOLLERANCE<float>)
        curvature_data.signGauss = 1;
    else if (curvature_data.principal_curvatures[0] * curvature_data.principal_curvatures[1] < -TOLLERANCE<float>)
        curvature_data.signGauss = -1;
    // Store absolute values of principal curvatures
    if (fabs(curvature_data.principal_curvatures[0]) < fabs(curvature_data.principal_curvatures[1])) {
        curvature_data.absKmin = float(fabs(curvature_data.principal_curvatures[0]));
        curvature_data.absKmax = float(fabs(curvature_data.principal_curvatures[1]));
    }
    else {
        curvature_data.absKmax = float(fabs(curvature_data.principal_curvatures[0]));
        curvature_data.absKmin = float(fabs(curvature_data.principal_curvatures[1]));
    }
}

/*
Convert the 2D principal curvature directions (in the local tangent plane basis) into 3D vectors in world space.
-  dmin_2d: 2D direction corresponding to minimum curvature, written as coordinates in the (u_p, v_p) tangent basis.
-  dmax_2d: 2D direction corresponding to maximum curvature, also in (u_p, v_p) coordinates.
-  u_p, v_p: two orthonormal 3D tangent vectors at the vertex (the vertex’s local tangent frame).
-  dmin, dmax: resulting 3D principal direction vectors (outputs).

calculatePrincipalDirections computes 2D eigenvectors of the Weingarten matrix W (principal directions in 2D).
principalDirections2Dto3D converts those 2D eigenvectors into 3D using the vertex tangent basis (u_p, v_p).
*/
void principalDirections2Dto3D(float* dmin_2d, float* dmax_2d, fvec3& u_p, fvec3& v_p, fvec3& dmin, fvec3& dmax)
{
    dmin = fvec3(
        dmin_2d[0] * u_p.x() + dmin_2d[1] * v_p.x(),
        dmin_2d[0] * u_p.y() + dmin_2d[1] * v_p.y(),
        dmin_2d[0] * u_p.z() + dmin_2d[1] * v_p.z()
    );
    dmax = fvec3(
        dmax_2d[0] * u_p.x() + dmax_2d[1] * v_p.x(),
        dmax_2d[0] * u_p.y() + dmax_2d[1] * v_p.y(),
        dmax_2d[0] * u_p.z() + dmax_2d[1] * v_p.z()
    );
    dmin.normalize();
    dmax.normalize();
}

/*
It computes principal curvature directions at a vertex:
-  Input W is the 2?2 Weingarten (shape) matrix in the vertex tangent frame:
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
void calculatePrincipalDirections(float* W, float* princ, fvec3& u_p, fvec3& v_p, fvec3& dmin, fvec3& dmax)
{
    // Compute eigenvector components with numerical stability checks
    /*
    For a symmetric 2x2 matrix
    [ W = [ a  b ][ b  c ] ]
    and eigenvalue ë, an eigenvector [x, y]^T satisfies:
    [ (W - ë I)[x  y] = 0 ]
    One row is:
    [ (a - ë) x + b y = 0 ]
    They choose to set y = 1, and solve for x:
    [ x = -b/(a - ë) ]
    In the code:
    -  a = W[0]
    -  b = W[1]
    -  ë is one of princ[0] or princ[1]
    -  denom_1 = a - princ[0], denom_2 = a - princ[1]
    So these are the denominators for the formula x = -b / (a - ë).
    */
    float denom_1 = W[0] - princ[0];
    float denom_2 = W[0] - princ[1];
    
    // Handle near-zero denominators to prevent division by zero
    /*
    If a is very close to ë, then a - ë is ~0 and the division could explode numerically.
    They clamp denom_1 and denom_2 away from zero by enforcing a minimum magnitude of TOLLERANCE<float>.
    This avoids division by zero / huge eigenvector components.
    */
    if (fabs(denom_1) < TOLLERANCE<float>) denom_1 = TOLLERANCE<float>;
    if (fabs(denom_2) < TOLLERANCE<float>) denom_2 = TOLLERANCE<float>;
    
    // Compute the 2D eigenvector components (in UV coordinates)
    /*
    For each eigenvalue ë:
    -  They define the eigenvector in 2D as [x, 1] where:
    -  x = -b / (a - ë) which is implemented as -W[1] / denom.
    So:
    -  For eigenvalue princ[0]: 2D eigenvector ~ [eigvec_1, 1]
    -  For eigenvalue princ[1]: 2D eigenvector ~ [eigvec_2, 1]
    These are still in the 2D tangent coordinate system (u_p, v_p basis).
    */
    float eigvec_1 = -W[1] / denom_1;
    float eigvec_2 = -W[1] / denom_2;
    
    // Determine which eigenvector corresponds to min/max curvature and convert to 3D
    if (fabs(princ[0]) <= fabs(princ[1]))
    {
        float princ_dir_min[] = { eigvec_1, 1 };
        float princ_dir_max[] = { eigvec_2, 1 };
        principalDirections2Dto3D(princ_dir_min, princ_dir_max, u_p, v_p, dmin, dmax);
    }
    else
    {
        float princ_dir_max[] = { eigvec_1, 1 };
        float princ_dir_min[] = { eigvec_2, 1 };
        principalDirections2Dto3D(princ_dir_min, princ_dir_max, u_p, v_p, dmin, dmax);
    }
}

// calculate_vertex_curvatures computes curvature information at a single vertex v of a triangle mesh:
// -  It aggregates the second fundamental form coefficients from all faces incident to v.
// -  It converts those per-face coefficients into the vertex’s tangent frame.
// -  It averages them with area weights.
// -  From that, it later derives:
// -    Principal curvatures (min/max curvature values)
// -    Principal directions (3D directions of min/max curvature)
// -  It finally stores the principal direction vectors on v.
// This function is called from compute_vertex_curvatures for each vertex in the mesh.
//  
// Inputs and Outputs
// Inputs:
// -  half_edge_mesh<float>* mesh
// Full mesh structure. Needed to:
// -  Access incident faces of the vertex (via incident_faces + mesh->faces).
// -  Compute per-face weights via voronoi_based_weighting.
// -  vertex<float>& v
// The specific vertex for which we want curvature:
// -  normal is used to define the vertex tangent frame.
// -  incident_faces is the list of face indices touching this vertex.
// -  id is used by voronoi_based_weighting.
// Outputs:
// -  principal_dir_min.
// -  principal_dir_max – 3D direction of maximum principal curvature at the vertex.
static void calculate_vertex_curvatures(half_edge_mesh<float>* mesh, vertex<float>& v) {
    // -  sum_ee, sum_ff, sum_gg -> running area - weighted sums of second fundamental form components in the vertex frame :
    // -  ee -> curvature in local u direction
    // -  ff -> mixed term
    // -  gg -> curvature in local v direction
    // -  sum_of_V_areas -> sum of per - face vertex areas(weights) for normalization.
    float sum_ee = 0;
    float sum_ff = 0;
    float sum_gg = 0;
    float sum_of_V_areas = 0;

    // -  Given the vertex normal normal, create_uv_reference_plane constructs an orthonormal tangent basis:
    // -  u_vertex, v_vertex are tangent to the surface at the vertex and orthogonal to normal.
    // -  This (u_vertex, v_vertex) is the vertex tangent frame in which we ultimately want curvatures and directions.
    fvec3 u_vertex, v_vertex;
    create_uv_reference_plane(v.normal, u_vertex, v_vertex);

    // Iterate over all faces incident to the vertex
    for (auto& i : v.incident_faces) {
        // Get the face and its normal
        face<float>* n_face = mesh->faces[i];
        fvec3& face_normal = n_face->normal;

        // Create local tangent basis for the face
        // (u_face, v_face) is a 3D tangent basis aligned with the face
        fvec3 u_face, v_face;
        create_uv_reference_plane(face_normal, u_face, v_face);
        // -  ref_map encodes a 2x3 matrix mapping 3D vectors into the face's 2D (u_f, v_f) coordinate system:
        // -  Row 0: components of u_face
        // -  Row 1: components of v_face
        // So ref_map * [x y z]^T gives 2D coordinates of a 3D vector in the face’s tangent plane.
        float ref_map[] = { u_face.x(), u_face.y(), u_face.z(), v_face.x(), v_face.y(), v_face.z() };

        // Compute per face second fundamental form (face frame)
        // -  This function:
        // -    Projects the triangle into the face’s 2D coordinates.
        // -    Looks at normal differences between vertices on that face.
        // -    Uses a least squares fit to get [e, f, g] for that face, in the face frame (u_face, v_face).
        // After this call:
        // -  ee_ff_gg[0]: curvature along u_face.
        // -  ee_ff_gg[1]: mixed curvature term in face frame.
        // -  ee_ff_gg[2]: curvature along v_face.
        // Note: they are still in the face tangent frame, not yet in the vertex tangent frame.
        float ee_ff_gg[3];
        computeSecondFundamentalFormComponentsPerElem(mesh, *n_face, v, u_face, v_face, ee_ff_gg);
      
        // Transform face-frame [e, f, g] into vertex-frame [e_v, f_v, g_v]
        // -  convertSecondFundamentalFormComponentsToVertexFrame:
        // -    Expresses the vertex tangent basis u_vertex, v_vertex in the face’s coordinates (or a rotated version if face and vertex planes are not coplanar).
        // -    Applies the quadratic form basis change to re-express [e, f, g] in the vertex frame.
        // -  On output:
        // -    ee_ff_gg_v[0] = e_v (curvature along u_vertex)
        // -    ee_ff_gg_v[1] = f_v (mixed term in vertex frame)
        // -    ee_ff_gg_v[2] = g_v (curvature along v_vertex)
        // So these are per-face contributions already expressed in the vertex’s tangent frame.
        float ee_ff_gg_v[3];
        convertSecondFundamentalFormComponentsToVertexFrame(ref_map, u_vertex, v_vertex, u_face, v_face, v.normal, face_normal, ee_ff_gg, ee_ff_gg_v);

        // Area weighting and accumulation

        // -  v_area is an area-like weight for the pair (vertex v, face n_face):
        // -    Typically a Voronoi area or similar local area assigned to v from that face.
        // -    Encodes how much of the surface "near" v is contributed by this face.
        float v_area = voronoi_based_weighting<float>(*mesh, v.id, *n_face);

        // -  Each coefficient is summed with this area weight.
        // -  Conceptually: you’re computing an area-weighted average of the second fundamental form over the 1-ring of faces around v.        sum_ee += v_area * ee_ff_gg_v[0];
        sum_ff += v_area * ee_ff_gg_v[1];
        sum_gg += v_area * ee_ff_gg_v[2];

        sum_of_V_areas += v_area;
    }
    // Prevent division by zero in normalization
    if (!sum_of_V_areas)
        sum_of_V_areas = 1;
   
    // Normalize to get averaged second fundamental form at vertex
    // If a sum is exactly zero, it’s replaced by a small TOLLERANCE_SQ<float> 
    // to avoid degenerate zero matrices in later computations.
    float normalized_ee = (!sum_ee) ? TOLLERANCE_SQ<float> : (sum_ee / sum_of_V_areas);
    float normalized_ff = (!sum_ff) ? TOLLERANCE_SQ<float> : (sum_ff / sum_of_V_areas);
    float normalized_gg = (!sum_gg) ? TOLLERANCE_SQ<float> : (sum_gg / sum_of_V_areas);

    // Transform second fundamental form components into curvature data
    // W represents the entries of the 2x2 Weingarten (shape) matrix in the vertex tangent frame:
    // W = [ e  f ]
    //     [ f  g ]
    // 
    // where:
    //   e = W[0]
    //   f = W[1]
    //   g = W[2]
    float W[] = { normalized_ee,normalized_ff,normalized_gg };

    // Derive principal curvature values and curvature signs
    VertexCurvatureData<float>& curvature_data = v.curvature_data;

    // -  transformComponentsToCurvatureData:
    // -    Computes Gaussian curvature K = e*g - f².
    // -    Computes Mean curvature H = 0.5*(e + g).
    // -    Solves the eigenvalue problem of W to get the 2 principal curvatures (principal_curvs[0], principal_curvs[1]).
    // -    Determines sign information and absolute magnitudes.
    // At this point, we have the eigenvalues (principal curvature magnitudes) at the vertex.
    transformComponentsToCurvatureData(W, curvature_data);

    //Get position of Node
    //fvec3 Node_3d = v.coords;

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
void compute_vertex_curvatures(half_edge_mesh<float>* mesh) {
    // Precompute face properties and vertex normals
    mesh->compute_face_properties();
    mesh->compute_vertex_normals();

    // Compute curvature for each vertex in the mesh
    for (auto& v : mesh->vertices) {
        calculate_vertex_curvatures(mesh, *(v.second));
    }
    mesh->curvatures_computed() = true;
}