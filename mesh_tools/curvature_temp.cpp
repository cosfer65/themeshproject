#include "base_definitions.h"

#include "vector.h"

#include "mesh.h"
#include "mesh_curvature.h"

using namespace base_math;

static void WorldCoordsToFaceCoords(const fvec3& u_f, const fvec3& v_f, const fvec3* vertices, const fvec3& centroid, fvec2* results)
{

    fvec3 proj_n1 = vertices[0] - centroid;
    fvec3 proj_n2 = vertices[1] - centroid;
    fvec3 proj_n3 = vertices[2] - centroid;

    float uvf[] = { u_f.x(), u_f.y(), u_f.z(), v_f.x(), v_f.y(), v_f.z() };

    results[0][0] = uvf[0] * proj_n1.x() + uvf[1] * proj_n1.y() + uvf[2] * proj_n1.z();
    results[0][1] = uvf[3] * proj_n1.x() + uvf[4] * proj_n1.y() + uvf[5] * proj_n1.z();

    results[1][0] = uvf[0] * proj_n2.x() + uvf[1] * proj_n2.y() + uvf[2] * proj_n2.z();
    results[1][1] = uvf[3] * proj_n2.x() + uvf[4] * proj_n2.y() + uvf[5] * proj_n2.z();

    results[2][0] = uvf[0] * proj_n3.x() + uvf[1] * proj_n3.y() + uvf[2] * proj_n3.z();
    results[2][1] = uvf[3] * proj_n3.x() + uvf[4] * proj_n3.y() + uvf[5] * proj_n3.z();
}

static void convert3dTo2dCoords(const fmat2x3& uv_map, const fvec2(&nodes_local)[3], const fvec3(&normals)[3], float* result)
{

    fvec3 diff_n2n1(normals[2] - normals[1]);
    fvec3 diff_n0n2(normals[0] - normals[2]);
    fvec3 diff_n1n0(normals[1] - normals[0]);

    fvec3 diff_n2_n1_3d(diff_n2n1.x(), diff_n2n1.y(), diff_n2n1.z());
    fvec3 diff_n0_n2_3d(diff_n0n2.x(), diff_n0n2.y(), diff_n0n2.z());
    fvec3 diff_n1_n0_3d(diff_n1n0.x(), diff_n1n0.y(), diff_n1n0.z());

    fvec2 diff_n2_n1_2d(0);
    fvec2 diff_n0_n2_2d(0);
    fvec2 diff_n1_n0_2d(0);

    diff_n2_n1_2d[0] = uv_map[0] * diff_n2_n1_3d.x() + uv_map[1] * diff_n2_n1_3d.y() + uv_map[2] * diff_n2_n1_3d.z();
    diff_n2_n1_2d[1] = uv_map[3] * diff_n2_n1_3d.x() + uv_map[4] * diff_n2_n1_3d.y() + uv_map[5] * diff_n2_n1_3d.z();

    diff_n0_n2_2d[0] = uv_map[0] * diff_n0_n2_3d.x() + uv_map[1] * diff_n0_n2_3d.y() + uv_map[2] * diff_n0_n2_3d.z();
    diff_n0_n2_2d[1] = uv_map[3] * diff_n0_n2_3d.x() + uv_map[4] * diff_n0_n2_3d.y() + uv_map[5] * diff_n0_n2_3d.z();

    diff_n1_n0_2d[0] = uv_map[0] * diff_n1_n0_3d.x() + uv_map[1] * diff_n1_n0_3d.y() + uv_map[2] * diff_n1_n0_3d.z();
    diff_n1_n0_2d[1] = uv_map[3] * diff_n1_n0_3d.x() + uv_map[4] * diff_n1_n0_3d.y() + uv_map[5] * diff_n1_n0_3d.z();

    fvec2 diff_nodes_n2n1 = nodes_local[2] - nodes_local[1];
    fvec2 diff_nodes_n0n2 = nodes_local[0] - nodes_local[2];
    fvec2 diff_nodes_n1n0 = nodes_local[1] - nodes_local[0];

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
    fvec3 world_coords[3] = { mesh.vertices[cur_face.v1]->coords, 
	                          mesh.vertices[cur_face.v2]->coords,
							  mesh.vertices[cur_face.v3]->coords };
    fvec2 local_coords[3];
    WorldCoordsToFaceCoords(u_f, v_f, world_coords, elemCentroid, local_coords);
    const fvec3& node_normal_1 = mesh.vertices[cur_face.v1]->normal;
    const fvec3& node_normal_2 = mesh.vertices[cur_face.v2]->normal;
    const fvec3& node_normal_3 = mesh.vertices[cur_face.v3]->normal;
    fmat2x3 uv_map({ u_f.x(), u_f.y(), u_f.z(), v_f.x(), v_f.y(), v_f.z() });
    convert3dTo2dCoords(uv_map, local_coords, { node_normal_1, node_normal_2, node_normal_3 }, results);
}

static void leastSquares(const float* second_fundamental, float* ee_ff_gg)
{

    basematrix<float, 6, 1> N({ second_fundamental[0],second_fundamental[1],second_fundamental[2],
                    second_fundamental[3],second_fundamental[4],second_fundamental[5] });

    basematrix<float, 6, 3> E({ second_fundamental[6], second_fundamental[7],0,0,second_fundamental[6],second_fundamental[7] ,
                   second_fundamental[8], second_fundamental[9],0,0,second_fundamental[8],second_fundamental[9] ,
                  second_fundamental[10], second_fundamental[11],0,0,second_fundamental[10],second_fundamental[11] });

    basematrix<float, 3, 6> Trans_E = E.transpose();

    basematrix<float, 3, 3> tr_x_e = Trans_E * E;           

    basematrix<float, 3, 1> rhs = Trans_E * N;              

    basematrix<float, 3, 1> sol;                            

    if (!cholesky_solve_3x3<float>(tr_x_e, rhs, sol)) {

        basematrix<float, 3, 3> reg = tr_x_e;
        for (int i = 0; i < 3; ++i) reg[i * 3 + i] += TOLLERANCE<float>;
        cholesky_solve_3x3<float>(reg, rhs, sol);
    }

    ee_ff_gg[0] = sol[0];
    ee_ff_gg[1] = sol[1];
    ee_ff_gg[2] = sol[2];
}

static void project_uv_to_face_2d(const basematrix<float, 2, 3>& map_mat, const fvec3& u_p, const fvec3& v_p,
    basematrix<float, 2, 1>& res1, basematrix<float, 2, 1>& res2)
{

    basematrix<float, 2, 1> res_mat1(map_mat * u_p);
    basematrix<float, 2, 1> res_mat2(map_mat * v_p);
    for (int r = 0; r < 2; ++r) {
        res1[r] = res_mat1[r];
        res2[r] = res_mat2[r];
    }
}

static void rotateElemPlane(fvec3& u_f, fvec3& v_f, fvec3& vertex_normal, fvec3& face_normal, basematrix<float, 2, 3>& result)
{

    fvec3 rotation_axis = (vertex_normal * face_normal).normalize();
    fvec3 rotated_u_f;
    fvec3 rotated_v_f;

    float norm_1 = vertex_normal.length();
    float norm_2 = face_normal.length();
    float theta = -acos((vertex_normal.dot(face_normal)) / (norm_1 * norm_2));
    float cth = cos(theta);
    float sth = sin(theta);

    float r11 = cth + rotation_axis.x() * rotation_axis.x() * (1 - cth);
    float r12 = rotation_axis.x() * rotation_axis.y() * (1 - cth) - rotation_axis.z() * sth;
    float r13 = rotation_axis.x() * rotation_axis.z() * (1 - cth) + rotation_axis.y() * sth;
    float r21 = rotation_axis.x() * rotation_axis.y() * (1 - cth) + rotation_axis.z() * sth;
    float r22 = cth + rotation_axis.y() * rotation_axis.y() * (1 - cth);
    float r23 = rotation_axis.z() * rotation_axis.y() * (1 - cth) - rotation_axis.x() * sth;
    float r31 = rotation_axis.z() * rotation_axis.x() * (1 - cth) - rotation_axis.y() * sth;
    float r32 = rotation_axis.z() * rotation_axis.y() * (1 - cth) + rotation_axis.x() * sth;
    float r33 = cth + rotation_axis.z() * rotation_axis.z() * (1 - cth);

    rotated_u_f[0] = u_f.x() * r11 + u_f.y() * r12 + u_f.z() * r13;
    rotated_u_f[1] = u_f.x() * r21 + u_f.y() * r22 + u_f.z() * r23;
    rotated_u_f[2] = u_f.x() * r31 + u_f.y() * r32 + u_f.z() * r33;
    rotated_v_f[0] = v_f.x() * r11 + v_f.y() * r12 + v_f.z() * r13;
    rotated_v_f[1] = v_f.x() * r21 + v_f.y() * r22 + v_f.z() * r23;
    rotated_v_f[2] = v_f.x() * r31 + v_f.y() * r32 + v_f.z() * r33;
    rotated_u_f.normalize();
    rotated_v_f.normalize();

    result[0] = rotated_u_f.x();
    result[1] = rotated_u_f.y();
    result[2] = rotated_u_f.z();
    result[3] = rotated_v_f.x();
    result[4] = rotated_v_f.y();
    result[5] = rotated_v_f.z();
}

static void computeSecondFundamentalFormComponentsPerElem(half_edge_mesh<float>* mesh,
    face<float>& n_face,            

    vertex<float>& currNode,        

    fvec3& u_f, fvec3& v_f,         

    float* ee_ff_gg)                

{

    float form_II[12];

    prepareSeconFundamentalCoefficients(*mesh, n_face, u_f, v_f, form_II);

    leastSquares(form_II, ee_ff_gg);
}

static void convertSecondFundamentalFormComponentsToVertexFrame(
    float* iso_map,          

    fvec3& u_p, fvec3& v_p,  

    fvec3& u_f, fvec3& v_f,  

    fvec3& node_normal,      

    fvec3& elem_normal,      

    float* ee_ff_gg,         

    float* ee_ff_gg_V        

)
{

    basematrix<float, 2, 1> u_p_in_face_coords(0);
    basematrix<float, 2, 1> v_p_in_face_coords(0);

    if (isCoplanar<float>(u_f, v_f, u_p, v_p)) {

        project_uv_to_face_2d(iso_map, u_p, v_p, u_p_in_face_coords, v_p_in_face_coords);
    }
    else {

        basematrix<float, 2, 3> iso_map_new(0);

        rotateElemPlane(u_f, v_f, node_normal, elem_normal, iso_map_new);

        project_uv_to_face_2d(iso_map_new, u_p, v_p, u_p_in_face_coords, v_p_in_face_coords);
    }

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

    float ee = ee_ff_gg[0];
    float ff = ee_ff_gg[1];
    float gg = ee_ff_gg[2];

    ee_ff_gg_V[0] = u_p_f.x() * u_p_f.x() * ee + 2 * u_p_f.x() * u_p_f.y() * ff + u_p_f.y() * u_p_f.y() * gg;
    ee_ff_gg_V[1] = v_p_f.x() * u_p_f.x() * ee
        + (u_p_f.y() * v_p_f.x() + u_p_f.x() * v_p_f.y()) * ff
        + (u_p_f.y() * v_p_f.y()) * gg;
    ee_ff_gg_V[2] = v_p_f.x() * v_p_f.x() * ee + 2 * v_p_f.x() * v_p_f.y() * ff + v_p_f.y() * v_p_f.y() * gg;
}

void transformComponentsToCurvatureData(const float* accum_II, VertexCurvatureData<float>& curvature_data)

{

    float Gauss_curvature = accum_II[0] * accum_II[2] - accum_II[1] * accum_II[1];
    curvature_data.gaussCurvature = Gauss_curvature;

    float Mean_curvature = float(0.5 * (accum_II[0] + accum_II[2]));
    curvature_data.meanCurvature = Mean_curvature;

    std::tuple<float, float> res;
    solve_quadratic<float>(1, -2 * Mean_curvature, Gauss_curvature, res);
    curvature_data.principal_curvatures[0] = std::get<0>(res);
    curvature_data.principal_curvatures[1] = std::get<1>(res);
    if (curvature_data.principal_curvatures[0] > curvature_data.principal_curvatures[1]) {
        std::swap(curvature_data.principal_curvatures[0], curvature_data.principal_curvatures[1]);
    }

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

    if (fabs(curvature_data.principal_curvatures[0]) < fabs(curvature_data.principal_curvatures[1])) {
        curvature_data.absKmin = float(fabs(curvature_data.principal_curvatures[0]));
        curvature_data.absKmax = float(fabs(curvature_data.principal_curvatures[1]));
    }
    else {
        curvature_data.absKmax = float(fabs(curvature_data.principal_curvatures[0]));
        curvature_data.absKmin = float(fabs(curvature_data.principal_curvatures[1]));
    }
}

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

void calculatePrincipalDirections(float* W, float* princ, fvec3& u_p, fvec3& v_p, fvec3& dmin, fvec3& dmax)
{

    float denom_1 = W[0] - princ[0];
    float denom_2 = W[0] - princ[1];

    if (fabs(denom_1) < TOLLERANCE<float>) denom_1 = TOLLERANCE<float>;
    if (fabs(denom_2) < TOLLERANCE<float>) denom_2 = TOLLERANCE<float>;

    float eigvec_1 = -W[1] / denom_1;
    float eigvec_2 = -W[1] / denom_2;

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

static void calculate_vertex_curvatures(half_edge_mesh<float>* mesh, vertex<float>& v) {

    float sum_ee = 0;
    float sum_ff = 0;
    float sum_gg = 0;
    float sum_of_V_areas = 0;

    fvec3 u_vertex, v_vertex;
    create_uv_reference_plane(v.normal, u_vertex, v_vertex);

    for (auto& i : v.incident_faces) {

        face<float>* n_face = mesh->faces[i];
        fvec3& face_normal = n_face->normal;

        fvec3 u_face, v_face;
        create_uv_reference_plane(face_normal, u_face, v_face);

        float ref_map[] = { u_face.x(), u_face.y(), u_face.z(), v_face.x(), v_face.y(), v_face.z() };

        float ee_ff_gg[3];
        computeSecondFundamentalFormComponentsPerElem(mesh, *n_face, v, u_face, v_face, ee_ff_gg);

        float ee_ff_gg_v[3];
        convertSecondFundamentalFormComponentsToVertexFrame(ref_map, u_vertex, v_vertex, u_face, v_face, v.normal, face_normal, ee_ff_gg, ee_ff_gg_v);

        float v_area = voronoi_based_weighting<float>(*mesh, v.id, *n_face);

        sum_ff += v_area * ee_ff_gg_v[1];
        sum_gg += v_area * ee_ff_gg_v[2];

        sum_of_V_areas += v_area;
    }

    if (!sum_of_V_areas)
        sum_of_V_areas = 1;

    float normalized_ee = (!sum_ee) ? TOLLERANCE_SQ<float> : (sum_ee / sum_of_V_areas);
    float normalized_ff = (!sum_ff) ? TOLLERANCE_SQ<float> : (sum_ff / sum_of_V_areas);
    float normalized_gg = (!sum_gg) ? TOLLERANCE_SQ<float> : (sum_gg / sum_of_V_areas);

    float W[] = { normalized_ee,normalized_ff,normalized_gg };

    VertexCurvatureData<float>& curvature_data = v.curvature_data;

    transformComponentsToCurvatureData(W, curvature_data);

    calculatePrincipalDirections(W, curvature_data.principal_curvatures, u_vertex, v_vertex, 
        curvature_data.principal_directions[0], curvature_data.principal_directions[1]);
}

void compute_vertex_curvatures(half_edge_mesh<float>* mesh) {

    mesh->compute_face_properties();
    mesh->compute_vertex_normals();

    for (auto& v : mesh->vertices) {
        calculate_vertex_curvatures(mesh, *(v.second));
    }
}