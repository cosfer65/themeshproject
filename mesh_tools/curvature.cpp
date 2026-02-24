#include "base_definitions.h"

#include "vector.h"

#include "mesh.h"
#include "mesh_curvature.h"

using namespace base_math;

static void worldCoordsToFaceCoords(const dvec3& u_f, const dvec3& v_f, const dvec3* vertices, const dvec3& centroid, dvec2* results)
{

    dvec3 proj_n1 = vertices[0] - centroid;
    dvec3 proj_n2 = vertices[1] - centroid;
    dvec3 proj_n3 = vertices[2] - centroid;

    double uvf[] = { u_f.x(), u_f.y(), u_f.z(), v_f.x(), v_f.y(), v_f.z() };

    results[0][0] = uvf[0] * proj_n1.x() + uvf[1] * proj_n1.y() + uvf[2] * proj_n1.z();
    results[0][1] = uvf[3] * proj_n1.x() + uvf[4] * proj_n1.y() + uvf[5] * proj_n1.z();

    results[1][0] = uvf[0] * proj_n2.x() + uvf[1] * proj_n2.y() + uvf[2] * proj_n2.z();
    results[1][1] = uvf[3] * proj_n2.x() + uvf[4] * proj_n2.y() + uvf[5] * proj_n2.z();

    results[2][0] = uvf[0] * proj_n3.x() + uvf[1] * proj_n3.y() + uvf[2] * proj_n3.z();
    results[2][1] = uvf[3] * proj_n3.x() + uvf[4] * proj_n3.y() + uvf[5] * proj_n3.z();
}

static void convert3dTo2dCoords(const dmat2x3& uv_map, const dvec2(&nodes_local)[3], const dvec3(&normals)[3], double* result)
{

    dvec3 diff_n2n1(normals[2] - normals[1]);
    dvec3 diff_n0n2(normals[0] - normals[2]);
    dvec3 diff_n1n0(normals[1] - normals[0]);

    dvec3 diff_n2_n1_3d(diff_n2n1.x(), diff_n2n1.y(), diff_n2n1.z());
    dvec3 diff_n0_n2_3d(diff_n0n2.x(), diff_n0n2.y(), diff_n0n2.z());
    dvec3 diff_n1_n0_3d(diff_n1n0.x(), diff_n1n0.y(), diff_n1n0.z());

    dvec2 diff_n2_n1_2d(0);
    dvec2 diff_n0_n2_2d(0);
    dvec2 diff_n1_n0_2d(0);

    diff_n2_n1_2d[0] = uv_map[0] * diff_n2_n1_3d.x() + uv_map[1] * diff_n2_n1_3d.y() + uv_map[2] * diff_n2_n1_3d.z();
    diff_n2_n1_2d[1] = uv_map[3] * diff_n2_n1_3d.x() + uv_map[4] * diff_n2_n1_3d.y() + uv_map[5] * diff_n2_n1_3d.z();

    diff_n0_n2_2d[0] = uv_map[0] * diff_n0_n2_3d.x() + uv_map[1] * diff_n0_n2_3d.y() + uv_map[2] * diff_n0_n2_3d.z();
    diff_n0_n2_2d[1] = uv_map[3] * diff_n0_n2_3d.x() + uv_map[4] * diff_n0_n2_3d.y() + uv_map[5] * diff_n0_n2_3d.z();

    diff_n1_n0_2d[0] = uv_map[0] * diff_n1_n0_3d.x() + uv_map[1] * diff_n1_n0_3d.y() + uv_map[2] * diff_n1_n0_3d.z();
    diff_n1_n0_2d[1] = uv_map[3] * diff_n1_n0_3d.x() + uv_map[4] * diff_n1_n0_3d.y() + uv_map[5] * diff_n1_n0_3d.z();

    dvec2 diff_nodes_n2n1 = nodes_local[2] - nodes_local[1];
    dvec2 diff_nodes_n0n2 = nodes_local[0] - nodes_local[2];
    dvec2 diff_nodes_n1n0 = nodes_local[1] - nodes_local[0];

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

static void prepareSecondFundamentalCoefficients(mesh<double>& fmesh, const meshFace<double>& cur_face, const dvec3& u_f, const dvec3& v_f, double* results)
{
    dvec3 elemCentroid = cur_face.center;
    dvec3 world_coords[3] = { cur_face.vertices[0]->position,
                                cur_face.vertices[1]->position, 
                                cur_face.vertices[2]->position };
    dvec2 local_coords[3];
    worldCoordsToFaceCoords(u_f, v_f, world_coords, elemCentroid, local_coords);
    const dvec3& node_normal_1 = cur_face.vertices[0]->normal;
    const dvec3& node_normal_2 = cur_face.vertices[1]->normal;
    const dvec3& node_normal_3 = cur_face.vertices[2]->normal;
    dmat2x3 uv_map({ u_f.x(), u_f.y(), u_f.z(), v_f.x(), v_f.y(), v_f.z() });
    convert3dTo2dCoords(uv_map, local_coords, { node_normal_1, node_normal_2, node_normal_3 }, results);
}

static void project_uv_to_face_2d(const basematrix<double, 2, 3>& map_mat, const dvec3& u_p, const dvec3& v_p,
    basematrix<double, 2, 1>& res1, basematrix<double, 2, 1>& res2)
{

    basematrix<double, 2, 1> res_mat1(map_mat * u_p);
    basematrix<double, 2, 1> res_mat2(map_mat * v_p);
    for (int r = 0; r < 2; ++r) {
        res1[r] = res_mat1[r];
        res2[r] = res_mat2[r];
    }
}

static void rotateElemPlane(dvec3& u_f, dvec3& v_f, dvec3& vertex_normal, dvec3& face_normal, basematrix<double, 2, 3>& result)
{

    dvec3 rotation_axis = vertex_normal.cross(face_normal).normalize();
    dvec3 rotated_u_f;
    dvec3 rotated_v_f;

    double norm_1 = vertex_normal.length();
    double norm_2 = face_normal.length();
    double theta = -acos((vertex_normal.dot(face_normal)) / (norm_1 * norm_2));
    double cth = cos(theta);
    double sth = sin(theta);

    double r11 = cth + rotation_axis.x() * rotation_axis.x() * (1 - cth);
    double r12 = rotation_axis.x() * rotation_axis.y() * (1 - cth) - rotation_axis.z() * sth;
    double r13 = rotation_axis.x() * rotation_axis.z() * (1 - cth) + rotation_axis.y() * sth;
    double r21 = rotation_axis.x() * rotation_axis.y() * (1 - cth) + rotation_axis.z() * sth;
    double r22 = cth + rotation_axis.y() * rotation_axis.y() * (1 - cth);
    double r23 = rotation_axis.z() * rotation_axis.y() * (1 - cth) - rotation_axis.x() * sth;
    double r31 = rotation_axis.z() * rotation_axis.x() * (1 - cth) - rotation_axis.y() * sth;
    double r32 = rotation_axis.z() * rotation_axis.y() * (1 - cth) + rotation_axis.x() * sth;
    double r33 = cth + rotation_axis.z() * rotation_axis.z() * (1 - cth);

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

static void computeSecondFundamentalFromComponents(mesh<double>* fmesh,
    meshFace<double>& n_face,            

    meshVertex<double>& currNode,        

    dvec3& u_f, dvec3& v_f,          

    double* e_f_g)                

{

    double fdifferences[12];

    prepareSecondFundamentalCoefficients(*fmesh, n_face, u_f, v_f, fdifferences);

    leastSquares(fdifferences, e_f_g);
}

static void convertSecondFundamentalFormComponentsToVertexFrame(
    double* iso_map,          

    dvec3& u_p, dvec3& v_p,  

    dvec3& u_f, dvec3& v_f,  

    dvec3& node_normal,      

    dvec3& elem_normal,      

    double* e_f_g,         

    double* e_f_g_V        

)
{

    basematrix<double, 2, 1> u_p_in_face_coords(0);
    basematrix<double, 2, 1> v_p_in_face_coords(0);

    if (isCoplanar<double>(u_f, v_f, u_p, v_p)) {

        project_uv_to_face_2d(iso_map, u_p, v_p, u_p_in_face_coords, v_p_in_face_coords);
    }
    else {

        basematrix<double, 2, 3> iso_map_new(0);

        rotateElemPlane(u_f, v_f, node_normal, elem_normal, iso_map_new);

        project_uv_to_face_2d(iso_map_new, u_p, v_p, u_p_in_face_coords, v_p_in_face_coords);
    }

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

    double ee = e_f_g[0];
    double ff = e_f_g[1];
    double gg = e_f_g[2];

    e_f_g_V[0] = u_p_f.x() * u_p_f.x() * ee + 2 * u_p_f.x() * u_p_f.y() * ff + u_p_f.y() * u_p_f.y() * gg;
    e_f_g_V[1] = v_p_f.x() * u_p_f.x() * ee
        + (u_p_f.y() * v_p_f.x() + u_p_f.x() * v_p_f.y()) * ff
        + (u_p_f.y() * v_p_f.y()) * gg;
    e_f_g_V[2] = v_p_f.x() * v_p_f.x() * ee + 2 * v_p_f.x() * v_p_f.y() * ff + v_p_f.y() * v_p_f.y() * gg;
}

void transformComponentsToCurvatureData(const double* efg, curvatureData<double>& curvature_data)
{

    double Gauss_curvature = efg[0] * efg[2] - efg[1] * efg[1];
    curvature_data.gaussCurvature = Gauss_curvature;
    curvature_data.absGaussCurvature = double(fabs(Gauss_curvature));

    double Mean_curvature = double(0.5 * (efg[0] + efg[2]));
    curvature_data.meanCurvature = Mean_curvature;

    std::tuple<double, double> res;
    solve_quadratic<double>(1, -2 * Mean_curvature, Gauss_curvature, res);
    curvature_data.k_min = std::get<0>(res);
    curvature_data.k_max = std::get<1>(res);
    if (curvature_data.k_min > curvature_data.k_max) {
        std::swap(curvature_data.k_min, curvature_data.k_max);
    }

    double mean = double(0.5 * (curvature_data.k_min + curvature_data.k_max));
    curvature_data.signGauss = 0;
    if (mean > TOLLERANCE<double>)
        curvature_data.signMean = 1;
    else if (mean < -TOLLERANCE<double>)
        curvature_data.signMean = -1;
    if (curvature_data.k_min * curvature_data.k_max > TOLLERANCE<double>)
        curvature_data.signGauss = 1;
    else if (curvature_data.k_min * curvature_data.k_max < -TOLLERANCE<double>)
        curvature_data.signGauss = -1;

    if (fabs(curvature_data.k_min) < fabs(curvature_data.k_max)) {
        curvature_data.absKmin = double(fabs(curvature_data.k_min));
        curvature_data.absKmax = double(fabs(curvature_data.k_max));
    }
    else {
        curvature_data.absKmax = double(fabs(curvature_data.k_min));
        curvature_data.absKmin = double(fabs(curvature_data.k_max));
    }
}

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

void calculatePrincipalDirections(double* W, double k_min, double k_max, dvec3& u_p, dvec3& v_p, dvec3& dmin, dvec3& dmax)
{

    double denom_1 = W[0] - k_min;
    double denom_2 = W[0] - k_max;

    if (fabs(denom_1) < TOLLERANCE<double>) denom_1 = TOLLERANCE<double>;
    if (fabs(denom_2) < TOLLERANCE<double>) denom_2 = TOLLERANCE<double>;

    double eigval_1 = -W[1] / denom_1;
    double eigval_2 = -W[1] / denom_2;

    if (fabs(k_min) <= fabs(k_max))
    {
        double princ_dir_min[] = { eigval_1, 1 };
        double princ_dir_max[] = { eigval_2, 1 };
        principalDirections2Dto3D(princ_dir_min, princ_dir_max, u_p, v_p, dmin, dmax);
    }
    else
    {
        double princ_dir_min[] = { eigval_2, 1 };
        double princ_dir_max[] = { eigval_1, 1 };
        principalDirections2Dto3D(princ_dir_min, princ_dir_max, u_p, v_p, dmin, dmax);
    }
}

static void calculate_vertex_curvatures(mesh<double>* fmesh, meshVertex<double>& v) {

    double sum_ee = 0;
    double sum_ff = 0;
    double sum_gg = 0;
    double sum_of_V_areas = 0;

    dvec3 u_vertex, v_vertex;
    create_uv_reference_plane(v.normal, u_vertex, v_vertex);

    for (auto& i : v.adjacentFaces) {

        meshFace<double>* n_face = i;// fmesh->faces[i->id];
        if (!n_face)
            continue;
        dvec3& face_normal = n_face->normal;

        dvec3 u_face, v_face;
        create_uv_reference_plane(face_normal, u_face, v_face);

        double ref_map[] = { u_face.x(), u_face.y(), u_face.z(), v_face.x(), v_face.y(), v_face.z() };

        double e_f_g[3];
        computeSecondFundamentalFromComponents(fmesh, *n_face, v, u_face, v_face, e_f_g);

        double e_f_g_v[3];
        convertSecondFundamentalFormComponentsToVertexFrame(ref_map, u_vertex, v_vertex, u_face, v_face, v.normal, face_normal, e_f_g, e_f_g_v);

        double v_area = voronoi_based_weighting<double>(*fmesh, v.id, *n_face);

        sum_ff += v_area * e_f_g_v[1];
        sum_gg += v_area * e_f_g_v[2];

        sum_of_V_areas += v_area;
    }

    if (!sum_of_V_areas)
        sum_of_V_areas = 1;

    double normalized_ee = (!sum_ee) ? TOLLERANCE_SQ<double> : (sum_ee / sum_of_V_areas);
    double normalized_ff = (!sum_ff) ? TOLLERANCE_SQ<double> : (sum_ff / sum_of_V_areas);
    double normalized_gg = (!sum_gg) ? TOLLERANCE_SQ<double> : (sum_gg / sum_of_V_areas);

    double W[] = { normalized_ee,normalized_ff,normalized_gg };

    curvatureData<double>& curvature_data = v.curvature_info;

    transformComponentsToCurvatureData(W, v.curvature_info);

    calculatePrincipalDirections(W, v.curvature_info.k_min, v.curvature_info.k_max, u_vertex, v_vertex,
        v.curvature_info.k_min_dir, v.curvature_info.k_max_dir);

    v.curvature_info.k_min_dir.normalize();
    v.curvature_info.k_max_dir.normalize();
}

void compute_vertex_curvatures(mesh<double>* fmesh) {

    fmesh->computeFaceProperties();
    fmesh->computeVertexNormals();

    fmesh->curvatures_computed() = false;

    for (auto& v : fmesh->vertices) {
        calculate_vertex_curvatures(fmesh, *(v.second));
    }
    fmesh->curvatures_computed() = true;
}

