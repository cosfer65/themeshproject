#include <vector>
#include <cmath>
#include "mesh.h"
#include "vector.h"

using namespace base_math;

static double computeCotangentWeight(meshVertex<double>* v1, meshVertex<double>* v2) {
    // Find the two triangles that share edge (i,j)
    double cot_sum = 0.0;

    for (meshFace<double>* face : v1->adjacentFaces) {
        // Check if this triangle contains both vertices i and j
        bool has_i = false, has_j = false;
        meshVertex<double>* other_vertex = nullptr;

        for (int k = 0; k < 3; ++k) {
            if (face->vertices[k] == v1) has_i = true;
            else if (face->vertices[k] == v2) has_j = true;
            else other_vertex = face->vertices[k];
        }

        if (has_i && has_j && other_vertex != nullptr) {
            dvec3 vi = v1->position;
            dvec3 vj = v2->position;
            dvec3 vk = other_vertex->position;

            dvec3 e1 = vk - vi;
            dvec3 e2 = vk - vj;
            double cos_angle = e1.normalize().dot(e2.normalize());
            double sin_angle = e1.cross(e2).length() / (e1.length() * e2.length());
            cot_sum += cos_angle / sin_angle;
        }
    }

    return cot_sum / 2.0;
}

static void computeVertexCurvature(meshVertex<double>* v) {
    // fetch angle defect and mixed area
    double area_mixed = v->voronoiArea;
    double angle_defect = 2.0 * PI<double> -v->angle_sum;
    if (v->is_boundary) {
        // For boundary vertices, the angle defect is different
        angle_defect = PI<double> -v->angle_sum;
    }

    v->curvature_info.gaussCurvature = angle_defect / area_mixed;
    v->curvature_info.absGaussCurvature = fabs(v->curvature_info.gaussCurvature);
    if (v->curvature_info.absGaussCurvature < 0.02) {
        v->curvature_info.signGauss = 0;
    }
    else
        v->curvature_info.signGauss = (v->curvature_info.gaussCurvature > 0) ? 1 : ((v->curvature_info.gaussCurvature < 0) ? -1 : 0);

    // T absMeanCurvature;
    // T absGaussCurvature;
    // int signGauss;
    // int signMean;
    // Compute mean curvature using discrete Laplace-Beltrami operator
    dvec3 laplace(0, 0, 0);

    for (meshVertex<double>* neighbor : v->neighbourVertices) {
        double cot_weight = computeCotangentWeight(v, neighbor);
        laplace += cot_weight * (neighbor->position - v->position);
    }

    // Mean curvature = 0.5 * ||Laplace-Beltrami||
    v->curvature_info.meanCurvature = 0.5 * laplace.dot(v->normal) / area_mixed;
    v->curvature_info.meanCurvatureDir = v->normal * v->curvature_info.meanCurvature * 2;
}

/////////////////////////////////////////////////////////////
// symmetric eigen decomposition
// S = [a b; b c], compute eigenvalues/eigenvectors
static void eigenSymmetric2x2(double a, double b, double c,
    double& l1, double& l2,
    double& v1x, double& v1y,
    double& v2x, double& v2y)
{
    // trace and determinant
    double tr = a + c;
    double det = a * c - b * b;
    double disc = tr * tr - 4.0 * det;
    if (disc < 0.0) disc = 0.0;
    double s = std::sqrt(disc);

    l1 = 0.5 * (tr + s);
    l2 = 0.5 * (tr - s);

    // eigenvector for l1
    if (std::fabs(b) > 1e-12) {
        v1x = l1 - c;
        v1y = b;
    }
    else {
        // matrix is diagonal or nearly
        v1x = 1.0;
        v1y = 0.0;
    }
    double n1 = std::sqrt(v1x * v1x + v1y * v1y);
    if (n1 > 0.0) { v1x /= n1; v1y /= n1; }

    // eigenvector for l2: orthogonal to v1
    v2x = -v1y;
    v2y = v1x;
}


static void computePrincipalDirections(meshVertex<double>* v)    // input vertex for which to compute principal directions
{
    double& k1 = v->curvature_info.k_min;
    double& k2 = v->curvature_info.k_max;
    basevector<double, 3>& dir1 = v->curvature_info.k_min_dir;
    basevector<double, 3>& dir2 = v->curvature_info.k_max_dir;

    // Build tangent frame (t1, t2, n)
    basevector<double, 3>& pi = v->position;
    basevector<double, 3>& ni = v->normal;

    basevector<double, 3> ref = (std::fabs(ni.x()) > std::fabs(ni.z()))
        ? basevector<double, 3>(-ni.y(), ni.x(), 0.0)
        : basevector<double, 3>(0.0, -ni.z(), ni.y());

    basevector<double, 3> t1(ref); t1.normalize();
    basevector<double, 3> t2(cross(ni, t1)); t2.normalize();

    // Accumulate least-squares system A^T A x = A^T b
    // Unknowns: x = [s11, s12, s22]^T
    double ATA[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
    double ATb[3] = { 0,0,0 };

    std::vector<meshVertex<double>*>& nbrs = v->neighbourVertices; //

    if (nbrs.size() < 3) {
        // Not enough neighbors to fit a tensor
        k1 = k2 = 0.0;
        dir1 = t1;
        dir2 = t2;
        return;
    }

    for (meshVertex<double>* vert : nbrs) {
        const basevector<double, 3>& pj = vert->position;
        const basevector<double, 3>& nj = vert->normal;

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
    } // (for (meshVertex<double>* vert : nbrs))

    // Symmetrize ATA
    ATA[1][0] = ATA[0][1];
    ATA[2][0] = ATA[0][2];
    ATA[2][1] = ATA[1][2];

    // Solve 3x3 system ATA * x = ATb (simple Gaussian elimination)
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
        double s11 = M[0][3];
        double s12 = M[1][3];
        double s22 = M[2][3];

        // Eigen decomposition of S = [s11 s12; s12 s22]
        double l1, l2, v1x, v1y, v2x, v2y;
        eigenSymmetric2x2(s11, s12, s22, l1, l2, v1x, v1y, v2x, v2y);

        // Map eigenvectors back to 3D
        basevector<double, 3> d1 = v1x * t1 + v1y * t2;
        basevector<double, 3> d2 = v2x * t1 + v2y * t2;

        // Store (you may want to reorder so |k1| >= |k2|)
        k1 = l1;
        k2 = l2;
        dir1 = d1.normalize();
        dir2 = d2.normalize();
    }

}

bool computeMeshCurvatures(mesh<double>& fmesh) {
    for (auto& v : fmesh.vertices) {
        computeVertexCurvature(v.second);
        computePrincipalDirections(v.second);
    }
    return true;
}

bool computeMeshCurvatures(mesh<double>& fmesh);
void compute_vertex_curvatures(mesh<double>* fmesh) {
    fmesh->computeFaceProperties();
    fmesh->computeVertexNormals();

    fmesh->curvatures_computed() = false;

    // for (auto& v : fmesh->vertices) {
    //     calculate_vertex_curvatures(fmesh, *(v.second));
    // }
    computeMeshCurvatures(*fmesh);
    fmesh->curvatures_computed() = true;
}