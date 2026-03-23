#ifndef __nurbs_core_h__
#define __nurbs_core_h__

#include "vector.h"
#include "mesh.h"
#include "string_utils.h"

namespace base_math {

    template <typename T>
    struct NURBSSurface {
        std::vector<basevector<T, 4>> controlPoints;
        std::vector<T> knotsU, knotsV;
        int degreeU;
        int degreeV;

        int numCtrlPtsU, numCtrlPtsV;
        basematrix<T, 4, 4> world_matrix;
        basevector<T, 3> translation;
        basevector<T, 3> rotation; // Euler angles in radians
        basevector<T, 3> scale;
        const basevector<T, 4>& ctrl(int i, int j) const {
            return controlPoints[j * numCtrlPtsU + i];
        }
        T getWeight(int i, int j) const {
            return ctrl(i, j).w();
        }
        int numControlPointsU() const { return numCtrlPtsU; }
        int numControlPointsV() const { return numCtrlPtsV; }
    };

    template <typename T>
    std::vector<mesh<T>*> create_from_nurbs_file(const std::string& fnm, int resolutionU, int resolutionV) {
        std::vector<mesh<T>*> he_meshes;
        if (file_extension(fnm) != "nurbs")
            return he_meshes;

        std::vector<NURBSSurface<T>*> patches;
        // NURBSSurface<T> surf;
        load_nurbs<T>(fnm, patches);
        for (NURBSSurface<T>* surf : patches) {
            TessellatedMesh<T>* m = tessellateNURBSSurface<T>(*surf, resolutionU, resolutionV);// Convert NURBS surface to mesh
            mesh<T>* he_mesh = parse_mesh<T>(*m);
            // he_mesh->transform(surf->world_matrix); // Apply the original transformation of the NURBS surface to the mesh
            he_mesh->translate(surf->translation * 4); // Apply the original translation of the NURBS surface to the mesh);
            he_meshes.push_back(he_mesh);
            delete m; // Clean up the tessellated mesh
        }
        return he_meshes;
    }

    template <typename T>
    bool load_nurbs(const std::string& fnm, std::vector<NURBSSurface<T>*>& patches) {
        std::string line;
        std::ifstream mdl_file(fnm);
        NURBSSurface<T>* surf = nullptr;

        if (mdl_file.is_open()) {
            while (std::getline(mdl_file, line)) {
                if (line.empty())
                    continue;
                auto tokens = splitString(line);
                if (tokens.empty())
                    continue;
                if (tokens[0] == "END_OBJECT") {
                    surf = nullptr;
                    continue;
                }
                if (tokens[0] == "START_OBJECT") {
                    surf = new NURBSSurface<T>;
                    patches.push_back(surf);
                    continue;
                }
                if (tokens[0] == "matrix_world") {
                    // Parse world matrix
                    for (int i = 0; i < 4; ++i) {
                        for (int j = 0; j < 4; ++j) {
                            surf->world_matrix(i, j) = std::stod(tokens[1 + i * 4 + j]);
                        }
                    }
                    continue;
                }
                if (tokens[0] == "location") {
                    surf->translation.x() = std::stod(tokens[1]);
                    surf->translation.y() = std::stod(tokens[2]);
                    surf->translation.z() = std::stod(tokens[3]);
                    continue;
                }
                if (tokens[0] == "rotation") {
                    surf->rotation.x() = std::stod(tokens[1]);
                    surf->rotation.y() = std::stod(tokens[2]);
                    surf->rotation.z() = std::stod(tokens[3]);
                    continue;
                }
                if (tokens[0] == "scale") {
                    surf->scale.x() = std::stod(tokens[1]);
                    surf->scale.y() = std::stod(tokens[2]);
                    surf->scale.z() = std::stod(tokens[3]);
                    continue;
                }
                if (tokens[0] == "degree") {
                    surf->degreeU = std::stoul(tokens[1]);
                    surf->degreeV = std::stoul(tokens[2]);
                }
                else if (tokens[0] == "knots") {
                    // U/V index followed by knot values
                    char dir = tokens[1][0]; // 'U' or 'V'
                    std::vector<T> knots;
                    for (int i = 2; i < tokens.size(); ++i) {
                        knots.push_back(std::stod(tokens[i]));
                    }
                    if (dir == 'u')
                        surf->knotsU = std::move(knots);
                    else if (dir == 'v')
                        surf->knotsV = std::move(knots);
                }
                else if (tokens[0] == "num_ctrl_pts") {
                    // U V count
                    surf->numCtrlPtsU = std::stoul(tokens[1]);
                    surf->numCtrlPtsV = std::stoul(tokens[2]);
                }
                else if (tokens[0] == "point") {
                    // U V X Y Z W
                    int u = std::stoul(tokens[1]);
                    int v = std::stoul(tokens[2]);
                    T x = std::stod(tokens[3]);
                    T y = std::stod(tokens[4]);
                    T z = std::stod(tokens[5]);
                    T w = std::stod(tokens[6]);
                    surf->controlPoints.push_back(basevector<T, 4>(x, y, z, w));
                }
            }
            mdl_file.close();
        }
        else {
            return false;
        }

        return true;
    }

    // ------------------------------------------------------------
    // Cox–de Boor basis function
    // ------------------------------------------------------------
    template <typename T>
    T bsplineBasis(int i, int p, T u, const std::vector<T>& knots)
    {
        if (p == 0) {
            return (u >= knots[i] && u < knots[i + 1]) ? T(1) : T(0);
        }

        T left = T(0), right = T(0);

        T denom1 = knots[i + p] - knots[i];
        T denom2 = knots[i + p + 1] - knots[i + 1];
        if (denom1 > 1e-14)
            left = (u - knots[i]) / denom1 * bsplineBasis(i, p - 1, u, knots);

        if (denom2 > 1e-14)
            right = (knots[i + p + 1] - u) / denom2 * bsplineBasis(i + 1, p - 1, u, knots);

        return left + right;
    }

    // ------------------------------------------------------------
    // Map normalized [0,1] → true knot domain [k_p, k_{m-p-1}]
    // ------------------------------------------------------------
    template <typename T>
    T mapToKnotDomain(T t, const std::vector<T>& knots, int order)
    {
        int p = order - 1;
        int m = static_cast<int>(knots.size()) - 1;

        T umin = knots[p];
        T umax = knots[m - p - 1];

        if (t < T(0)) t = T(0);
        if (t > T(1)) t = T(1);

        return umin + t * (umax - umin);
    }

    // ------------------------------------------------------------
    // Clamp inside last knot span to avoid evaluating at u = umax
    // ------------------------------------------------------------
    template <typename T>
    T clampInsideDomain(T u, const std::vector<T>& knots, int order)
    {
        int p = order - 1;
        int m = static_cast<int>(knots.size()) - 1;

        T umin = knots[p];
        T umax = knots[m - p - 1];

        if (u >= umax) {
            T span = umax - umin;
            u = umax - T(1e-6) * span;
        }
        if (u < umin) {
            u = umin;
        }
        return u;
    }

    // ------------------------------------------------------------
    // Full NURBS surface evaluation
    // ------------------------------------------------------------
    template <typename T>
    basevector<T, 3> evaluateSurface(
        T u_norm, T v_norm,
        const std::vector<basevector<T, 4>>& ctrl,
        int nu, int nv,
        int orderU, int orderV,
        const std::vector<T>& knotU,
        const std::vector<T>& knotV)
    {
        // Map normalized parameters to true knot domain
        T u = mapToKnotDomain(u_norm, knotU, orderU);
        T v = mapToKnotDomain(v_norm, knotV, orderV);
        // Clamp inside domain to avoid u=umax or v=vmax
        u = clampInsideDomain(u, knotU, orderU);
        v = clampInsideDomain(v, knotV, orderV);

        int p = orderU - 1;
        int q = orderV - 1;

        T denom = T(0);
        T x = T(0), y = T(0), z = T(0);

        // Blender stores control points in row-major order: P[i * nv + j]
        for (int i = 0; i < nu; ++i) {
            T Nu = bsplineBasis(i, p, u, knotU);

            for (int j = 0; j < nv; ++j) {
                T Mv = bsplineBasis(j, q, v, knotV);

                const basevector<T, 4>& P = ctrl[i * nv + j];
                T w = P.w();

                T B = Nu * Mv * w;

                denom += B;
                x += B * P.x();
                y += B * P.y();
                z += B * P.z();
            }
        }

        return { x / denom, y / denom, z / denom };
    }

    template <typename T>
    TessellatedMesh<T>* tessellateNURBSSurface(const NURBSSurface<T>& surf, int samplesU, int samplesV)
    {
        TessellatedMesh<T>* mesh = new TessellatedMesh<T>;

        // ------------------------------------------------------------
        // 1. Generate UV grid and evaluate surface
        // ------------------------------------------------------------
        for (int i = 0; i <= samplesU; ++i) {
            T u = T(i) / T(samplesU);

            for (int j = 0; j <= samplesV; ++j) {
                T v = T(j) / T(samplesV);

                basevector<T, 3> p;

                if (i == 0)      p = evaluateSurface(0.0, v, surf.controlPoints,
                    surf.numCtrlPtsU, surf.numCtrlPtsV,
                    surf.degreeU, surf.degreeV,
                    surf.knotsU, surf.knotsV);
                else if (i == samplesU) p = evaluateSurface(1.0, v, surf.controlPoints,
                    surf.numCtrlPtsU, surf.numCtrlPtsV,
                    surf.degreeU, surf.degreeV,
                    surf.knotsU, surf.knotsV);
                else if (j == 0)      p = evaluateSurface(u, 0.0, surf.controlPoints,
                    surf.numCtrlPtsU, surf.numCtrlPtsV,
                    surf.degreeU, surf.degreeV,
                    surf.knotsU, surf.knotsV);
                else if (j == samplesV) p = evaluateSurface(u, 1.0, surf.controlPoints,
                    surf.numCtrlPtsU, surf.numCtrlPtsV,
                    surf.degreeU, surf.degreeV,
                    surf.knotsU, surf.knotsV);
                else
                    p = evaluateSurface(u, v, surf.controlPoints,
                        surf.numCtrlPtsU, surf.numCtrlPtsV,
                        surf.degreeU, surf.degreeV,
                        surf.knotsU, surf.knotsV);
                /*
                basevector<T, 3> p = evaluateSurface(
                    u, v,
                    surf.controlPoints,
                    surf.numCtrlPtsU, surf.numCtrlPtsV,
                    surf.degreeU, surf.degreeV,
                    surf.knotsU, surf.knotsV
                );
                */

                mesh->vertices.push_back(p);
                mesh->normals.push_back(basevector<T, 3>(0, 1, 0)); // Placeholder normal
            }
        }

        // ------------------------------------------------------------
        // 2. Build triangle connectivity
        // ------------------------------------------------------------
        auto idx = [&](int i, int j) {
            return i * (samplesV + 1) + j;
            };

        for (int i = 0; i < samplesU; ++i) {
            for (int j = 0; j < samplesV; ++j) {

                int v00 = idx(i, j);
                int v10 = idx(i + 1, j);
                int v01 = idx(i, j + 1);
                int v11 = idx(i + 1, j + 1);

                // Triangle 1
                mesh->indices.push_back(v00);
                mesh->indices.push_back(v10);
                mesh->indices.push_back(v11);

                // Triangle 2
                mesh->indices.push_back(v00);
                mesh->indices.push_back(v11);
                mesh->indices.push_back(v01);
            }
        }

        return mesh;
    }

} // namespace base_math

#endif // __nurbs_core_h__
