#ifndef __mesh_h__
#define __mesh_h__

// might lead to conflict with std::min and std::max in some environments, so we undefine them here to avoid issues
#undef min
#undef max

#include <vector>
#include <list>
#include <map>
#include <set>
#include <fstream>

#include "vector.h"
#include "geometry.h"
#include "string_utils.h"

namespace base_math {
    template <typename T>
    struct TessellatedMesh {
        std::vector<basevector<T, 3>> vertices;
        std::vector<basevector<T, 3>> normals;
        std::vector<int> indices;
    };

    struct chainLink {
        size_t source;
        size_t target;
        bool visited = false;
        chainLink(size_t s, size_t t) : source(s), target(t) {}
    };

    inline bool isChainClosed(const std::vector<chainLink>& chain) {
        if (chain.empty()) return false;
        std::set<size_t> sources;
        std::set<size_t> targets;
        for (const auto& edge : chain) {
            sources.insert(edge.source);
            targets.insert(edge.target);
        }
        if (sources.size() != targets.size()) return false; // must have same number of unique sources and targets
        if (sources != targets) return false; // must have same unique vertices as sources and targets
        return true; // all checks passed, it's a closed chain
    }

    inline bool orderChain(std::vector<chainLink>& chain) {
        const size_t many = chain.size();
        while (true) {
            bool changed = false;
            size_t requested_target = chain[0].source;
            for (size_t cur = 1; cur < many; ++cur) {
                if (chain[cur].target == requested_target && !chain[cur].visited) {
                    for (size_t i = cur; i > 0; --i) {
                        std::swap(chain[i], chain[i - 1]);
                    }
                    requested_target = chain[0].source;
                    chain[0].visited = true;
                    changed = true;
                }
            }
            if (!changed) {
                break; // no changes made, we are done
            }
        }
        return true; // successfully ordered into a chain
    }

    template<typename T>
    struct curvatureData {
        /// Principal curvatures at the vertex
        T k_min;
        T k_max;
        /// Principal directions corresponding to 'principal_curvatures'.
        basevector<T, 3> k_min_dir;
        basevector<T, 3> k_max_dir;
        /// Mean curvature H = (k_min + k_max) / 2.
        T meanCurvature;
        basevector<T, 3> meanCurvatureDir; // this can be computed as the normalized sum of the principal curvature directions weighted by their respective curvatures, i.e., (k_min * k_min_dir + k_max * k_max_dir) / (k_min + k_max), and it points in the direction of the surface normal if the mean curvature is positive (convex) and in the opposite direction if the mean curvature is negative (concave).
        /// Gaussian curvature K = k_min * k_max.
        T gaussCurvature;

        // cached derived quantities for convenience
        /// Absolute value of the minimum principal curvature |k_min|.
        T absKmin;
        /// Absolute value of the maximum principal curvature |k_max|.
        T absKmax;
        /// absolute value of Mean curvature.
        T absMeanCurvature;
        /// absolute value of Gaussian curvature.
        T absGaussCurvature;

        /// Sign of Gaussian curvature K (e.g., -1 for saddle, 0 for flat, 1 for elliptic).
        int signGauss;
        /// Sign of mean curvature H (e.g., -1 for concave, 0 for minimal, 1 for convex).
        int signMean;

        /// @brief Default-initialize all curvature values to zero and label to unknown.
        curvatureData() : k_min(0), k_max(0), meanCurvature(0), gaussCurvature(0), absKmin(0), absKmax(0),
            absMeanCurvature(0), absGaussCurvature(0), signGauss(0), signMean(0) {}
        void reset() {
            k_min = 0; k_max = 0; meanCurvature = 0; gaussCurvature = 0;
            absKmin = 0; absKmax = 0; absMeanCurvature = 0; absGaussCurvature = 0;
            signGauss = 0; signMean = 0;
        }
    };

    template <typename T>
    struct meshFace;

    template <typename T>
    struct meshEdge;

    template<typename T>
    struct meshVertex {
        size_t id;

        std::vector<meshVertex<T>*> neighbourVertices;
        std::list<meshFace<T>*> adjacentFaces;

        std::set<meshEdge<T>*> incomingEdges;
        std::set<meshEdge<T>*> outgoingEdges;

        basevector<T, 3> position;
        basevector<T, 3> normal;

        bool is_boundary = false; // this can be used to mark boundary vertices, which is useful for curvature estimation and other algorithms that need to treat boundaries differently

        meshVertex(size_t id, const basevector<T, 3>& position) : id(id), position(position) {}
        meshVertex(size_t id, const basevector<T, 3>& position, const basevector<T, 3>& normal) : id(id), position(position), normal(normal) {}
        meshVertex(size_t id, T x, T y, T z) : id(id), position(x, y, z) {}

        T angle_sum; // this can be used to store the angle sum around this vertex, which is useful for curvature estimation
        T area_sum; // this can be used to store the total area of faces adjacent to this vertex, which is also useful for curvature estimation
        /// @brief Mixed (Voronoi) area associated with this vertex.
        T voronoiArea;
        /// @brief Differential-geometric curvature info at this vertex.
        curvatureData<T> curvature_info;
        T curvature_map_value; // this can be used to store a single scalar curvature value for visualization (e.g. Gaussian curvature)

        void compute_area_sum() {
            area_sum = 0;
            for (meshFace<T>* face : adjacentFaces) {
                area_sum += face->area;
            }
        }
        void compute_angle_sum() {
            angle_sum = 0;
            for (meshFace<T>* face : adjacentFaces) {
                std::pair<meshVertex<T>*, meshVertex<T>*> other = face->getIncidentVertices(id);
                if (other.first && other.second) {
                    T angle = corner_angle(position, other.first->position, other.second->position);
                    angle_sum += angle;
                }
            }
        }
    };

    typedef std::pair<size_t, size_t> edge_descriptor;

    inline edge_descriptor edge_desc(size_t _s, size_t _t) {
        return edge_descriptor(_s, _t);
    }

    template <typename T>
    struct meshEdge {
        meshVertex<T>* source;
        meshVertex<T>* target;
        meshFace<T>* adjacentFace;
        meshEdge<T>* oppositeEdge;
        meshEdge<T>* nextEdge;
        meshEdge<T>* prevEdge;

        meshEdge(meshVertex<T>* s, meshVertex<T>* t) : source(s), target(t), adjacentFace(nullptr), oppositeEdge(nullptr), nextEdge(nullptr), prevEdge(nullptr) {}
    };

    template <typename T>
    struct meshFace {
        size_t id;
        std::vector<meshVertex<T>*> vertices;
        std::vector<T> corner_angles; // this can be used to store the corner angles at each vertex of the face, which is useful for curvature estimation
        meshEdge<T>* firstEdge;
        basevector<T, 3> normal;
        basevector<T, 3> normal_raw; // this can be used to store the unnormalized normal vector (cross product of two edge vectors) which is useful for curvature estimation
        basevector<T, 3> center;
        T area;

        meshFace(size_t _id, const std::vector<meshVertex<T>*>& vertices) : id(_id), vertices(vertices), firstEdge(nullptr) {}

        size_t vertexCount() const { return vertices.size(); }

        void computeNormal() {
            if (vertices.size() < 3) return;
            basevector<T, 3>& v0 = vertices[0]->position;
            basevector<T, 3>& v1 = vertices[1]->position;
            basevector<T, 3>& v2 = vertices[2]->position;
            normal_raw = (v1 - v0).cross(v2 - v0);
            normal = normalize(normal_raw);
        }

        void computeAreaandCenter() {
            if (vertices.size() < 3) return;
            basevector<T, 3>& v0 = vertices[0]->position;
            basevector<T, 3>& v1 = vertices[1]->position;
            basevector<T, 3>& v2 = vertices[2]->position;
            T div = 3.;
            basevector<T, 3> _center = (v0 + v1 + v2);
            if (vertices.size() == 4) {
                basevector<T, 3>& v3 = vertices[3]->position;
                _center += v3;
                div = 4.;
            }
            center = _center / div;

            area = (v1 - v0).cross(v2 - v0).length() * 0.5;
            if (vertices.size() == 4) {
                basevector<T, 3>& v3 = vertices[3]->position;
                area += (v2 - v0).cross(v3 - v0).length() * 0.5;
            }
        }

        void computeCornerAngles() {
            if (vertices.size() < 3) return;
            corner_angles.clear();
            size_t vCount = vertices.size();
            for (size_t i = 0; i < vCount; ++i) {
                basevector<T, 3>& v_prev = vertices[(i + vCount - 1) % vCount]->position;
                basevector<T, 3>& v_curr = vertices[i]->position;
                basevector<T, 3>& v_next = vertices[(i + 1) % vCount]->position;
                T angle = corner_angle(v_curr, v_prev, v_next);
                corner_angles.push_back(angle);
            }
        }

        std::pair<meshVertex<T>*, meshVertex<T>*> getIncidentVertices(size_t v_id) {
            size_t vCount = vertices.size();
            for (size_t i = 0; i < vCount; ++i) {
                meshVertex<T>* v = vertices[i];
                if (v->id == v_id) {
                    return std::pair<meshVertex<T>*, meshVertex<T>*>(vertices[(i + vCount - 1) % vCount], vertices[(i + 1) % vCount]);
                }
            }
            return std::pair<meshVertex<T>*, meshVertex<T>*>(nullptr, nullptr);
        }
    };

    template <typename T>
    class mesh {
        size_t lastFaceID = 0;
        T average_edge_length = T(0);
        T total_edge_length = T(0);
        bool m_curvatures_computed = false;

        bool checkFlipEdge(size_t v1, size_t v2) {
            auto it = edges.find(edge_desc(v1, v2));
            if (it != edges.end()) {
                return true;
            }
            return false;
        }

        void build_face_edges(meshFace<T>* face) {
            const std::vector<meshVertex<T>*>& verts = face->vertices;
            size_t vCount = verts.size();
            std::vector<meshEdge<T>*> faceEdges(vCount, nullptr);

            for (size_t i = 0; i < vCount; ++i) {
                // previous vertex is (i-1) mod vCount, current vertex is i, next vertex is (i+1) mod vCount to wrap around
                meshVertex<T>* v0 = verts[(i + vCount - 1) % vCount]; // previous vertex
                meshVertex<T>* v1 = verts[i];  // current vertex
                meshVertex<T>* v2 = verts[(i + 1) % vCount];  // next vertex

                v1->adjacentFaces.push_back(face); // add face to current vertex's adjacent faces
                // create edge from v1 to v2
                edge_descriptor ed = edge_desc(v1->id, v2->id);
                meshEdge<T>* e = new meshEdge<T>(v1, v2);
                total_edge_length += (basevector<T, 3>(v1->position - v2->position)).length();

                v1->outgoingEdges.insert(e); // add edge to current vertex's outgoing edges
                v2->incomingEdges.insert(e); // add edge to next vertex's incoming edges

                // check if opposite edge from v2 to v1 already exists
                edge_descriptor ed_opposite = edge_desc(v2->id, v1->id);
                auto it = edges.find(ed_opposite);
                meshEdge<T>* op_Edge = (it != edges.end()) ? it->second : nullptr;
                if (op_Edge != nullptr) {
                    op_Edge->oppositeEdge = e; // set opposite edge's opposite to current edge
                    e->oppositeEdge = op_Edge; // set current edge's opposite to opposite edge
                }

                // set adjacent face for the edge
                e->adjacentFace = face;
                faceEdges[i] = e; // store edge in temporary array

                // add edge to mesh edge map
                edges[ed] = e;
                if (face->firstEdge == nullptr) {
                    face->firstEdge = e; // set the first edge of the face
                }
            }
            // link the edges in a cycle
            for (size_t i = 0; i < vCount; ++i) {
                faceEdges[i]->nextEdge = faceEdges[(i + 1) % vCount];
                faceEdges[i]->prevEdge = faceEdges[(i + vCount - 1) % vCount];
            }
        }
    public:
        std::map<size_t, meshVertex<T>*> vertices;
        std::map<edge_descriptor, meshEdge<T>*> edges;
        std::map<size_t, meshFace<T>*> faces;

        mesh() {}

        ~mesh() { cleanUp(); }

        void cleanUp() {
            for (auto& v_pair : vertices) {
                delete v_pair.second;
            }
            for (auto& e_pair : edges) {
                delete e_pair.second;
            }
            for (auto& f_pair : faces) {
                delete f_pair.second;
            }
            vertices.clear();
            edges.clear();
            faces.clear();
        }

        void flip_face(meshFace<T>* f) {
            std::reverse(f->vertices.begin(), f->vertices.end());
            meshEdge<T>* e = f->firstEdge;
            do {
                std::swap(e->source, e->target);
                e = e->nextEdge;
            } while (e != f->firstEdge);
            f->normal = -f->normal; // flip the normal direction
        }

        void flip_face(size_t id) {
            auto it = faces.find(id);
            if (it != faces.end()) {
                meshFace<T>* f = it->second;
                flip_face(f);
            }
        }

        void flip_all_faces() {
            for (auto& f_pair : faces) {
                flip_face(f_pair.first);
            }
        }

        void addVertex(meshVertex<T>* vertex) { vertices[vertex->id] = const_cast<meshVertex<T>*>(vertex); }

        meshVertex<T>* addVertex(size_t id, const basevector<T, 3>& position) {
            vertices[id] = new meshVertex<T>(id, position);
            return vertices[id];
        }

        bool curvatures_computed() const {
            return m_curvatures_computed;
        }

        bool& curvatures_computed() {
            return m_curvatures_computed;
        }

        T getAverageEdgeLength() const { return average_edge_length; }

        void addFace(size_t v1, size_t v2, size_t v3) {
            bool flip = checkFlipEdge(v1, v2) || checkFlipEdge(v2, v3) || checkFlipEdge(v3, v1);
            if (flip) {
                size_t t = v2;
                v2 = v3;
                v3 = t;
            }

            meshFace<T>* f = new meshFace<T>(lastFaceID, { vertices[v1], vertices[v2], vertices[v3] });
            faces[lastFaceID] = f;
            build_face_edges(f);
            ++lastFaceID;
        }

        void addFace(size_t v1, size_t v2, size_t v3, size_t v4) {
            bool flip = checkFlipEdge(v1, v2) || checkFlipEdge(v2, v3) || checkFlipEdge(v3, v4) || checkFlipEdge(v4, v1);
            if (flip) {
                size_t t = v2;
                v2 = v4;
                v4 = t;
            }
            meshFace<T>* f = new meshFace<T>(lastFaceID, { vertices[v1], vertices[v2], vertices[v3], vertices[v4] });
            faces[lastFaceID] = f;
            build_face_edges(f);
            ++lastFaceID;
        }

        void removeFace(size_t id) {
            auto it = faces.find(id);
            if (it != faces.end()) {
                meshFace<T>* f = it->second;
                // Remove edges associated with this face
                meshEdge<T>* e = f->firstEdge;
                size_t vCount = f->vertexCount();
                for (size_t i = 0; i < vCount; ++i) {
                    meshVertex<T>* sourceVertex = e->source;
                    meshVertex<T>* targetVertex = e->target;

                    edge_descriptor ed(sourceVertex->id, targetVertex->id);
                    edges.erase(ed); // remove edge from mesh edge map
                    if (e->oppositeEdge) {
                        e->oppositeEdge->oppositeEdge = nullptr; // disconnect opposite edge
                    }
                    sourceVertex->outgoingEdges.erase(e); // remove edge from source vertex's outgoing edges
                    targetVertex->incomingEdges.erase(e); // remove edge from target vertex's incoming edges
                    sourceVertex->adjacentFaces.remove(f); // remove face from source vertex's adjacent faces
                    targetVertex->adjacentFaces.remove(f); // remove face from target vertex's adjacent faces
                    meshEdge<T>* next = e->nextEdge; // store next edge
                    delete e; // free edge memory
                    e = next; // move to next edge in the face
                }
                faces.erase(it); // remove face from mesh face map
                delete f; // free face memory
            }
        }

        void splitFace(meshFace<T>* face) {
            if (face->vertexCount() != 4) {
                // Only split quads
                return;
            }
            // Add two new triangular faces to replace the quad
            // better select the shorter diagonal to split the quad
            T diag1 = (basevector<T, 3>(face->vertices[2]->position - face->vertices[0]->position)).length();
            T diag2 = (basevector<T, 3>(face->vertices[3]->position - face->vertices[1]->position)).length();
            // Get vertex IDs to build new faces (face will be deleted in removeFace, so we need to store the IDs before that)
            size_t v1 = face->vertices[0]->id, v2 = face->vertices[1]->id;
            size_t v3 = face->vertices[2]->id, v4 = face->vertices[3]->id;
            removeFace(face->id); // Remove the original quad face from the mesh
            if (diag2 < diag1) {
                // Split along v2-v4
                addFace(v1, v2, v4);
                addFace(v2, v3, v4);
            }
            else {
                // Split along v1-v3
                addFace(v1, v2, v3);
                addFace(v1, v3, v4);
            }
        }

        void breakQuads() {
            std::vector<meshFace<T>*> quads;
            for (auto f_pair : faces) {
                meshFace<T>* f = f_pair.second;
                if (f->vertexCount() == 4) {
                    quads.push_back(f);
                }
            }
            size_t next_face_id = faces.size() + 1;

            for (meshFace<T>* f : quads) {
                splitFace(f);
            }
        }

        void computeFaceNormals() {
            for (auto& f_pair : faces) {
                meshFace<T>* f = f_pair.second;
                f->computeNormal();
            }
        }

        void computeFaceProperties() {
            for (auto& f_pair : faces) {
                meshFace<T>* f = f_pair.second;
                f->computeNormal();
                f->computeAreaandCenter();
                f->computeCornerAngles();
            }
            average_edge_length = total_edge_length / edges.size();
        }

        void computeVertexNormals() {
            // Initialize all vertex normals to zero
            for (auto v_pair : vertices) {
                meshVertex<T>* v = v_pair.second;
                v->normal = basevector<T, 3>(0, 0, 0);
            }
            // Accumulate face normals into vertex normals
            // For each face, add its normal to each of its vertices' normals
            // we assume that the face normals have already been computed and stored in the face's `normal` member
            for (auto f_pair : faces) {
                meshFace<T>* f = f_pair.second;
                for (meshVertex<T>* v : f->vertices) {
                    // Use area‑weighted face normals
                    v->normal += f->normal*f->area; // Add face normal to vertex normal
                }
            }
            // Normalize vertex normals
            for (auto v_pair : vertices) {
                meshVertex<T>* v = v_pair.second;
                v->normal = v->normal.normalize(); // Normalize the accumulated normal
            }
        }

        void computeVerticesVoronoiArea() {
            for (auto v_pair : vertices) {
                meshVertex<T>* v = v_pair.second;
                T area_mixed = 0.0;
                basevector<T, 3>& P = v->position;
                // Iterate over all faces adjacent to this vertex
                for (meshFace<T>* f : v->adjacentFaces) {
                    // Get the other two vertices of the face (assuming triangular faces)
                    meshVertex<T>* v1 = nullptr;
                    meshVertex<T>* v2 = nullptr;
                    std::pair< meshVertex<T>*, meshVertex<T>*> other = f->getIncidentVertices(v->id);
                    v1 = other.first;
                    v2 = other.second;
                    if (!v1 || !v2) continue; // Skip if we couldn't find two other vertices
                    // Edge vectors
                    basevector<T, 3> PR = v2->position - P;
                    basevector<T, 3> PQ = v1->position - P;
                    basevector<T, 3> QR = v2->position - v1->position;
                    // Total area of this triangle
                    T area_T = T(0.5) * PQ.cross(PR).length();
                    // Check if triangle is obtuse
                    bool is_obtuse = PQ.dot(PR) < 0 || (-PQ).dot(QR) < 0 || (-PR).dot(-QR) < 0;
                    if (!is_obtuse) {
                        // Voronoi Area using cotangent weights
                        T cot_alpha = PQ.dot(PR) / PQ.cross(PR).length();
                        T cot_beta = (-PQ).dot(QR) / (-PQ).cross(QR).length();
                        area_mixed += (cot_alpha * PR.length() * PR.length() + cot_beta * PQ.length() * PQ.length()) / T(8.0);
                    }
                    else {
                        // Angle at Q or R is obtuse
                        area_mixed += area_T / T(4.0);
                    }
                }
                // Store the computed Voronoi area in the vertex's normal.w() component as an example (you can choose a different member to store it)
                v->voronoiArea = area_mixed; // Store Voronoi
            }
        }

        T computeVoronoiArea(int vertex_idx, meshFace<T>* face) {
            int i;
            for (i = 0; i < 3; ++i) {
                if (face->vertices[i]->id == vertex_idx) {
                    // Found the vertex in this face
                    break;
                }
            }
            basevector<T, 3> v0 = face->vertices[i]->position;
            basevector<T, 3> v1 = face->vertices[(i + 1) % 3]->position;
            basevector<T, 3> v2 = face->vertices[(i + 2) % 3]->position;

            // Check if triangle is obtuse
            T a2 = (v1 - v2).length_sq();
            T b2 = (v0 - v2).length_sq();
            T c2 = (v0 - v1).length_sq();

            if (a2 > b2 + c2) { // angle at v0 is obtuse
                return face->area / T(2.0);
            }
            else if (b2 > a2 + c2) { // angle at v1 is obtuse
                return face->area / T(2.0);
            }
            else if (c2 > a2 + b2) { // angle at v2 is obtuse
                return face->area / T(2.0);
            }
            // Voronoi area formula for non-obtuse triangles
            return face->area / T(3.0);
        }

        void computeVertexVoronoiArea(meshVertex<T>* v) {
            v->voronoiArea = T(0); // Initialize Voronoi area to zero
            // Iterate over all faces adjacent to this vertex
            for (meshFace<T>* f : v->adjacentFaces) {
                v->voronoiArea += computeVoronoiArea(int(v->id), f);
            }
        }

        std::vector<size_t> getVertexNeighbours(meshVertex<T>& v) {
            v.neighbourVertices.clear();
            std::vector<chainLink> ring_edges;
            ring_edges.reserve(v.adjacentFaces.size());
            for (meshFace<T>* f : v.adjacentFaces) {
                auto ov = f->getIncidentVertices(v.id);
                ring_edges.emplace_back(ov.first->id, ov.second->id);
            }
            orderChain(ring_edges);
            std::vector<size_t> result;
            result.reserve(ring_edges.size());
            size_t last_second_vertex = 0;
            for (const auto& e : ring_edges) {
                result.push_back(e.source);
                last_second_vertex = e.target;
                v.neighbourVertices.push_back(vertices[e.source]); // Add the source vertex to the neighbor list
            }
            if (last_second_vertex != result.front()) {
                result.push_back(last_second_vertex); // Add the last vertex to close the loop if it's not already the first vertex
                v.neighbourVertices.push_back(vertices[last_second_vertex]); // Add the last vertex to the neighbor list
            }
            return result;
        }

        void mark_boundary_vertices() {
            for (auto& edge : edges) {
                if (edge.second->oppositeEdge == nullptr) {
                    edge.second->source->is_boundary = true;
                    edge.second->target->is_boundary = true;
                }
            }
        }

        void computeVertexProperties() {
            // T angle_sum; // this can be used to store the angle sum around this vertex, which is useful for curvature estimation
            // T area_sum; // this can be used to store the total area of faces adjacent to this vertex, which is also useful for curvature estimation

            computeVertexNormals();
            mark_boundary_vertices(); // we need this to calculate Voronoi areas correctly and to handle curvature estimation properly at boundaries
            for (auto& v_pair : vertices) {
                meshVertex<T>* v = v_pair.second;
                computeVertexVoronoiArea(v);
                getVertexNeighbours(*v);

                v->compute_area_sum();
                v->compute_angle_sum();

                v->curvature_info.reset(); // reset curvature info to default values (you can implement actual curvature estimation here)
                v->curvature_map_value = 0; // reset curvature map value to default (you can compute an actual curvature value here, e.g. Gaussian curvature)
                v->is_boundary = false;
            }
        }

        void recalculateMesh() {
            computeFaceProperties();
            computeVertexProperties();
        }

        void recenterMesh() {
            T sumx = T(0), sumy = T(0), sumz = T(0);
            for (auto& v_pair : vertices) {
                basevector<T, 3>& vc = v_pair.second->position;
                sumx += vc.x();
                sumy += vc.y();
                sumz += vc.z();
            }
            sumx /= T(vertices.size());
            sumy /= T(vertices.size());
            sumz /= T(vertices.size());
            for (auto& v_pair : vertices) {
                basevector<T, 3>& vc = v_pair.second->position;
                vc.x() -= sumx;
                vc.y() -= sumy;
                vc.z() -= sumz;
            }
        }

        basematrix<T, 2, 3> getBoundingBox() {
            if (vertices.size() == 0) {
                return basematrix<T, 2, 3>(); // empty matrix
            }
            T min_x = std::numeric_limits<T>::max();
            T min_y = std::numeric_limits<T>::max();
            T min_z = std::numeric_limits<T>::max();
            T max_x = std::numeric_limits<T>::lowest();
            T max_y = std::numeric_limits<T>::lowest();
            T max_z = std::numeric_limits<T>::lowest();

            for (const auto& v_pair : vertices) {
                const basevector<T, 3>& coords = v_pair.second->position;
                if (coords.x() < min_x) min_x = coords.x();
                if (coords.y() < min_y) min_y = coords.y();
                if (coords.z() < min_z) min_z = coords.z();
                if (coords.x() > max_x) max_x = coords.x();
                if (coords.y() > max_y) max_y = coords.y();
                if (coords.z() > max_z) max_z = coords.z();
            }
            basematrix<T, 2, 3> bbox({ min_x, min_y, min_z, max_x, max_y, max_z });
            return bbox;
        }

        void translate(const basevector<T, 3>& offset) {
            for (auto& v_pair : vertices) {
                basevector<T, 3>& coords = v_pair.second->position;
                coords = coords + offset;
            }
        }
        void scale(const basevector<T, 3>& factor) {
            for (auto& v_pair : vertices) {
                basevector<T, 3>& coords = v_pair.second->position;
                coords.x() *= factor.x();
                coords.y() *= factor.y();
                coords.z() *= factor.z();
            }
        }
        void scale(T factor) {
            for (auto& v_pair : vertices) {
                basevector<T, 3>& coords = v_pair.second->position;
                coords.x() *= factor;
                coords.y() *= factor;
                coords.z() *= factor;
            }
        }
        void rotate(const basevector<T, 3>& euler_angles) {
            // Convert Euler angles to radians
            T cx = std::cos(euler_angles.x());
            T sx = std::sin(euler_angles.x());
            T cy = std::cos(euler_angles.y());
            T sy = std::sin(euler_angles.y());
            T cz = std::cos(euler_angles.z());
            T sz = std::sin(euler_angles.z());
            // Rotation matrices for each axis
            basematrix<T, 3, 3> Rx({ 1, 0, 0, 0, cx, -sx, 0, sx, cx });
            basematrix<T, 3, 3> Ry({ cy, 0, sy, 0, 1, 0, -sy, 0, cy });
            basematrix<T, 3, 3> Rz({ cz, -sz, 0, sz, cz, 0, 0, 0, 1 });
            // Combined rotation matrix (Rz * Ry * Rx)
            basematrix<T, 3, 3> R = Rz * Ry * Rx;
            for (auto& v_pair : vertices) {
                basevector<T, 3>& coords = v_pair.second->position;
                basevector<T, 3> rotated_coords(
                    R[0] * coords.x() + R[1] * coords.y() + R[2] * coords.z(),
                    R[3] * coords.x() + R[4] * coords.y() + R[5] * coords.z(),
                    R[6] * coords.x() + R[7] * coords.y() + R[8] * coords.z()
                );
                coords = rotated_coords;
            }
        }
        void rotate(T angle_x, T angle_y, T angle_z) {
            rotate(basevector<T, 3>(angle_x, angle_y, angle_z));
        }
        void transform(const basematrix<T, 4, 4>& transform) {
            for (auto& v_pair : vertices) {
                basevector<T, 3>& coords = v_pair.second->position;
                basevector<T, 4> homogenous_coords(coords.x(), coords.y(), coords.z(), T(1));
                basevector<T, 4> transformed_coords = transform * homogenous_coords;
                coords.x() = transformed_coords.x() / transformed_coords.w();
                coords.y() = transformed_coords.y() / transformed_coords.w();
                coords.z() = transformed_coords.z() / transformed_coords.w();
            }
        }
    };

    template <typename T>
    mesh<T>* parse_mesh(const TessellatedMesh<T>& _mesh_def) {
        mesh<T>* _mesh = new mesh<T>;
        // add vertices
        size_t vertex_index = 0;
        for (const basevector<T, 3>&vert : _mesh_def.vertices) {
            meshVertex<T>* v = _mesh->addVertex(vertex_index, vert);
            v->normal = _mesh_def.normals[vertex_index];
            ++vertex_index;
        }

        // add faces
        size_t index_count = _mesh_def.indices.size();
        for (size_t idx = 0; idx < index_count; idx += 3) {
            _mesh->addFace(_mesh_def.indices[idx], _mesh_def.indices[idx + 1], _mesh_def.indices[idx + 2]);
        }
        //_mesh->recenterMesh();
        // _mesh->computeFaceProperties();
        // _mesh->computeVertexNormals(); //tesselated meshes from NURBS should already have good normals, so we can skip this step for now
        _mesh->recalculateMesh();

        return _mesh;
    }

    template <typename T>
    struct mesh_def {
        std::vector<std::tuple<size_t, T, T, T>> vertices;
        std::vector<std::tuple<size_t, size_t, size_t, size_t>> faces;
    };

    template <typename T, typename Tinput>
    struct mesh_def<T>* read_mesh_file(const std::string& fnm) {
        std::ifstream rf(fnm, std::ios::in | std::ios::binary);

        if (!rf)
            return nullptr;
        mesh_def<T>* mesh = new mesh_def<T>;
        size_t d, i;
        rf.read((char*)&d, sizeof(size_t));
        for (i = 0; i < d; ++i) {
            size_t id;
            Tinput x, y, z; // the file is built with 'float'!!!
            rf.read((char*)&id, sizeof(size_t));
            rf.read((char*)&x, sizeof(Tinput));
            rf.read((char*)&y, sizeof(Tinput));
            rf.read((char*)&z, sizeof(Tinput));
            mesh->vertices.push_back(std::tuple<size_t, T, T, T>(id, T(x), T(y), T(z)));
        }
        rf.read((char*)&d, sizeof(size_t));
        for (i = 0; i < d; ++i) {
            size_t id;
            size_t v1, v2, v3;
            rf.read((char*)&id, sizeof(size_t));
            rf.read((char*)&v1, sizeof(size_t));
            rf.read((char*)&v2, sizeof(size_t));
            rf.read((char*)&v3, sizeof(size_t));
            mesh->faces.push_back(std::tuple<size_t, size_t, size_t, size_t>(id, v1, v2, v3));
        }
        return mesh;
    }

    template <typename T>
    mesh<T>* parse_mesh(const mesh_def<T>& _mesh_def) {
        mesh<T>* _mesh = new mesh<T>;
        // add vertices
        for (const std::tuple<size_t, T, T, T>& vert : _mesh_def.vertices) {
            _mesh->addVertex(std::get<0>(vert), basevector<T, 3>(T(std::get<1>(vert)), T(std::get<2>(vert)), T(std::get<3>(vert))));
        }

        // add faces
        for (const std::tuple<size_t, size_t, size_t, size_t>& fc : _mesh_def.faces) {
            _mesh->addFace(std::get<1>(fc), std::get<2>(fc), std::get<3>(fc));
        }
        // mesh->average_edge_length = mesh->total_edge_length / mesh->half_edges.size();
        // mesh->orient_mesh();
        _mesh->recenterMesh();
        //_mesh->computeFaceNormals();
        // _mesh->computeFaceProperties();
        // _mesh->computeVertexNormals();
        _mesh->recalculateMesh();

        return _mesh;
    }

    template <typename T, typename Tinput>
    mesh<T>* create_from_mesh_file(const std::string& fnm) {
        if (file_extension(fnm) != "prim")
            return nullptr;
        struct mesh_def<T>* _mesh = read_mesh_file<T, Tinput>(fnm);
        if (_mesh == nullptr)
            return nullptr;
        mesh<T>* he_mesh = parse_mesh<T>(*_mesh);
        delete _mesh;
        return he_mesh;
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_cone() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_cone.prim");
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_cube() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_cube.prim");
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_cylinder() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_cylinder.prim");
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_dodecahedron() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_dodecahedron.prim");
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_icosahedron() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_icosahedron.prim");
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_octa() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_octa.prim");
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_penta() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_penta.prim");
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_plane() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_plane.prim");
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_sphere() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_sphere.prim");
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_tetra() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_tetra.prim");
    }

    template <typename T, typename Tinput>
    mesh<T>* create_unit_torus() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_torus.prim");
    }
} // namespace base_math

#endif // __mesh_h__
