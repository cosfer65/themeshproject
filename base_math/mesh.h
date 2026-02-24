#ifndef __mesh_h__
#define __mesh_h__

/// @file mesh.h
/// @brief Core half-edge-like mesh data structures and utilities for geometric processing.
///
/// This header defines a lightweight half-edge style mesh representation, curvature
/// data attached to vertices, and helpers for loading simple mesh files and
/// constructing basic primitive meshes.

// might lead to conflict with std::min and std::max in some environments, so we undefine them here to avoid issues
#undef min
#undef max

#include <vector>
#include <list>
#include <map>
#include <set>
#include <fstream>

#include "vector.h"
#include "string_utils.h"

namespace base_math {
    /// @brief Directed edge between two vertices in a 1D chain used for vertex neighborhood ordering.
    ///
    /// `chainLink` is not a topological half-edge; it is a lightweight helper type
    /// for building and ordering the 1-ring of neighbors around a vertex.
    struct chainLink {
        /// @brief Index of the source vertex.
        size_t source;
        /// @brief Index of the target vertex.
        size_t target;
        /// @brief Marks whether this link has already been used when constructing a chain.
        bool visited = false;

        /// @brief Construct a chain edge from @p s to @p t.
        chainLink(size_t s, size_t t) : source(s), target(t) {}
    };

    /// @brief Check whether a set of chain links forms a closed chain.
    ///
    /// A chain is considered *closed* if:
    /// - It is non-empty.
    /// - The set of unique sources equals the set of unique targets.
    ///
    /// The order of the links is not inspected here; only vertex incidence is checked.
    ///
    /// @param chain Collection of directed links (edges) forming the chain.
    /// @return `true` if the chain is closed, `false` otherwise.
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

    /// @brief Reorders a set of chain links to form a continuous chain.
    ///
    /// This function attempts to reorder the vector of `chainLink` entries such that
    /// consecutive elements connect head-to-tail, i.e. the target of one link matches
    /// the source of the next (up to wrapping). It uses a simple O(n²) insertion
    /// strategy and the `visited` flag on each link to avoid reusing links.
    ///
    /// @param chain In/out vector of links to order. The vector is modified in-place.
    /// @return Always returns `true` once the algorithm terminates.
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

    /// @brief Classification of vertices according to Gaussian curvature sign/type.
    ///
    /// This label is typically derived from the per-vertex curvature values stored
    /// inside `curvatureData`.
    enum vertexLabel {
        VERTEX_TYPE_UNKNOWN = 0,   ///< Vertex type is unknown.
        VERTEX_TYPE_ELLIPTIC,      ///< Elliptic (Sphere-like).
        VERTEX_TYPE_PARABOLIC,     ///< Parabolic (Cylinder-like).
        VERTEX_TYPE_PLANE,         ///< Plane.
        VERTEX_TYPE_HYPERBOLIC     ///< Hyperbolic (Saddle-like).
    };

    /// @brief Container for curvature-related data at a single mesh vertex.
    ///
    /// Stores principal curvatures, principal directions, derived scalar measures
    /// (mean, Gaussian, absolute values) and a categorical label for the vertex.
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
        /// Label indicating the type of vertex based on curvature.
        vertexLabel label;

        /// @brief Default-initialize all curvature values to zero and label to unknown.
        curvatureData() : k_min(0), k_max(0), meanCurvature(0), gaussCurvature(0), absKmin(0), absKmax(0),
            absMeanCurvature(0), absGaussCurvature(0), signGauss(0), signMean(0), label(VERTEX_TYPE_UNKNOWN) {
        }
    };

    template <typename T>
    struct meshFace;

    template <typename T>
    struct meshEdge;

    /// @brief Vertex of a polygonal mesh with adjacency and curvature information.
    ///
    /// Each `meshVertex` stores:
    /// - A unique identifier (`id`).
    /// - 3D position and surface normal.
    /// - Incident faces and half-edges (incoming/outgoing).
    /// - Optional Voronoi area and curvature descriptors.
    ///
    /// The adjacency collections are managed by the `mesh` class as faces and edges are added.
    template<typename T>
    struct meshVertex {
        /// @brief Unique vertex identifier (index).
        size_t id;

        /// @brief Cached 1-ring neighbors (may be optional / partially used).
        std::vector<meshVertex<T>*> neighbourVertices;
        /// @brief Incoming directed edges ending at this vertex.
        std::set<meshEdge<T>*> incomingEdges;
        /// @brief Outgoing directed edges starting at this vertex.
        std::set<meshEdge<T>*> outgoingEdges;
        /// @brief Faces that are incident to this vertex.
        std::list<meshFace<T>*> adjacentFaces;

        /// @brief 3D position of the vertex.
        basevector<T, 3> position;
        /// @brief Vertex normal, typically averaged from incident face normals.
        basevector<T, 3> normal;
        /// @brief Mixed (Voronoi) area associated with this vertex.
        T voronoiArea;

        /// @brief Differential-geometric curvature info at this vertex.
        curvatureData<T> curvature_info;

        /// @brief Construct a vertex with id and position.
        meshVertex(size_t id, const basevector<T, 3>& position) : id(id), position(position) {}
        /// @brief Construct a vertex with id, position and normal.
        meshVertex(size_t id, const basevector<T, 3>& position, const basevector<T, 3>& normal) : id(id), position(position), normal(normal) {}
        /// @brief Construct a vertex with id and explicit coordinates.
        meshVertex(size_t id, T x, T y, T z) : id(id), position(x, y, z) {}
    };

    /// @brief Compact identifier for a directed mesh edge given by (source, target) vertex indices.
    typedef std::pair<size_t, size_t> edge_descriptor;

    /// @brief Helper to construct an `edge_descriptor` from two vertex indices.
    ///
    /// @param _s Source vertex index.
    /// @param _t Target vertex index.
    /// @return A pair `( _s, _t )`.
    inline edge_descriptor edge_desc(size_t _s, size_t _t) {
        return edge_descriptor(_s, _t);
    }

    /// @brief Directed edge in the mesh (similar to a half-edge but minimal).
    ///
    /// Each edge stores:
    /// - Source and target vertices.
    /// - The face on its left-hand side (`adjacentFace`).
    /// - The opposite directed edge (if it exists).
    /// - Links to the previous and next edges in the face loop.
    template <typename T>
    struct meshEdge {
        /// @brief Source vertex from which this directed edge starts.
        meshVertex<T>* source;
        /// @brief Target vertex at which this directed edge ends.
        meshVertex<T>* target;
        /// @brief Face that lies to the left of this directed edge in the mesh embedding.
        meshFace<T>* adjacentFace;
        /// @brief Opposite directed edge sharing the same undirected edge but reversed orientation, if any.
        meshEdge<T>* oppositeEdge;
        /// @brief Next edge in the oriented boundary cycle of the incident face.
        meshEdge<T>* nextEdge;
        /// @brief Previous edge in the oriented boundary cycle of the incident face.
        meshEdge<T>* prevEdge;

        /// @brief Construct an edge from source @p s to target @p t.
        ///
        /// The face and connectivity pointers are initialized to `nullptr` and
        /// later populated by `mesh::build_face_edges`.
        meshEdge(meshVertex<T>* s, meshVertex<T>* t) : source(s), target(t), adjacentFace(nullptr), oppositeEdge(nullptr), nextEdge(nullptr), prevEdge(nullptr) {}
    };

    /// @brief Polygonal face of the mesh (triangles or quads are assumed here).
    ///
    /// The face stores a list of vertex pointers (ordered counterclockwise for
    /// consistent normals), a linked half-edge (`firstEdge`), its geometric center,
    /// area, and outward normal vector.
    template <typename T>
    struct meshFace {
        /// @brief Unique face identifier.
        size_t id;
        /// @brief Vertices of this face in boundary order.
        std::vector<meshVertex<T>*> vertices;
        /// @brief One of the edges bordering this face; forms a cycle via `nextEdge`.
        meshEdge<T>* firstEdge;
        /// @brief Normal vector of the face (unit-length when computed).
        basevector<T, 3> normal;
        /// @brief Geometric center (centroid) of the face.
        basevector<T, 3> center;
        /// @brief Surface area of the face.
        T area;

        /// @brief Construct a face with id and vertex list.
        ///
        /// @param _id Unique face identifier.
        /// @param vertices Ordered vertices defining this polygon.
        meshFace(size_t _id, const std::vector<meshVertex<T>*>& vertices) : id(_id), vertices(vertices), firstEdge(nullptr) {}

        /// @brief Number of vertices in this face.
        size_t vertexCount() const { return vertices.size(); }

        /// @brief Compute and cache the face normal from the first three vertices.
        ///
        /// The vertices are assumed to be non-degenerate and ordered such that
        /// the right-hand rule produces an outward normal.
        void computeNormal() {
            if (vertices.size() < 3) return;
            basevector<T, 3>& v0 = vertices[0]->position;
            basevector<T, 3>& v1 = vertices[1]->position;
            basevector<T, 3>& v2 = vertices[2]->position;
            normal = (v1 - v0).cross(v2 - v0).normalize();
        }

        /// @brief Compute both face area and centroid.
        ///
        /// For triangles, a standard triangle area and centroid are used.
        /// For quads, the face is treated as two triangles for area computation
        /// and the centroid is taken as the average of all four vertices.
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

        /// @brief Get the two vertices adjacent to the vertex with id @p v_id on this face.
        ///
        /// @param v_id Id of the vertex on this face.
        /// @return Pair `(prev, next)` of neighboring vertices in boundary order,
        ///         or `(nullptr, nullptr)` if @p v_id is not part of the face.
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

    /// @brief Half-edge like mesh container with basic geometry utilities.
    ///
    /// The mesh stores:
    /// - Vertices (`meshVertex<T>`).
    /// - Directed edges (`meshEdge<T>`) keyed by `(source, target)`.
    /// - Faces (`meshFace<T>`).
    ///
    /// It provides:
    /// - Construction and cleanup of topology.
    /// - Triangle/quad face insertion and quad splitting.
    /// - Computation of face and vertex normals, vertex Voronoi areas.
    /// - Mesh recentering and axis-aligned bounding box extraction.
    template <typename T>
    class mesh {
        /// @brief Next id that will be assigned to newly created faces.
        size_t lastFaceID = 0;
        // average_edge_length: is the mean length of edges, useful for scale - dependent algorithms
        // (e.g., curvature estimation, smoothing step sizes, etc.).
        /// @brief Mean length of all edges in the mesh.
        T average_edge_length = T(0);
        /// @brief Sum of all edge lengths; used to derive `average_edge_length`.
        T total_edge_length = T(0);
        /// @brief Flag indicating whether curvature information has been computed.
        bool m_curvatures_computed = false;

        /// @brief Check whether an edge (v1, v2) already exists and might require a face flip.
        ///
        /// Used when inserting new faces to enforce consistent orientation.
        ///
        /// @param v1 Id of the first vertex.
        /// @param v2 Id of the second vertex.
        /// @return `true` if such an edge is already present, `false` otherwise.
        bool checkFlipEdge(size_t v1, size_t v2) {
            auto it = edges.find(edge_desc(v1, v2));
            if (it != edges.end()) {
                return true;
            }
            return false;
        }

        /// @brief Build half-edge connectivity for a given face.
        ///
        /// This:
        /// - Creates one directed edge per boundary segment.
        /// - Links edges cyclically via `nextEdge` and `prevEdge`.
        /// - Sets `adjacentFace` on edges.
        /// - Binds opposite edges when they already exist in the mesh.
        /// - Updates per-vertex incoming/outgoing edge sets and incident face lists.
        ///
        /// @param face Face whose boundary should be turned into linked edges.
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
        /// @brief All vertices in the mesh, indexed by their id.
        std::map<size_t, meshVertex<T>*> vertices;
        /// @brief All directed edges, indexed by (source, target) descriptor.
        std::map<edge_descriptor, meshEdge<T>*> edges;
        /// @brief All faces, indexed by face id.
        std::map<size_t, meshFace<T>*> faces;

        /// @brief Construct an empty mesh.
        mesh() {}

        /// @brief Destructor that releases all heap-allocated mesh elements.
        ~mesh() { cleanUp(); }

        /// @brief Delete all vertices, edges, and faces, and clear the containers.
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

        /// @brief Insert an externally allocated vertex into the mesh.
        ///
        /// Ownership of @p vertex is transferred to the mesh.
        ///
        /// @param vertex Pointer to vertex to insert.
        void addVertex(meshVertex<T>* vertex) { vertices[vertex->id] = const_cast<meshVertex<T>*>(vertex); }

        /// @brief Create and insert a vertex with a position into the mesh.
        ///
        /// @param id Unique vertex identifier.
        /// @param position Position in 3D.
        void addVertex(size_t id, const basevector<T, 3>& position) { vertices[id] = new meshVertex<T>(id, position); }

        /// @brief Read-only accessor for the curvature-computed flag.
        bool curvatures_computed() const {
            return m_curvatures_computed;
        }

        /// @brief Mutable accessor for the curvature-computed flag.
        ///
        /// This allows client code to set the flag without exposing the member directly.
        bool& curvatures_computed() {
            return m_curvatures_computed;
        }

        /// @brief Get the mean edge length over all edges in the mesh.
        T getAverageEdgeLength() const { return average_edge_length; }

        /// @brief Add a triangular face defined by three vertex ids.
        ///
        /// The method:
        /// - Optionally flips the orientation when shared edges indicate inconsistent winding.
        /// - Constructs the `meshFace`.
        /// - Builds and links its edges.
        ///
        /// @param v1 Id of the first vertex.
        /// @param v2 Id of the second vertex.
        /// @param v3 Id of the third vertex.
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

        /// @brief Add a quadrilateral face defined by four vertex ids.
        ///
        /// Similar to the triangle version but accepts four vertices. The orientation
        /// may be flipped if a conflicting edge is detected.
        ///
        /// @param v1 Id of the first vertex.
        /// @param v2 Id of the second vertex.
        /// @param v3 Id of the third vertex.
        /// @param v4 Id of the fourth vertex.
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

        /// @brief Remove the face with id @p id and all its incident directed edges.
        ///
        /// This:
        /// - Deletes all edges bordering the face.
        /// - Unlinks opposite edges and removes adjacency from vertices.
        /// - Deletes the face object itself.
        ///
        /// @param id Face identifier to be removed.
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

        /// @brief Split a quad face into two triangles along its shorter diagonal.
        ///
        /// If the face is not a quad, the function does nothing.
        ///
        /// The original quad is removed and replaced by two triangular faces. The
        /// choice of diagonal is made by comparing the lengths of v0-v2 and v1-v3.
        ///
        /// @param face Quad face to be split.
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

        /// @brief Convert all quad faces into pairs of triangles.
        ///
        /// This iterates over the current set of faces, collects quads, and
        /// calls `splitFace` on each of them.
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

        /// @brief Compute face normals for all faces in the mesh.
        void computeFaceNormals() {
            for (auto& f_pair : faces) {
                meshFace<T>* f = f_pair.second;
                f->computeNormal();
            }
        }

        /// @brief Compute face normals, centers and areas, and update average edge length.
        ///
        /// This calls `computeNormal()` and `computeAreaandCenter()` on each face and
        /// then computes `average_edge_length` from `total_edge_length` and the
        /// number of edges.
        void computeFaceProperties() {
            for (auto& f_pair : faces) {
                meshFace<T>* f = f_pair.second;
                f->computeNormal();
                f->computeAreaandCenter();
            }
            average_edge_length = total_edge_length / edges.size();
        }

        /// @brief Compute per-vertex normals as the normalized sum of incident face normals.
        ///
        /// Assumes that `computeFaceNormals()` (or `computeFaceProperties()`) has
        /// been called beforehand to populate the face normals.
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
                    v->normal += f->normal; // Add face normal to vertex normal
                }
            }
            // Normalize vertex normals
            for (auto v_pair : vertices) {
                meshVertex<T>* v = v_pair.second;
                v->normal = v->normal.normalize(); // Normalize the accumulated normal
            }
        }

        /// @brief Translate the mesh so that its centroid is moved to the origin.
        ///
        /// The centroid is computed as the arithmetic mean of all vertex coordinates.
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

        /// @brief Compute and store the mixed (Voronoi) area for each vertex.
        ///
        /// The computation uses the classical cotangent weighting scheme for non-obtuse
        /// triangles and a quarter-area rule for obtuse ones, following standard
        /// discrete differential geometry practice.
        ///
        /// The result is stored in `meshVertex::voronoiArea`.
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

        /// @brief Get the ordered 1-ring neighbor ids around a vertex.
        ///
        /// The ordering is derived from incident faces and a chain-link ordering
        /// heuristic so that consecutive neighbors approximate a circular ordering.
        ///
        /// The final list generally starts from an arbitrary neighbor and may repeat
        /// the first vertex at the end to close the loop.
        ///
        /// @param v_id Id of the central vertex.
        /// @return Vector of neighboring vertex ids in ring order.
        std::vector<size_t> getVertexNeighbours(size_t v_id) {
            meshVertex<T>& v = *(vertices[v_id]);
            std::vector<chainLink> ring_edges;
            ring_edges.reserve(v.adjacentFaces.size());
            for (meshFace<T>* f : v.adjacentFaces) {
                auto ov = f->getIncidentVertices(v_id);
                ring_edges.emplace_back(ov.first->id, ov.second->id);
            }
            orderChain(ring_edges);
            std::vector<size_t> result;
            result.reserve(ring_edges.size());
            size_t last_second_vertex = 0;
            for (const auto& e : ring_edges) {
                result.push_back(e.source);
                last_second_vertex = e.target;
            }
            if (last_second_vertex != result.front()) {
                result.push_back(last_second_vertex); // Add the last vertex to close the loop if it's not already the first vertex
            }
            return result;
        }

        /// @brief Compute the axis-aligned bounding box (AABB) of the mesh.
        ///
        /// The returned 2x3 matrix encodes:
        /// - Row 0: (min_x, min_y, min_z)
        /// - Row 1: (max_x, max_y, max_z)
        ///
        /// @return Bounding box as a `basematrix<T, 2, 3>`. If the mesh is empty,
        ///         an empty/default matrix is returned.
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

        /// @brief Translate all vertices of the mesh by a fixed offset.
        ///
        /// @param offset 3D vector added to each vertex position.
        void translate(const basevector<T, 3>& offset) {
            for (auto& v_pair : vertices) {
                basevector<T, 3>& coords = v_pair.second->position;
                coords = coords + offset;
            }
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////
    /**
  * @brief Container for a lightweight mesh definition used for file I/O.
  *
  * The `vertices` vector stores tuples of (id, x, y, z).
  * The `faces` vector stores tuples of (id, v1, v2, v3).
  *
  * @tparam T numeric type for coordinates in the in-memory representation.
  */
    template <typename T>
    struct mesh_def {
        /// @brief List of vertices as (id, x, y, z).
        std::vector<std::tuple<size_t, T, T, T>> vertices;
        /// @brief List of faces as (id, v1, v2, v3).
        std::vector<std::tuple<size_t, size_t, size_t, size_t>> faces;
    };

    /**
     * @brief Read a simple binary mesh file into a `mesh_def` structure.
     *
     * File format (binary):
     *  - size_t vertex_count
     *  - repeated for each vertex: size_t id; float x; float y; float z;
     *  - size_t face_count
     *  - repeated for each face: size_t id; size_t v1; size_t v2; size_t v3;
     *
     * Note: the file is written with 'float' for vertex coordinates; this
     * function converts to template type `T`.
     *
     * @param fnm Path to the binary mesh file.
     * @return Pointer to a newly allocated `mesh_def<T>` on success, or nullptr on failure.
    */
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

    /// @brief Build a `mesh<T>` instance from a lightweight `mesh_def<T>` description.
    ///
    /// This:
    /// - Creates all vertices.
    /// - Adds all faces (triangles).
    /// - Recenters the mesh.
    /// - Computes face properties and vertex normals.
    ///
    /// @tparam T Coordinate type used in the in-memory mesh.
    /// @param _mesh_def Input mesh definition with vertex and face lists.
    /// @return Newly allocated `mesh<T>` instance.
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
        _mesh->computeFaceProperties();
        _mesh->computeVertexNormals();

        return _mesh;
    }

    /// @brief Load a `.prim` mesh file and convert it into a `mesh<T>` instance.
    ///
    /// The file is first parsed into a `mesh_def<T>` via `read_mesh_file`, then
    /// transformed into a fully connected mesh by `parse_mesh`.
    ///
    /// @tparam T Coordinate type used in the in-memory mesh.
    /// @tparam Tinput Floating-point type used in the file (usually `float`).
    /// @param fnm Path to the `.prim` file.
    /// @return Newly allocated `mesh<T>` on success, or `nullptr` on failure.
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

    /// @brief Create a unit cone mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_cone() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_cone.prim");
    }

    /// @brief Create a unit cube mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_cube() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_cube.prim");
    }

    /// @brief Create a unit cylinder mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_cylinder() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_cylinder.prim");
    }

    /// @brief Create a unit dodecahedron mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_dodecahedron() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_dodecahedron.prim");
    }

    /// @brief Create a unit icosahedron mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_icosahedron() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_icosahedron.prim");
    }

    /// @brief Create a unit octahedron mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_octa() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_octa.prim");
    }

    /// @brief Create a unit pentagonal solid mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_penta() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_penta.prim");
    }

    /// @brief Create a unit plane mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_plane() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_plane.prim");
    }

    /// @brief Create a unit sphere mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_sphere() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_sphere.prim");
    }

    /// @brief Create a unit tetrahedron mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_tetra() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_tetra.prim");
    }

    /// @brief Create a unit torus mesh from the bundled `.prim` resource.
    template <typename T, typename Tinput>
    mesh<T>* create_unit_torus() {
        return create_from_mesh_file<T, Tinput>("resources\\models\\unit_torus.prim");
    }
} // namespace base_math

#endif // __mesh_h__
