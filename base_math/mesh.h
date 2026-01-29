#ifndef __mesh_h__
#define __mesh_h__

/**
 * @file mesh.h
 * @brief Lightweight half-edge mesh implementation and helper functions.
 *
 * This header provides a compact half-edge mesh representation suitable for
 * simple geometry processing tasks: storing vertices, triangular faces and
 * half-edges, computing face/vertex normals, centers and areas, reading a
 * simple binary mesh format and constructing predefined unit meshes from
 * resource files.
 *
 * The implementation favors clarity over performance and uses std::map for
 * lookups by id/edge descriptors. Faces are assumed to be triangles.
 */

#undef min
#undef max

#ifdef _DEBUG
#include <iostream>
#endif
#include <fstream>
#include <map>
#include <set>
#include "string_utils.h"
// 
#include "geometry.h"
#include "vector.h"

namespace base_math {
    inline bool order_chain(std::vector<std::pair<size_t, size_t>>& chain) {
        size_t cur;
        size_t many = chain.size();
        for (cur = 1; cur < many; ++cur) {
            size_t i;
            for (i = cur; i < many; ++i) {
                if (chain[i].first == chain[cur - 1].second) {
                    if (i != cur) {
                        auto temp = chain[cur];
                        chain[cur] = chain[i];
                        chain[i] = temp;
                    }
                    break;
                }
            }
        }

        return true;
    }

    enum vertex_label {
        VERTEX_TYPE_UNKNOWN = 0,   ///< Vertex type is unknown.
        VERTEX_TYPE_ELLIPTIC,      ///< Elliptic (Sphere-like).
        VERTEX_TYPE_PARABOLIC,     ///< Parabolic (Cylinder-like).
        VERTEX_TYPE_PLANE,         ///< Plane.
        VERTEX_TYPE_HYPERBOLIC     ///< Hyperbolic (Saddle-like).
    };

    /// Edge key type used to index half-edges in maps (source vertex id, target vertex id).
    typedef std::pair<size_t, size_t> edge_descriptor;

    template <typename T>
    struct face;

    template <typename T>
    struct half_edge;

    /**
  * @brief Vertex stored in the half-edge mesh.
  *
  * Stores a unique identifier, its 3D coordinates, a per-vertex normal and a
  * list of incident face ids. Normals are expected to be unit-length after
  * calling the mesh convenience routines.
  *
  * @tparam T numeric type (float, double, etc.)
  */
    template <typename T>
    struct vertex {
        /// Unique vertex id (user-provided).
        size_t id;

        /// List of incident face ids (each face is identified by its id).
        std::vector<size_t> incident_faces;

        /// 3D coordinates of the vertex.
        basevector<T, 3> coords;

        /// Per-vertex normal (accumulated and normalized by compute_vertex_normals()).
        basevector<T, 3> normal;

        // Voronoi area associated with this vertex
        T voronoi_area;

        /// Maximum principal curvature at this vertex (largest eigenvalue of shape operator).
        T curvature_max;

        /// Minimum principal curvature at this vertex (smallest eigenvalue of shape operator).
        T curvature_min;

        /// Principal direction corresponding to minimum curvature (eigenvector of shape operator).
        basevector<T, 3> principal_dir_min;

        /// Principal direction corresponding to maximum curvature (eigenvector of shape operator).
        basevector<T, 3> principal_dir_max;

        /// Gaussian curvature (product of principal curvatures: K = k1 * k2).
        T gauss_curvature;

        /// Mean curvature (average of principal curvatures: H = (k1 + k2) / 2).
        T mean_curvature;

        /// Surface type label for this vertex.
        vertex_label label;

        std::set<size_t> neighbour_vertices;

        /**
         * @brief Construct a vertex with id and coordinates.
         * @param _id Vertex identifier.
         * @param x X coordinate.
         * @param y Y coordinate.
         * @param z Z coordinate.
         */
        vertex(size_t _id, T x, T y, T z)
            : id(_id), coords(basevector<T, 3>(x, y, z)) {
        }
    };

    /**
     * @brief Triangle face stored in the mesh.
     *
     * The face stores three vertex ids (v1,v2,v3), a pointer to one of its
     * half-edges, a computed geometric normal, face center and triangle area.
     *
     * @tparam T numeric type for geometric quantities.
     */
    template <typename T>
    struct face {
        /// Unique face id.
        size_t id;

        /// Vertex ids (indices into the mesh vertex map).
        size_t v1, v2, v3;

        /// Pointer to one of the three half-edges that bound this face.
        half_edge<T>* first_edge;

        /// Unit-length face normal (computed by compute_face_normals()).
        basevector<T, 3> normal;

        /// Face centroid (computed by compute_face_properties()).
        basevector<T, 3> center;

        /// Face area (computed by compute_face_properties()).
        T area;

        /**
         * @brief Construct a triangular face with the provided id and vertices.
         * @param id_ Face identifier.
         * @param v1_ First vertex id.
         * @param v2_ Second vertex id.
         * @param v3_ Third vertex id.
         */
        face(size_t id_, size_t v1_, size_t v2_, size_t v3_)
            : id(id_), v1(v1_), v2(v2_), v3(v3_), first_edge(nullptr) {
        }

        std::pair<size_t, size_t> get_other_vertices(size_t v_id) {
            if (v_id == v1)
                return std::pair<size_t, size_t>(v2, v3);
            else if (v_id == v2)
                return std::pair<size_t, size_t>(v3, v1);
            else // v_id == v3
                return std::pair<size_t, size_t>(v1, v2);
        }
    };

    /**
     * @brief Half-edge element for representing directed edge between two vertices.
     *
     * Each triangular face has three half-edges linked cyclically by `next`.
     * Each half-edge may also have a `twin` pointer to the opposite half-edge
     * on an adjacent face. `m_face` points to the face that this half-edge borders.
     *
     * @tparam T numeric type (unused directly but kept templated for symmetry).
     */
    template <typename T>
    struct half_edge {
        /// Source vertex id for this directed half-edge.
        size_t source;

        /// Target vertex id for this directed half-edge.
        size_t target;

        /// Pointer to twin half-edge (opposite direction along same undirected edge).
        half_edge* twin;

        /// Pointer to the face that this half-edge borders.
        face<T>* m_face;

        /// Next half-edge in the face cycle (counter-clockwise around the face).
        half_edge<T>* next;

        /**
         * @brief Construct a half-edge that goes from s to t.
         * @param s source vertex id.
         * @param t target vertex id.
         */
        half_edge(size_t s, size_t t)
            : source(s), target(t), twin(nullptr), m_face(nullptr), next(nullptr) {
        }
    };

    /**
    * @brief Simple half-edge mesh class storing vertices, faces and half-edges.
    *
    * The class uses std::map keyed by ids and edge descriptors for lookups.
    * It provides utilities to add geometry, orient faces consistently, compute
    * normals/centers/areas and read simple binary mesh files.
    *
    * Note: memory ownership is raw pointers allocated with new in this header.
    * The class destructor currently does not free allocated memory; the user is
    * responsible for lifetime management or extending the destructor.
    *
    * @tparam T numeric type for coordinates and geometric computations.
    */
    template <typename T>
    class half_edge_mesh {
        /**
        * @brief Flip the ordering/orientation of the triangle face `f`.
        *
        * This function inverts each half-edge in the face by swapping their
        * source and target fields and then updates the face's vertex indices
        * to maintain consistency with the new edge directions.
        *
        * @param f Pointer to the face to flip.
        */
        void flip_face(face<T>* f) {
            half_edge<T>* he = f->first_edge;
            size_t te[3];
            for (int i = 0; i < 3; ++i) {
                // Store target vertices to temporary array
                te[i] = he->target;
                // Swap source and target vertices
                size_t temp = he->source;
                he->source = he->target;
                he->target = temp;
                he = he->next;
            }
            // Update face vertex indices
            f->v1 = te[2];
            f->v2 = te[1];
            f->v3 = te[0];
        }
    public:
        T total_edge_length = T(0);
        // average_edge_length: is the mean length of edges, useful for scale - dependent algorithms
        // (e.g., curvature estimation, smoothing step sizes, etc.).
        T average_edge_length = T(0);

        /// Map from directed edge descriptor (source,target) to half-edge pointer.
        std::map<edge_descriptor, half_edge<T>*> half_edges;

        /// Map of vertex id to vertex pointer.
        std::map<size_t, vertex<T>*> vertices;

        /// Map of face id to face pointer.
        std::map<size_t, face<T>*> faces;

        /// Default constructor.
        half_edge_mesh() {};
        /// Default destructor (does not free heap-allocated vertices/edges/faces).
        ~half_edge_mesh() {};

        const std::map<size_t, vertex<T>*>& get_vertices() const {
            return vertices;
        }

        /**
        * @brief Add a vertex to the mesh (no-op if id already exists).
        * @param id Vertex identifier.
        * @param x X coordinate.
        * @param y Y coordinate.
        * @param z Z coordinate.
        * @return The id of the added or existing vertex.
        */
        size_t add_vertex(size_t id, T x, T y, T z) {
            if (vertices.find(id) != vertices.end())
                return id; // vertex already exists
            vertices[id] = new vertex(id, x, y, z); // add vertex to map
            return id;
        }

        /**
        * @brief Add a triangular face and its half-edges to the mesh.
        *
        * This function:
        * - Optionally flips the ordering of the provided vertices to avoid
        *   creating a duplicated directed edge (simple heuristic).
        * - Allocates three half-edge objects, links them in a cycle and tries
        *   to connect twins to already-existing opposite half-edges.
        * - Records the face id in each incident vertex's `incident_faces`.
        *
        * @param id Face identifier.
        * @param i1 First vertex id.
        * @param i2 Second vertex id.
        * @param i3 Third vertex id.
        * @return The id of the added face.
        */
        size_t add_face(size_t id, size_t i1, size_t i2, size_t i3) {
            bool flip = false;
            // TODO: check if edge is triple bound
            if (half_edges.find(edge_descriptor(i1, i2)) != half_edges.end())
                flip = true;
            if (half_edges.find(edge_descriptor(i2, i3)) != half_edges.end())
                flip = true;
            if (half_edges.find(edge_descriptor(i3, i1)) != half_edges.end())
                flip = true;

            if (flip) {
                size_t t = i2;
                i2 = i3;
                i3 = t;
            }

            face<T>* f = new face<T>(id, i1, i2, i3);
            faces[id] = f;

            half_edge<T>* he1 = new half_edge<T>(i1, i2);// half-edge from v1 to v2, starting half-edge
            half_edge<T>* he2 = new half_edge<T>(i2, i3); // half-edge from v2 to v3
            half_edge<T>* he3 = new half_edge<T>(i3, i1); // half-edge from v3 to v1
            half_edges[edge_descriptor(i1, i2)] = he1;
            total_edge_length += (basevector<T, 3>(vertices[i2]->coords - vertices[i1]->coords)).length();
            half_edges[edge_descriptor(i2, i3)] = he2;
            total_edge_length += (basevector<T, 3>(vertices[i3]->coords - vertices[i2]->coords)).length();
            half_edges[edge_descriptor(i3, i1)] = he3;
            total_edge_length += (basevector<T, 3>(vertices[i1]->coords - vertices[i3]->coords)).length();

            he1->m_face = f;
            he2->m_face = f;
            he3->m_face = f;
            f->first_edge = he1;
            // Link half-edges in a cycle
            he1->next = he2;
            he2->next = he3;
            he3->next = he1;

            // Attempt to find and link twins (opposite directed half-edges).
            he1->twin = get_half_edge(i2, i1);
            if (he1->twin) he1->twin->twin = he1;
            he2->twin = get_half_edge(i3, i2);
            if (he2->twin) he2->twin->twin = he2;
            he3->twin = get_half_edge(i1, i3);
            if (he3->twin) he3->twin->twin = he3;

            // Update incident faces for each vertex
            vertices[i1]->incident_faces.push_back(id);
            vertices[i2]->incident_faces.push_back(id);
            vertices[i3]->incident_faces.push_back(id);

            vertices[i1]->neighbour_vertices.insert(i2);
            vertices[i1]->neighbour_vertices.insert(i3);
            vertices[i2]->neighbour_vertices.insert(i1);
            vertices[i2]->neighbour_vertices.insert(i3);
            vertices[i3]->neighbour_vertices.insert(i1);
            vertices[i3]->neighbour_vertices.insert(i2);

            return id;
        }

        /**
        * @brief Retrieve a half-edge given source and target vertex ids.
        * @param v1 Source vertex id.
        * @param v2 Target vertex id.
        * @return Pointer to the half-edge if present, or nullptr otherwise.
        */
        half_edge<T>* get_half_edge(size_t v1, size_t v2) {
            auto it = half_edges.find(edge_descriptor(v1, v2));
            if (it == half_edges.end())
                return nullptr;
            return it->second;
        }

        /**
        * @brief Finalize the mesh by orienting faces consistently.
        *
        * Currently this calls `orient_mesh()` and returns 1 on completion.
        * Can be extended to perform validation / cleanup steps.
        *
        * @return 1 on success (current implementation).
        */
        int finalize_mesh() {
            orient_mesh();
            return 1;
        }

        /**
        * @brief Orient all faces so they have consistent winding across the mesh.
        *
        * The algorithm repeatedly scans faces and their half-edges; when an
        * orientation mismatch with a twin half-edge is detected the adjacent
        * face is flipped. The process repeats until no flips are required.
        *
        * @return 1 on completion.
        */
        int orient_mesh() {
            bool done = false;
            while (!done) {
                done = true;
                for (auto f_pair : faces) {
                    face<T>* f = f_pair.second;
                    half_edge<T>* he = f->first_edge;
                    for (int i = 0; i < 3; ++i) {
                        if (he->twin) {
                            // Check if the orientation is consistent
                            if (he->source != he->twin->target ||
                                he->target != he->twin->source) {
                                // Orientation mismatch, flip the whole face
                                face<T>* fface = he->twin->m_face;
                                flip_face(fface);
                                done = false; // Need another pass
                            }
                        }
                        he = he->next;
                    }
                }
            }
            return 1;
        }

        /**
        * @brief Compute and store unit-length face normals for every face.
        *
        * For each triangular face, the cross product of two edges is computed
        * and normalized to form the face normal. The normal is stored in the
        * corresponding face's `normal` member.
        */
        void compute_face_normals() {
            for (auto f_pair : faces) {
                face<T>* f = f_pair.second;
                f->normal = calc_normal<T>(vertices[f->v1]->coords, vertices[f->v2]->coords, vertices[f->v3]->coords);
            }
        }

        /**
        * @brief Translate mesh so its vertex centroid is at the origin.
        *
        * Computes the average of all vertex coordinates and subtracts it from
        * each vertex position.
        */
        void recenter() {
            T sumx = T(0), sumy = T(0), sumz = T(0);
            for (auto& v : vertices) {
                basevector<T, 3>& vc = v.second->coords;
                sumx += vc.x();
                sumy += vc.y();
                sumz += vc.z();
            }
            sumx /= T(vertices.size());
            sumy /= T(vertices.size());
            sumz /= T(vertices.size());
            for (auto& v : vertices) {
                basevector<T, 3>& vc = v.second->coords;
                vc.x() -= sumx;
                vc.y() -= sumy;
                vc.z() -= sumz;
            }
        }

        /**
        * @brief Compute and cache geometric properties for each face.
        *
        * For every face in the mesh this function:
        * - Fetches the three vertex positions `v1`, `v2` and `v3` from the
        *   global `vertices` map using the face's vertex indices.
        * - Computes the face centroid as the arithmetic mean of the three
        *   vertex positions and stores it in `f->center`.
        * - Computes the triangle area using `triangle_area` on the three edge
        *   vectors `(v1 - v2)`, `(v2 - v3)` and `(v3 - v1)` and stores it in
        *   `f->area`.
        *
        * The resulting `center` and `area` values can be used for geometric
        * processing routines such as integration over the surface, mass
        * properties, or visualization of per-face quantities.
        */
        void compute_face_properties() {
            for (auto f_pair : faces) {
                face<T>* f = f_pair.second;
                basevector<T, 3> v1 = vertices[f->v1]->coords;
                basevector<T, 3> v2 = vertices[f->v2]->coords;
                basevector<T, 3> v3 = vertices[f->v3]->coords;
                basevector<T, 3> center = (v1 + v2 + v3) / 3.0f;
                f->center = center;
                f->area = triangle_area<T>(v1 - v2, v2 - v3, v3 - v1);
            }
        }

        /**
        * @brief Compute per-vertex normals by area-weighted accumulation of face normals.
        *
        * The algorithm:
        *  - Initializes vertex normals to zero.
        *  - Adds each face's normal to its three incident vertices.
        *  - Normalizes each vertex normal to unit-length.
        *
        * This produces smooth shading normals suitable for rendering.
        */
        void compute_vertex_normals() {
            // Initialize all vertex normals to zero
            for (auto v_pair : vertices) {
                vertex<T>* v = v_pair.second;
                v->normal = basevector<T, 3>(0.0f, 0.0f, 0.0f);
            }
            // Accumulate face normals to vertex normals
            for (auto f_pair : faces) {
                face<T>* f = f_pair.second;
                basevector<T, 3> face_normal = f->normal;
                vertices[f->v1]->normal += face_normal;
                vertices[f->v2]->normal += face_normal;
                vertices[f->v3]->normal += face_normal;
            }
            // Normalize vertex normals
            for (auto v_pair : vertices) {
                vertex<T>* v = v_pair.second;
                v->normal.normalize();
            }
        }

        T compute_vertex_voronoi_area(size_t v_id) {
            T area_mixed = 0.0;
            basevector<T, 3>& P = vertices[v_id]->coords;

            vertex<T>& v = *(vertices[v_id]);

            // Iterate over all triangles sharing vertex v_id
            for (size_t f_id : v.incident_faces) {
                face<T>& f = *(faces[f_id]);
                std::pair<size_t, size_t> other = f.get_other_vertices(v_id);
                // Get vertices of the triangle (P, Q, R)
                basevector<T, 3>& Q = vertices[other.first]->coords;
                basevector<T, 3>& R = vertices[other.second]->coords;

                // Edge vectors
                basevector<T, 3> PR = R - P;
                basevector<T, 3> PQ = Q - P;
                basevector<T, 3> QR = R - Q;

                // Total area of this triangle
                T area_T = T(0.5) * PQ.cross(PR).length();

                // Check if triangle is obtuse
                // Angle at P: PQ . PR
                // Angle at Q: (-PQ) . QR
                // Angle at R: (-PR) . (-QR)
                bool is_obtuse = PQ.dot(PR) < 0 || (-PQ).dot(QR) < 0 || (-PR).dot(-QR) < 0;

                if (!is_obtuse) {
                    // 1. Voronoi Area using cotangent weights
                    // Area = 1/8 * (|PQ|^2 * cot(angle_at_R) + |PR|^2 * cot(angle_at_Q))
                    T cot_Q = PQ.dot(-QR) / PQ.cross(-QR).length();
                    T cot_R = PR.dot(QR) / PR.cross(QR).length();
                    area_mixed += (PQ.dot(PQ) * cot_R + PR.dot(PR) * cot_Q) / T(8.0);
                }
                else {
                    // 2. Obtuse case
                    if (PQ.dot(PR) < 0) {
                        // Angle at P is obtuse
                        area_mixed += area_T / T(2.0);
                    }
                    else {
                        // Angle at Q or R is obtuse
                        area_mixed += area_T / T(4.0);
                    }
                }
            }
            return area_mixed;
        }

        void compute_vertex_voronoi_areas() {
            // Initialize all vertex normals to zero
            for (auto v_pair : vertices) {
                vertex<T>* v = v_pair.second;
                T v_area = compute_vertex_voronoi_area(v->id);
                v->voronoi_area = v_area;
            }
        }

        /**
        * @brief Get the ordered 1-ring neighbor vertices of a given vertex.
        *
        * For the vertex identified by `v_id`, this function collects all incident
        * faces, extracts from each face the opposite edge (the two vertices that
        * are not `v_id`), and stores these edges in `ring_edges`. The helper
        * function `order_chain` is then used to order these edges into a
        * continuous chain around the central vertex. Finally, the first vertex
        * of each ordered edge is copied into the result vector.
        *
        * The returned vector thus contains the neighboring vertex ids of `v_id`
        * in a (locally) consistent circular order, which is useful for
        * differential geometry operators, curvature estimation, or any operation
        * that requires traversing the 1-ring in order.
        *
        * @param v_id The id of the central vertex whose neighbors are requested.
        * @return A vector of neighboring vertex ids ordered around `v_id`.
        */
        std::vector<size_t> get_vertex_neighbours(size_t v_id) {
            vertex<T>& v = *(vertices[v_id]);
            std::vector<std::pair<size_t, size_t>> ring_edges;
            ring_edges.reserve(v.incident_faces.size());
            for (size_t face_id : v.incident_faces) {
                face<T>* f = faces[face_id];
                auto ov = f->get_other_vertices(v_id);
                ring_edges.emplace_back(ov.first, ov.second);
            }
            order_chain(ring_edges);
            std::vector<size_t> result;
            result.reserve(ring_edges.size());
            for (const auto& e : ring_edges)
                result.push_back(e.first);
            return result;
        }

        basematrix<T, 2, 3> get_bounding_box() {
            if (vertices.size() == 0) {
                return basematrix<T,2,3>(); // empty matrix
            }
            T min_x = std::numeric_limits<T>::max();
            T min_y = std::numeric_limits<T>::max();
            T min_z = std::numeric_limits<T>::max();
            T max_x = std::numeric_limits<T>::lowest();
            T max_y = std::numeric_limits<T>::lowest();
            T max_z = std::numeric_limits<T>::lowest();

            for (const auto& v_pair : vertices) {
                const basevector<T, 3>& coords = v_pair.second->coords;
                if (coords.x() < min_x) min_x = coords.x();
                if (coords.y() < min_y) min_y = coords.y();
                if (coords.z() < min_z) min_z = coords.z();
                if (coords.x() > max_x) max_x = coords.x();
                if (coords.y() > max_y) max_y = coords.y();
                if (coords.z() > max_z) max_z = coords.z();
            }
            basematrix<T, 2, 3> bbox({min_x, min_y, min_z, max_x, max_y, max_z});
            return bbox;
        }

        void translate(const basevector<T, 3>& offset) {
            for (auto& v_pair : vertices) {
                basevector<T, 3>& coords = v_pair.second->coords;
                coords = coords + offset;
            }
        }

        /**
        * @brief Debug print of the mesh details to stdout (enabled only in _DEBUG).
        *
        * Prints per-face vertex indices and half-edge information, then lists
        * each vertex's incident faces. Useful for small meshes during
        * development.
        */
        void print_mesh() {
#ifdef _DEBUG
            std::cout << "Mesh built with " << vertices.size() << " vertices and " << faces.size() << " faces." << std::endl;
            for (auto f : faces) {
                std::cout << "Face " << f.first << ": vertices (" << f.second->v1 << ", " << f.second->v2 << ", " << f.second->v3 << ")" << std::endl;
#if 0
                std::cout << "\t" << f.second->v1 << "(" << vertices[f.second->v1]->coords[0] << "," << vertices[f.second->v1]->coords[1] << "," << vertices[f.second->v1]->coords[2] << "), ";
                std::cout << f.second->v2 << "(" << vertices[f.second->v2]->coords[0] << "," << vertices[f.second->v2]->coords[1] << "," << vertices[f.second->v2]->coords[2] << "), ";
                std::cout << f.second->v3 << "(" << vertices[f.second->v3]->coords[0] << "," << vertices[f.second->v3]->coords[1] << "," << vertices[f.second->v3]->coords[2] << ")\n";
#endif
                std::cout << "  Half-edges: \n";
                half_edge* he = f.second->first_edge;
                for (int i = 0; i < 3; ++i) {
                    std::cout << "\t(" << he->source << " -> " << he->target << ") " << "twin:" << he->twin << "\n";
                    he = he->next;
                }
            }
            for (auto v : vertices) {
                std::cout << "Vertex " << v.first << ": incident faces ";
                for (auto fid : v.second->incident_faces) {
                    std::cout << fid << " ";
                }
                std::cout << std::endl;
            }
#endif
        }

        /**
   * @brief Minimal debug print (counts only).
   *
   * Prints the number of vertices and faces when in _DEBUG mode.
   */
        void print_mesh_quick() {
#ifdef _DEBUG
            std::cout << "Mesh built with " << vertices.size() << " vertices and " << faces.size() << " faces." << std::endl;
#endif
        }
    };

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
        std::vector<std::tuple<size_t, T, T, T>> vertices;
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
    template <typename T>
    struct mesh_def<T>* read_mesh_file(const std::string& fnm) {
        std::ifstream rf(fnm, std::ios::in | std::ios::binary);

        if (!rf)
            return nullptr;
        mesh_def<T>* mesh = new mesh_def<T>;
        size_t d, i;
        rf.read((char*)&d, sizeof(size_t));
        for (i = 0; i < d; ++i) {
            size_t id;
            float x, y, z; // the file is built with 'float'!!!
            rf.read((char*)&id, sizeof(size_t));
            rf.read((char*)&x, sizeof(float));
            rf.read((char*)&y, sizeof(float));
            rf.read((char*)&z, sizeof(float));
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

    /**
  * @brief Construct a `half_edge_mesh` from an in-memory `mesh_def`.
  *
  * This function:
  *  - Allocates an empty half_edge_mesh.
  *  - Adds all vertices via `add_vertex`.
  *  - Adds all faces via `add_face`.
  *  - Computes face normals, centers and vertex normals.
  *
  * @param _mesh_def Reference to the mesh definition to parse.
  * @return Pointer to a newly allocated `half_edge_mesh<T>`.
  */
    template <typename T>
    half_edge_mesh<T>* parse_mesh(const mesh_def<T>& _mesh_def) {
        half_edge_mesh<T>* mesh = new half_edge_mesh<T>;
        // add vertices
        for (const std::tuple<size_t, T, T, T>& vert : _mesh_def.vertices) {
            mesh->add_vertex(std::get<0>(vert), T(std::get<1>(vert)), T(std::get<2>(vert)), T(std::get<3>(vert)));
        }

        // add faces
        for (const std::tuple<size_t, size_t, size_t, size_t>& fc : _mesh_def.faces) {
            mesh->add_face(std::get<0>(fc), std::get<1>(fc), std::get<2>(fc), std::get<3>(fc));
        }
        mesh->average_edge_length = mesh->total_edge_length / mesh->half_edges.size();
        // mesh->orient_mesh();
        mesh->recenter();
        mesh->compute_face_normals();
        mesh->compute_face_properties();
        mesh->compute_vertex_normals();

        return mesh;
    }

    /**
  * @brief Convenience: read a mesh file and parse it into a half-edge mesh.
  * @param fnm Path to the binary mesh file.
  * @return Pointer to newly created `half_edge_mesh<T>` or nullptr on failure.
  */
    template <typename T>
    half_edge_mesh<T>* create_from_mesh_file(const std::string& fnm) {
        if (file_extension(fnm) != "prim")
            return nullptr;
        struct mesh_def<T>* mesh = read_mesh_file<T>(fnm);
        if (mesh == nullptr)
            return nullptr;
        half_edge_mesh<T>* he_mesh = parse_mesh<T>(*mesh);
        delete mesh;
        return he_mesh;
    }

    // The following factory functions provide shortcuts to create common unit meshes
    // from the resources folder. They all call create_from_mesh_file with a fixed path.

    /**
  * @brief Create a unit cone mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_cone() {
        return create_from_mesh_file<T>("resources\\models\\unit_cone.prim");
    }

    /**
  * @brief Create a unit cube mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_cube() {
        return create_from_mesh_file<T>("resources\\models\\unit_cube.prim");
    }

    /**
  * @brief Create a unit cylinder mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_cylinder() {
        return create_from_mesh_file<T>("resources\\models\\unit_cylinder.prim");
    }

    /**
  * @brief Create a unit dodecahedron mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_dodecahedron() {
        return create_from_mesh_file<T>("resources\\models\\unit_dodecahedron.prim");
    }

    /**
  * @brief Create a unit icosahedron mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_icosahedron() {
        return create_from_mesh_file<T>("resources\\models\\unit_icosahedron.prim");
    }

    /**
  * @brief Create a unit octahedron mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_octa() {
        return create_from_mesh_file<T>("resources\\models\\unit_octa.prim");
    }
    /**
  * @brief Create a unit penta mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_penta() {
        return create_from_mesh_file<T>("resources\\models\\unit_penta.prim");
    }

    /**
  * @brief Create a unit plane mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_plane() {
        return create_from_mesh_file<T>("resources\\models\\unit_plane.prim");
    }

    /**
  * @brief Create a unit-diameter sphere mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_sphere() {
        return create_from_mesh_file<T>("resources\\models\\unit_sphere.prim");
    }

    /**
  * @brief Create a unit tetrahedron mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_tetra() {
        return create_from_mesh_file<T>("resources\\models\\unit_tetra.prim");
    }
    /**
  * @brief Create a unit torus mesh from resources.
  * @return Newly allocated half_edge_mesh or nullptr if file missing.
  */
    template <typename T>
    half_edge_mesh<T>* create_unit_torus() {
        return create_from_mesh_file<T>("resources\\models\\unit_torus.prim");
    }
}
#endif // __mesh_h__
