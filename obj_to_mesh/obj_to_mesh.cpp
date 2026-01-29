// Converts Wavefront OBJ mesh files into a compact binary `.prim` format.
//
// This tool parses a subset of the OBJ specification, extracting only
// geometric vertex positions (`v` lines) and triangular faces (`f` lines),
// and then writes them into a binary file alongside simple integer IDs.
//
// The binary layout of the generated `.prim` file is:
//
//   1. `size_t` vertex_count
//      Repeated `vertex_count` times:
//        a. `size_t` vertex_id
//        b. `float`  x
//        c. `float`  y
//        d. `float`  z
//
//   2. `size_t` face_count
//      Repeated `face_count` times:
//        a. `size_t` face_id
//        b. `size_t` v1   // index of first vertex (1-based, as in OBJ)
//        c. `size_t` v2   // index of second vertex (1-based, as in OBJ)
//        d. `size_t` v3   // index of third vertex (1-based, as in OBJ)
//
// The mapping between OBJ constructs and internal structures is:
//
//   - Each `v` record becomes an entry in `vertex`:
//       * `id` is a sequential 1-based identifier assigned in load order.
//       * `x`, `y`, `z` are parsed as 32-bit floating point coordinates.
//   - Each `f` record becomes an entry in `face`:
//       * `id` is a sequential 1-based identifier assigned in load order.
//       * `v1`, `v2`, `v3` store vertex indices parsed from the face
//         definition (only the positional index before any `/` is used).
//
// Limitations and assumptions:
//
//   - Only lines beginning with `v` (vertex positions) and `f` (faces) are
//     interpreted; all other OBJ statements are ignored.
//   - Faces are assumed to be triangles and must contain exactly three
//     vertex references. Polygons with more than three vertices are not
//     supported and will be partially or incorrectly interpreted if present.
//   - For face vertex references, only the first component before a `/`
//     is used (e.g. `f 1/2/3 4/5/6 7/8/9` uses only `1`, `4`, `7`).
//   - Vertex and face indices are treated as positive 1-based indices.
//   - Input is read line by line; malformed lines generate warnings or
//     errors to `std::cerr` and are skipped when possible without
//     terminating the entire conversion.
//
// Error handling and robustness:
//
//   - If the input OBJ file cannot be opened, `parse_obj_model` logs an
//     error to `std::cerr` and returns `false`.
//   - When parsing vertices (`v`):
//       * Lines with fewer than three coordinate values emit a warning and
//         are skipped.
//       * Conversion errors on coordinates are caught via `std::exception`;
//         the offending line is reported and skipped.
//   - When parsing faces (`f`):
//       * Lines with fewer than three vertex references emit a warning and
//         are skipped.
//       * Face vertex tokens that cannot be split or converted into a valid
//         index emit a warning and are skipped.
//   - If the output `.prim` file cannot be created, an error is logged and
//     `false` is returned.
//   - After writing, the stream state is checked; if writing failed, an
//     error is logged and `false` is returned.
//
// Utility helpers:
//
//   - `splitString(const std::string& str, char delimiter = ' ')`
//       Splits the input string into non-empty tokens using the specified
//       delimiter. This is used for:
//         * Splitting input lines by spaces into OBJ tokens.
//         * Splitting face components (e.g. `1/2/3`) by `'/'` to extract
//           the vertex index.
//       The function pre-reserves a small token capacity to reduce
//       allocations for the common use cases.
//
//   - `change_file_extension(const std::string& fname, const std::string& new_ext)`
//       Returns a copy of `fname` with its extension replaced by `new_ext`.
//       If `fname` has no extension, `.` and `new_ext` are appended.
//       This is used to derive the `.prim` output filename from the
//       `.obj` input filename.
//
// Data structures:
//
//   - `struct vertex`
//       Represents a single geometric vertex:
//         * `id`   : Sequential 1-based identifier assigned at load time.
//         * `x,y,z`: Position components in object space.
//
//   - `struct face`
//       Represents a triangular face by vertex indices:
//         * `id`   : Sequential 1-based identifier assigned at load time.
//         * `v1,v2,v3`: Indices of the three vertices forming the triangle.
//
// Main conversion routine:
//
//   - `bool parse_obj_model(const std::string& fnm)`
//       Steps:
//         1. Open the OBJ file identified by `fnm`.
//         2. Initialize vertex and face containers with reserved capacity
//            to reduce reallocations for typical model sizes.
//         3. Read the file line by line, tracking the current line number.
//         4. For each non-empty line, tokenize and interpret only `v` and
//            `f` records, constructing `vertex` and `face` entries.
//         5. After the entire file is parsed, open the corresponding `.prim`
//            output file and serialize the vertex and face data in the
//            binary layout described above.
//         6. On success, log a summary (`<file> parsed (N vertices, M faces)`)
//            and return `true`.
//
// Program entry point:
//
//   - `int main(int argc, char** argv)`
//       Currently hardcodes a set of canonical unit geometry OBJ files
//       located under `resources\models\` (e.g. `unit_cone.obj`,
//       `unit_cube.obj`, etc.). For each file, `parse_obj_model` is
//       invoked to generate its corresponding `.prim` representation.
//       The return code is `0` regardless of per-file parse success;
//       errors are reported individually to standard error.

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cstdlib>

struct vertex {
    size_t id;
    float x;
    float y;
    float z;
    vertex(size_t _id, float _x, float _y, float _z) :id(_id), x(_x), y(_y), z(_z) {}
};

struct face {
    size_t id;
    size_t v1, v2, v3;
    face(size_t _id, size_t _v1, size_t _v2, size_t _v3) :id(_id), v1(_v1), v2(_v2), v3(_v3) {}
};

inline std::vector<std::string> splitString(const std::string& str, char delimiter = ' ') {
    std::vector<std::string> tokens;
    tokens.reserve(4);
    std::stringstream ss(str);
    std::string token;
    while (getline(ss, token, delimiter)) {
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }
    return tokens;
}

inline std::string change_file_extension(const std::string& fname, const std::string& new_ext) {
    std::size_t dot = fname.find_last_of(".");
    if (dot == std::string::npos)
        return fname + "." + new_ext;
    return fname.substr(0, dot + 1) + new_ext;
}

bool parse_obj_model(const std::string& fnm) {
    std::ifstream mdl(fnm);
    if (!mdl.is_open())
    {
        std::cerr << "Error: Could not open file " << fnm << "\n";
        return false;
    }

    std::vector<vertex> vertices;
    std::vector<face> faces;
    vertices.reserve(1024);
    faces.reserve(2048);

    size_t vertex_count = 1;
    size_t face_count = 1;
    std::string line;
    size_t line_num = 0;

    while (std::getline(mdl, line))
    {
        ++line_num;
        if (line.empty()) continue;

        auto tokens = splitString(line);
        if (tokens.empty()) continue;

        try {
            if (tokens[0] == "v") {
                if (tokens.size() < 4) {
                    std::cerr << "Warning: Invalid vertex at line " << line_num << "\n";
                    continue;
                }
                char* end;
                float x = std::strtof(tokens[1].c_str(), &end);
                float y = std::strtof(tokens[2].c_str(), &end);
                float z = std::strtof(tokens[3].c_str(), &end);
                vertices.emplace_back(vertex_count, x, y, z);
                ++vertex_count;
            }
            else if (tokens[0] == "f") {
                if (tokens.size() < 4) {
                    std::cerr << "Warning: Invalid face at line " << line_num << "\n";
                    continue;
                }
                auto v1 = splitString(tokens[1], '/');
                auto v2 = splitString(tokens[2], '/');
                auto v3 = splitString(tokens[3], '/');
                
                if (v1.empty() || v2.empty() || v3.empty()) {
                    std::cerr << "Warning: Invalid face indices at line " << line_num << "\n";
                    continue;
                }
                
                size_t idx1 = std::stoull(v1[0]);
                size_t idx2 = std::stoull(v2[0]);
                size_t idx3 = std::stoull(v3[0]);
                faces.emplace_back(face_count, idx1, idx2, idx3);
                ++face_count;
            }
        }
        catch (const std::exception& e) {
            std::cerr << "Error parsing line " << line_num << ": " << e.what() << "\n";
            continue;
        }
    }
    mdl.close();

    std::ofstream odl(change_file_extension(fnm, "prim"), std::ios::out | std::ios::binary);
    if (!odl.is_open())
    {
        std::cerr << "Error: Could not create output file for " << fnm << "\n";
        return false;
    }

    size_t s = vertices.size();
    odl.write((const char*)&s, sizeof(size_t));
    for (const auto& vtx : vertices) {
        odl.write((const char*)&vtx.id, sizeof(size_t));
        odl.write((const char*)&vtx.x, sizeof(float));
        odl.write((const char*)&vtx.y, sizeof(float));
        odl.write((const char*)&vtx.z, sizeof(float));
    }
    
    s = faces.size();
    odl.write((const char*)&s, sizeof(size_t));
    for (const auto& f : faces) {
        odl.write((const char*)&f.id, sizeof(size_t));
        odl.write((const char*)&f.v1, sizeof(size_t));
        odl.write((const char*)&f.v2, sizeof(size_t));
        odl.write((const char*)&f.v3, sizeof(size_t));
    }

    if (!odl.good()) {
        std::cerr << "Error: Failed to write output file for " << fnm << "\n";
        return false;
    }

    std::cout << fnm << " parsed (" << vertices.size() << " vertices, " << faces.size() << " faces)\n";

    return true;
}

int main(int argc, char** argv) {
    parse_obj_model("resources\\models\\unit_cone.obj");
    parse_obj_model("resources\\models\\unit_cube.obj");
    parse_obj_model("resources\\models\\unit_cylinder.obj");
    parse_obj_model("resources\\models\\unit_dodecahedron.obj");
    parse_obj_model("resources\\models\\unit_icosahedron.obj");
    parse_obj_model("resources\\models\\unit_octa.obj");
    parse_obj_model("resources\\models\\unit_penta.obj");
    parse_obj_model("resources\\models\\unit_plane.obj");
    parse_obj_model("resources\\models\\unit_sphere.obj");
    parse_obj_model("resources\\models\\unit_tetra.obj");
    parse_obj_model("resources\\models\\unit_torus.obj");
    return 0;
}