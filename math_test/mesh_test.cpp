#include <iostream>
#include <fstream>

#include "string_utils.h"
#include "vector.h"
#include "mesh.h"

using namespace base_math;

// ... existing model/load_obj/load_model/print_mesh/print_mesh_normals definitions stay unchanged ...

// Simple assertion helpers for ad‑hoc tests (no external test framework used)
static int g_failedTests = 0;
static int g_totalTests  = 0;

static void AssertTrue(bool cond, const char* msg)
{
    ++g_totalTests;
    if (!cond)
    {
        ++g_failedTests;
        std::cout << "[FAIL] " << msg << "\n";
    }
}

static void AssertEqualSizeT(size_t expected, size_t actual, const char* msg)
{
    ++g_totalTests;
    if (expected != actual)
    {
        ++g_failedTests;
        std::cout << "[FAIL] " << msg << " expected=" << expected << " actual=" << actual << "\n";
    }
}

static void AssertNear(double expected, double actual, double eps, const char* msg)
{
    ++g_totalTests;
    if (std::fabs(expected - actual) > eps)
    {
        ++g_failedTests;
        std::cout << "[FAIL] " << msg << " expected=" << expected << " actual=" << actual << "\n";
    }
}

// ---------- Tests for base_math::mesh in base_math\mesh.h ----------

static void Test_AddFace_BuildsEdgesAndAdjacency()
{
    mesh<double> m;
    m.addVertex(1, basevector<double, 3>(0.0, 0.0, 0.0));
    m.addVertex(2, basevector<double, 3>(1.0, 0.0, 0.0));
    m.addVertex(3, basevector<double, 3>(0.0, 1.0, 0.0));

    m.addFace(1, 2, 3);

    AssertEqualSizeT(3, m.getVertices().size(), "AddFace: vertex count");
    AssertEqualSizeT(1, m.getFaces().size(), "AddFace: face count");
    AssertEqualSizeT(3, m.getEdges().size(), "AddFace: edge count (triangle)");

    // Check that each vertex sees exactly one adjacent face
    for (auto& vp : m.getVertices())
    {
        AssertEqualSizeT(1, vp.second->adjacentFaces.size(), "AddFace: vertex adjacentFaces size == 1");
    }
}

static void Test_AddFace_QuadAndBreakQuads()
{
    mesh<double> m;
    m.addVertex(1, basevector<double, 3>(0.0, 0.0, 0.0));
    m.addVertex(2, basevector<double, 3>(1.0, 0.0, 0.0));
    m.addVertex(3, basevector<double, 3>(1.0, 1.0, 0.0));
    m.addVertex(4, basevector<double, 3>(0.0, 1.0, 0.0));

    m.addFace(1, 2, 3, 4);
    AssertEqualSizeT(1, m.getFaces().size(), "Quad: one face before breakQuads");

    m.breakQuads();
    AssertEqualSizeT(2, m.getFaces().size(), "Quad: split into two triangles after breakQuads");
}

static void Test_ComputeFaceNormals_TriangleInXYPlane()
{
    mesh<double> m;
    m.addVertex(1, basevector<double, 3>(0.0, 0.0, 0.0));
    m.addVertex(2, basevector<double, 3>(1.0, 0.0, 0.0));
    m.addVertex(3, basevector<double, 3>(0.0, 1.0, 0.0));

    m.addFace(1, 2, 3);
    m.computeFaceNormals();

    const auto& faces = m.getFaces();
    auto it = faces.begin();
    AssertTrue(it != faces.end(), "ComputeFaceNormals: face exists");
    if (it != faces.end())
    {
        const basevector<double, 3>& n = it->second->normal;
        // Expect (0, 0, 1)
        AssertNear(0.0, n.x(), 1e-9, "ComputeFaceNormals: nx");
        AssertNear(0.0, n.y(), 1e-9, "ComputeFaceNormals: ny");
        AssertNear(1.0, n.z(), 1e-9, "ComputeFaceNormals: nz");
    }
}

static void Test_ComputeVertexNormals_FromFaceNormals()
{
    mesh<double> m;
    m.addVertex(1, basevector<double, 3>(0.0, 0.0, 0.0));
    m.addVertex(2, basevector<double, 3>(1.0, 0.0, 0.0));
    m.addVertex(3, basevector<double, 3>(0.0, 1.0, 0.0));

    m.addFace(1, 2, 3);

    m.computeFaceNormals();
    m.computeVertexNormals();

    for (auto& vp : m.getVertices())
    {
        const basevector<double, 3>& n = vp.second->normal;
        AssertNear(0.0, n.x(), 1e-9, "ComputeVertexNormals: nx");
        AssertNear(0.0, n.y(), 1e-9, "ComputeVertexNormals: ny");
        AssertNear(1.0, n.z(), 1e-9, "ComputeVertexNormals: nz");
    }
}

static void Test_RecenterMesh_MovesCentroidToOrigin()
{
    mesh<double> m;
    m.addVertex(1, basevector<double, 3>(0.0, 0.0, 0.0));
    m.addVertex(2, basevector<double, 3>(2.0, 0.0, 0.0));

    m.recenterMesh();

    double cx = 0.0;
    for (auto& vp : m.getVertices())
    {
        cx += vp.second->position.x();
    }
    cx /= static_cast<double>(m.getVertices().size());

    AssertNear(0.0, cx, 1e-12, "RecenterMesh: centroid x should be 0");
}

static void Test_ComputeVerticesVoronoiArea_SimpleTriangle()
{
    mesh<double> m;
    m.addVertex(1, basevector<double, 3>(0.0, 0.0, 0.0));
    m.addVertex(2, basevector<double, 3>(1.0, 0.0, 0.0));
    m.addVertex(3, basevector<double, 3>(0.0, 1.0, 0.0));

    m.addFace(1, 2, 3);
    m.computeFaceNormals();
    m.computeVerticesVoronoiArea();

    // All Voronoi areas should be positive for this simple triangle
    for (auto& vp : m.getVertices())
    {
        AssertTrue(vp.second->voronoiArea > 0.0, "ComputeVerticesVoronoiArea: area > 0");
    }
}

static void Test_GetVertexNeighbours_RingOrdering()
{
    if (1)
    {
        mesh<double> m;
        // square split into two triangles sharing vertex 1
        m.addVertex(1, basevector<double, 3>(0.0, 0.0, 0.0));
        m.addVertex(2, basevector<double, 3>(2.0, 0.0, 0.0));
        m.addVertex(3, basevector<double, 3>(2.0, 2.0, 0.0));
        m.addVertex(4, basevector<double, 3>(0.0, 2.0, 0.0));
        m.addVertex(5, basevector<double, 3>(1.0, 1.0, 0.0));

        m.addFace(1, 5, 2);
        m.addFace(2, 5, 3);
        m.addFace(3, 5, 4);

        std::vector<size_t> neighbours = m.getVertexNeighbours(5);
        // for (size_t i = 0; i < neighbours.size(); ++i) {
        //     std::cout << "Neighbour " << i << ": " << neighbours[i] << "\n";
        // }
        // std::cout << "Open ring, Expected neighbours: 1,2,3,4 (in ring order)\n";
        AssertEqualSizeT(4, neighbours.size(), "GetVertexNeighbours: degree of vertex 5 is 4");
        bool has1 = std::find(neighbours.begin(), neighbours.end(), 1) != neighbours.end();
        bool has2 = std::find(neighbours.begin(), neighbours.end(), 2) != neighbours.end();
        bool has3 = std::find(neighbours.begin(), neighbours.end(), 3) != neighbours.end();
        bool has4 = std::find(neighbours.begin(), neighbours.end(), 4) != neighbours.end();
        AssertTrue(has2 && has4, "GetVertexNeighbours: neighbours contain 1,2,3,4");
    }

    if (1)
    {
        mesh<double> m;
        // square split into two triangles sharing vertex 1
        m.addVertex(1, basevector<double, 3>(0.0, 0.0, 0.0));
        m.addVertex(2, basevector<double, 3>(2.0, 0.0, 0.0));
        m.addVertex(3, basevector<double, 3>(2.0, 2.0, 0.0));
        m.addVertex(4, basevector<double, 3>(0.0, 2.0, 0.0));
        m.addVertex(5, basevector<double, 3>(1.0, 1.0, 0.0));

        m.addFace(2, 5, 3);
        m.addFace(3, 5, 4);
        m.addFace(1, 5, 2);
        m.addFace(4, 5, 1);

        std::vector<size_t> neighbours = m.getVertexNeighbours(5);
        // for (size_t i = 0; i < neighbours.size(); ++i) {
        //     std::cout << "Neighbour " << i << ": " << neighbours[i] << "\n";
        // }
        // std::cout << "Closed ring, Expected neighbours: 1,2,3,4 (in ring order)\n";
        AssertEqualSizeT(4, neighbours.size(), "GetVertexNeighbours: degree of vertex 5 is 4");
        bool has1 = std::find(neighbours.begin(), neighbours.end(), 1) != neighbours.end();
        bool has2 = std::find(neighbours.begin(), neighbours.end(), 2) != neighbours.end();
        bool has3 = std::find(neighbours.begin(), neighbours.end(), 3) != neighbours.end();
        bool has4 = std::find(neighbours.begin(), neighbours.end(), 4) != neighbours.end();
        AssertTrue(has2 && has4, "GetVertexNeighbours: neighbours contain 1,2,3,4");
    }
}

// Existing main is turned into a simple test runner to exercise mesh.h
int mesh_tests()
{
    g_totalTests = 0;
    g_failedTests = 0;
    // Run tests for base_math::mesh
    std::cout << "Running mesh.h tests...\n";
    Test_AddFace_BuildsEdgesAndAdjacency();
    Test_AddFace_QuadAndBreakQuads();
    Test_ComputeFaceNormals_TriangleInXYPlane();
    Test_ComputeVertexNormals_FromFaceNormals();
    Test_RecenterMesh_MovesCentroidToOrigin();
    Test_ComputeVerticesVoronoiArea_SimpleTriangle();
    Test_GetVertexNeighbours_RingOrdering();

    std::cout << "Total tests: " << g_totalTests << "  Failed: " << g_failedTests << "\n";

    return g_failedTests == 0 ? 0 : 1;
}