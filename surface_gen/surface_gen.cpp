// a surface generator for a 3D graphing application. 
// The program will take user input for the dimensions of the surface and generate a 3D mesh accordingly.
#include <iostream>
#include <string>
#include <cmath>

#include "vector.h"
#include "mesh.h"

using namespace base_math;

static void obj_file_dump(const std::string& fnm, const mesh<double>* mdl) {
    std::ofstream out(fnm);
    if (!out.is_open()) {
        return;
    }
    out << "o Part" << std::endl;
    for (const auto& v_pair : mdl->vertices) {
        const basevector<double, 3>& pos = v_pair.second->position;
        out << "v " << pos.x() << " " << pos.y() << " " << pos.z() << std::endl;
    }
    for (const auto& f_pair : mdl->faces) {
        const meshFace<double>* face = f_pair.second;
        out << "f";
        for (const auto& vert : face->vertices) {
            out << " " << vert->id;
        }
        out << std::endl;
    }
    out.close();
}


std::vector<dvec2> sinus_curve(double xstart, double xend, double step) {
    std::vector<dvec2> points;
    for (double x = xstart; x < xend; x += step) {
        double y = sin(x)*cos(x);
        points.push_back(dvec2(x, y));
    }
    return points;
}

std::vector<dvec2> hyperbolic_curve(double xstart, double xend, double step) {
    std::vector<dvec2> points;
    for (double x = xstart; x < xend; x += step) {
        double y = x * x;
        points.push_back(dvec2(x, y));
    }
    return points;
}

mesh<double>* generate_surface(double x, double y, double z) {
    mesh<double>* surface_mesh = new mesh<double>;
    // Generate vertices for the surface based on the sinus curve
    std::vector<dvec2> curve_points = hyperbolic_curve(-3, 3, 0.1);

    int count = 1; // OBJ format is 1-based
    for (double z_val = -3; z_val < 3; z_val += 0.1) {
        for (const dvec2& point : curve_points) {
            double offset = z_val * z_val;
            surface_mesh->addVertex(count, dvec3(point.x(), point.y()+offset, z_val)); // Extrude along z-axis
            ++count;
        }
    }

    // Create faces for the surface by extruding the curve along the y-axis
    size_t num_points = curve_points.size();
    size_t total_vertices = surface_mesh->vertices.size();
    size_t sripes = total_vertices / num_points; // Number of extruded layers
    for (size_t i = 0; i < sripes - 1; ++i) {
        size_t base = i * num_points;
        for (size_t j = 0; j < num_points-1; ++j) {
            size_t id0 = 1 + base + j;
            size_t id1 = 1 + base + j + 1;
            size_t id2 = 1 + base + num_points + j + 1;
            size_t id3 = 1 + base + num_points  + j ;
            surface_mesh->addFace(id0, id1, id2);
            surface_mesh->addFace(id0, id2, id3);
        }
    }
    return surface_mesh;
}

base_math::mesh<double>* create_sphere_mesh(double radius) {
    base_math::mesh<double>* ms = new base_math::mesh<double>;

    int sectorCount = 48;
    int stackCount = 48;

    double x, y, z, xz;                              // vertex position
    double lengthInv = 1.0f / radius;    // normal
    double sectorStep = 2 * PI<double> / sectorCount;
    double stackStep = PI<double> / stackCount;
    double sectorAngle, stackAngle;
    int vertexCount = 0;
    for (int i = 0; i <= stackCount; ++i)
    {
        stackAngle = -(PI<double> / 2 - i * stackStep);        // starting from -pi/2 to pi/2
        xz = radius * cos(stackAngle);
        y = radius * sin(stackAngle);

        // add (sectorCount+1) vertices per stack
        // the first and last vertices have same position and normal, but different tex coords
        for (int j = 0; j <= sectorCount; ++j)
        {
            sectorAngle = j * sectorStep;           // starting from 0 to 2pi

            // vertex position
            x = xz * cos(sectorAngle);             // r * cos(u) * cos(v)
            z = xz * sin(sectorAngle);             // r * cos(u) * sin(v)
            ms->addVertex(vertexCount, dvec3(x, y, z));

            vertexCount++;
        }
    }

    // indices
    //  k1--k1+1
    //  |  / |
    //  | /  |
    //  k2--k2+1
    int k1, k2;
    for (int i = 0; i < stackCount; ++i)
    {
        k1 = i * (sectorCount + 1);     // beginning of current stack
        k2 = k1 + sectorCount + 1;      // beginning of next stack

        for (int j = 0; j < sectorCount; ++j, ++k1, ++k2)
        {
            // 2 triangles per sector excluding 1st and last stacks
            if (i != 0)
                ms->addFace(k1, k2, k1 + 1);   // k1---k2---k1+1

            if (i != (stackCount - 1))
                ms->addFace(k1 + 1, k2, k2 + 1); // k1+1---k2---k2+1
        }
    }

    return ms;
}


int main(int argc, char** argv)
{
    // mesh<double>* mesh = create_sphere_mesh(1.0); //  generate_surface(0, 0, 0);
    mesh<double>* mesh = generate_surface(0, 0, 0);
    obj_file_dump("c:\\temp\\surface.obj", mesh);
    delete mesh;
    return 0;
    
}
