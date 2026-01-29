#pragma once

#include "predefined_meshes.h"

namespace base_math {
    template <typename T>
    const static_mesh_def static_cube(
        {
            std::tuple<T, T, T>(0.5f, 0.5f, -0.5f),
            std::tuple<T, T, T>(0.5f, -0.5f, -0.5f),
            std::tuple<T, T, T>(0.5f, 0.5f, 0.5f),
            std::tuple<T, T, T>(0.5f, -0.5f, 0.5f),
            std::tuple<T, T, T>(-0.5f, 0.5f, -0.5f),
            std::tuple<T, T, T>(-0.5f, -0.5f, -0.5f),
            std::tuple<T, T, T>(-0.5f, 0.5f, 0.5f),
            std::tuple<T, T, T>(-0.5f, -0.5f, 0.5f)
        },
        {
            std::tuple<size_t, size_t, size_t>(5, 3, 1),
            std::tuple<size_t, size_t, size_t>(3, 8, 4),
            std::tuple<size_t, size_t, size_t>(7, 6, 8),
            std::tuple<size_t, size_t, size_t>(2, 8, 6),
            std::tuple<size_t, size_t, size_t>(1, 4, 2),
            std::tuple<size_t, size_t, size_t>(5, 2, 6),
            std::tuple<size_t, size_t, size_t>(5, 7, 3),
            std::tuple<size_t, size_t, size_t>(3, 7, 8),
            std::tuple<size_t, size_t, size_t>(7, 5, 6),
            std::tuple<size_t, size_t, size_t>(2, 4, 8),
            std::tuple<size_t, size_t, size_t>(1, 3, 4),
            std::tuple<size_t, size_t, size_t>(5, 1, 2)
        }
    );
}