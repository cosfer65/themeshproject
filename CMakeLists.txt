cmake_minimum_required(VERSION 3.20)
project(meshproject LANGUAGES C CXX RC)

# Add the framework submodule
add_subdirectory(external/btm-framework)

# add projects
add_subdirectory(projects)
