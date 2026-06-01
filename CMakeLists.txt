cmake_minimum_required(VERSION 3.20)
project(meshproject LANGUAGES C CXX RC)

# Add the framework submodule
add_subdirectory(external/btm-framework)

# Main executable
add_executable(meshproject)

target_compile_features(meshproject PUBLIC cxx_std_20)
target_compile_features(meshproject PUBLIC c_std_11)


set_target_properties(meshproject PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/runtime"
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/runtime"
)
target_include_directories(meshproject
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_sources(meshproject
    PRIVATE
        # source files
        ${CMAKE_CURRENT_SOURCE_DIR}/src/themeshproject.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/themeshview.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/visual_objects.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/curvature.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/themeshframe.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/themeshmodel.cpp
        # header files
        ${CMAKE_CURRENT_SOURCE_DIR}/include/resource.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/themeshframe.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/themeshmodel.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/themeshproject.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/themeshview.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/visual_objects.h
        # resource files
        ${CMAKE_CURRENT_SOURCE_DIR}/src/themeshproject.rc
)

target_link_libraries(meshproject
    PUBLIC btm-framework
)

add_definitions(-DGLEW_STATIC)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
add_definitions(-DWIN32_LEAN_AND_MEAN)

set_target_properties(meshproject PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/runtime"
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/runtime"
)


if(MSVC)
    target_link_options(meshproject PRIVATE /SUBSYSTEM:WINDOWS)
    target_link_libraries(${PROJECT_NAME} PRIVATE opengl32)
    target_link_libraries(${PROJECT_NAME} PRIVATE glu32.lib)
endif()
