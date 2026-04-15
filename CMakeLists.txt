cmake_minimum_required(VERSION 3.20)
project(themeshproject LANGUAGES C CXX RC)

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
        ${CMAKE_CURRENT_SOURCE_DIR}/src/themeshproject.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/themeshview.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/visual_objects.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/curvature.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/themeshframe.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/themeshmodel.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/themeshproject.rc
)

target_link_libraries(meshproject
    PUBLIC btm-framework
)

if(MSVC)
    target_link_options(meshproject PRIVATE /SUBSYSTEM:WINDOWS)
endif()
