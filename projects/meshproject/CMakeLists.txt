set(NAME "meshproject")

message(STATUS "-----> Configuring project: ${NAME}")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# message (STATUS "-----> CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
# message (STATUS "-----> CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
# message (STATUS "-----> CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")

# colect source files in src directory and its subdirectories
set(HEADER_FILES)
file(GLOB_RECURSE HEADER_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
)

set(SOURCE_FILES)
file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

set(RESOURCE_FILES)
file(GLOB_RECURSE RESOURCE_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.rc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.ico"
)
# foreach (RESOURCE_FILE ${RESOURCE_FILES})
#     message(STATUS "Resource file: ${RESOURCE_FILE}")
# endforeach()

set(SHADER_FILES)
set(SHADER_FILES
    "${CMAKE_SOURCE_DIR}/runtime/resources/shaders/meshprojectVertexShader.glsl"
    "${CMAKE_SOURCE_DIR}/runtime/resources/shaders/meshprojectFragmentShader.glsl"
)
foreach (SHADER_FILE ${SHADER_FILES})
    message(STATUS "Shader file: ${SHADER_FILE}")
endforeach()

add_executable(${NAME} ${HEADER_FILES} ${SOURCE_FILES} ${RESOURCE_FILES} ${SHADER_FILES})

# we group files in Visual Studio solution explorer
source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${HEADER_FILES})
source_group("Resource Files" FILES ${RESOURCE_FILES})
source_group("Shader Files" FILES ${SHADER_FILES})

set_property(
    SOURCE ${SHADER_FILES}
    PROPERTY VS_SETTINGS "ExcludedFromBuild=true"
)

target_include_directories(${NAME}
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(${NAME} PUBLIC btm-framework)

add_definitions(-DGLEW_STATIC)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
add_definitions(-DWIN32_LEAN_AND_MEAN)

if(MSVC)
    target_link_options(${NAME} PRIVATE "/SUBSYSTEM:WINDOWS")
else()
    target_link_options(${NAME} PRIVATE "-mwindows")
endif()

set_target_properties(${NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/runtime"
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/runtime"
)
