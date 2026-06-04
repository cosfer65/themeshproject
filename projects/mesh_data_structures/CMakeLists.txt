set(NAME "mesh_data_structures")

message(STATUS "-----> Configuring project: ${NAME}")

add_executable(${NAME} main.cpp)

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
