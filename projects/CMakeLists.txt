#collect project directories
file(GLOB SUBDIRS
    LIST_DIRECTORIES true
    "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECTS_DIRECTORY}/*"
)

foreach(project ${SUBDIRS})
    if (IS_DIRECTORY ${project})
        # strip directory path to get project name
        get_filename_component(project_name ${project} NAME)
        message(STATUS "-----> Adding project: ${project_name}")
        add_subdirectory(${project_name})
    endif()
endforeach()
