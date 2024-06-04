if (TARGET qhull::qhull)
    return()
endif()

include(CPM)

set(QHULL_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_APPLICATIONS OFF CACHE BOOL "" FORCE)

CPMAddPackage(
    NAME qhull
    GITHUB_REPOSITORY qhull/qhull
    GIT_TAG v8.1-alpha1
)

MESSAGE(STATUS "QHULL SOURCE DIR: ${qhull_SOURCE_DIR}")

set_target_properties(qhullcpp PROPERTIES FOLDER third_party)
add_library(qhull_ INTERFACE)
target_link_libraries(qhull_ INTERFACE qhullcpp qhullstatic_r)
target_include_directories(qhull_ INTERFACE ${qhull_SOURCE_DIR}/src)
add_library(qhull::qhull ALIAS qhull_)
