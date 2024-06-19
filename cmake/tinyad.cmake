if (TARGET tinyad::TinyAD)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME TinyAD
    GITHUB_REPOSITORY patr-schm/TinyAD
    # VERSION 0
    GIT_TAG main
)


MESSAGE(STATUS "TINYAD SOURCE DIR: ${TinyAD_SOURCE_DIR}")

set_target_properties(TinyAD PROPERTIES FOLDER third_party)
add_library(tinyad_ INTERFACE)
target_link_libraries(tinyad_ INTERFACE TinyAD)
target_include_directories(tinyad_ INTERFACE ${TinyAD_SOURCE_DIR}/include/TinyAD)
add_library(tinyad::tinyad ALIAS tinyad_)