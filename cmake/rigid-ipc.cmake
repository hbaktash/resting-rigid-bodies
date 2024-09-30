if (TARGET rigid-ipc::rigid-ipc)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME rigid-ipc
    GITHUB_REPOSITORY ipc-sim/rigid-ipc
    GIT_TAG v1.0.0 #master
)

MESSAGE(STATUS "rigid-ipc SOURCE DIR: ${rigid-ipc_SOURCE_DIR}")

set_target_properties(ipc_rigid PROPERTIES FOLDER third_party)
add_library(rigid-ipc_ INTERFACE)
target_link_libraries(rigid-ipc_ INTERFACE ipc_rigid)
# target_include_directories(rigid-ipc_ INTERFACE ${rigid-ipc_SOURCE_DIR}/include)
add_library(rigid-ipc::rigid-ipc ALIAS rigid-ipc_)