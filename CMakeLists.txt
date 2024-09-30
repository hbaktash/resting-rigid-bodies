cmake_minimum_required(VERSION 3.14.0)

project(rolling_dragons)

# Add cmake recipe directory
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# Dependencies
include(autodiff)
include(geometry-central)
include(polyscope)
include(autodiff)
include(qhull)
include(tinyAD)
include(args)
include(stan)
include(osqp-cpp)
# include(ipc-toolkit)
# include(bullet3)

# Library source files
set(SRC_FILES
    src/utils.cpp
    src/arc_algebra.cpp
    src/forward3D.cpp
    src/mesh_factory.cpp
    src/geometry_utils.cpp
    src/coloring.cpp
    src/visual_utils.cpp
    src/implot.cpp
    src/implot_items.cpp
    src/boundary_tools.cpp
    src/inv_design.cpp
    src/convex_hull.cpp
    src/deformation.cpp
    src/optimization.cpp
    src/remesh_tools.cpp
    # src/bullet_sim.cpp
)

# Add rolling-dragons library, which will be shared by all executables.
add_library(rolling_dragons ${SRC_FILES})
target_include_directories(rolling_dragons PUBLIC include)
target_link_libraries(rolling_dragons
PUBLIC
geometry-central::geometry-central
polyscope::polyscope
autodiff::autodiff
qhull::qhull
tinyad::tinyad
osqp-cpp::osqp-cpp
)
add_library(rolling_dragons::rolling_dragons ALIAS rolling_dragons)

# Add executables

# add_executable(adDebug src/main_ad_debug.cpp)
# target_link_libraries(adDebug
# PRIVATE
# rolling_dragons::rolling_dragons
# taywee::args
# stan::stan
# )
# set_property(TARGET adDebug PROPERTY CXX_STANDARD 20)

# Add executables

add_executable(convexFill src/main_convex_fill.cpp)
target_link_libraries(convexFill
PRIVATE
rolling_dragons::rolling_dragons
taywee::args
stan::stan
)
set_property(TARGET convexFill PROPERTY CXX_STANDARD 20)

# Add executables
include(json)
# include(rigid-ipc)
add_executable(scriptIPC src/main_IPC_script.cpp)
target_link_libraries(scriptIPC
PRIVATE
rolling_dragons::rolling_dragons
taywee::args
nlohmann::json
# ipc::toolkit
)
set_property(TARGET scriptIPC PROPERTY CXX_STANDARD 20)

# Add executables

# add_executable(bulletSim src/main_bullet.cpp)
# target_link_libraries(bulletSim
# PRIVATE
# rolling_dragons::rolling_dragons
# taywee::args
# bullet3::bullet3
# )
# set_property(TARGET bulletSim PROPERTY CXX_STANDARD 20)

