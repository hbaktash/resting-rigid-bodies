cmake_minimum_required(VERSION 3.16.0)

project(rolling_dragons)

set(CMAKE_CXX_STANDARD 20)

# (Optional) toggle heavy deps later if desired
option(ROLLING_DRAGONS_USE_BULLET "Enable Bullet physics features" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# Dependencies (as before)
include(Eigen3)
include(geometry-central)
include(polyscope)
include(qhull)
include(tinyAD)
include(args)
include(stan)
include(osqp-cpp)
# include(autodiff)
# include(ipc-toolkit)
include(json)
if(ROLLING_DRAGONS_USE_BULLET)
  include(bullet3)
endif()

# Core library sources (exclude any old main_* files)
set(RD_CORE_SOURCES
    src/utils.cpp
    src/arc_algebra.cpp
    src/forward3D.cpp
    src/mesh_factory.cpp
    src/geometry_utils.cpp
    src/coloring.cpp
    src/visual_utils.cpp
    src/file_IO.cpp
    src/ambient_conversions.cpp
    src/implot.cpp
    src/implot_items.cpp
    src/boundary_tools.cpp
    src/inv_design.cpp        # remove later if not needed
    src/convex_hull.cpp
    src/remesh_tools.cpp
    src/bullet_sim.cpp        # keep only if bullet-dependent code is still referenced from main.cpp
)

if(NOT ROLLING_DRAGONS_USE_BULLET)
  list(REMOVE_ITEM RD_CORE_SOURCES src/bullet_sim.cpp)
endif()

add_library(rolling_dragons ${RD_CORE_SOURCES})
add_library(rolling_dragons::rolling_dragons ALIAS rolling_dragons)
target_include_directories(rolling_dragons PUBLIC include)
target_link_libraries(rolling_dragons
  PUBLIC
    geometry-central::geometry-central
    polyscope::polyscope
    qhull::qhull
    tinyad::tinyad
    osqp-cpp::osqp-cpp
    stan::stan
    nlohmann::json
)

if(ROLLING_DRAGONS_USE_BULLET)
  target_link_libraries(rolling_dragons PUBLIC bullet3::bullet3)
  target_compile_definitions(rolling_dragons PUBLIC ROLLING_DRAGONS_USE_BULLET)
endif()

# Single unified executable (assumes all previous functionality merged into src/main.cpp)
add_executable(main_app src/main.cpp)
target_link_libraries(main_app
  PRIVATE
    rolling_dragons::rolling_dragons
    taywee::args
)

set_property(TARGET main_app PROPERTY CXX_STANDARD 20)
