cmake_minimum_required(VERSION 3.16.0)

project(Resting_Rigid_Bodies)

set(CMAKE_CXX_STANDARD 20)

# Feature toggles
option(RESTING_RIGID_BODIES_USE_BULLET "Enable Bullet physics features" OFF)
option(RESTING_RIGID_BODIES_USE_POLYSCOPE "Enable Polyscope visualization" ON)
option(RESTING_RIGID_BODIES_ENABLE_INVERSE_DESIGN "Build inverse design libs/apps" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
# Dependencies
include(Eigen3)
include(geometry-central)
include(qhull)
include(args)
include(json)

if(RESTING_RIGID_BODIES_USE_POLYSCOPE)
  include(polyscope)
endif()

if(RESTING_RIGID_BODIES_USE_BULLET)
  include(bullet3)
endif()

# Core library sources
set(RD_CORE_SOURCES
    src/utils.cpp
    src/arc_algebra.cpp
    src/forward3D.cpp
    src/mesh_factory.cpp
    src/geometry_utils.cpp
    src/file_IO.cpp
    src/ambient_conversions.cpp
    src/boundary_tools.cpp
    src/convex_hull.cpp
)

# Add viz sources only when Polyscope is enabled
if(RESTING_RIGID_BODIES_USE_POLYSCOPE)
  list(APPEND RD_CORE_SOURCES
    src/visual_utils.cpp
  )
endif()

# Add Bullet sources only when Bullet is enabled
if(RESTING_RIGID_BODIES_USE_BULLET)
  list(APPEND RD_CORE_SOURCES src/bullet_sim.cpp)
endif()

add_library(resting_rigid_bodies ${RD_CORE_SOURCES})
add_library(resting_rigid_bodies::resting_rigid_bodies ALIAS resting_rigid_bodies)
target_include_directories(resting_rigid_bodies PUBLIC include)
target_link_libraries(resting_rigid_bodies
  PUBLIC
    geometry-central::geometry-central
    qhull::qhull
    nlohmann::json
)


# Link optional deps
target_compile_definitions(resting_rigid_bodies PUBLIC RESTING_RIGID_BODIES_ENABLE_INVERSE_DESIGN)
if(RESTING_RIGID_BODIES_USE_POLYSCOPE)
  target_link_libraries(resting_rigid_bodies PUBLIC polyscope::polyscope)
  target_compile_definitions(resting_rigid_bodies PUBLIC RESTING_RIGID_BODIES_USE_POLYSCOPE)
endif()
if(RESTING_RIGID_BODIES_USE_BULLET)
  target_link_libraries(resting_rigid_bodies PUBLIC bullet3::bullet3)
  target_compile_definitions(resting_rigid_bodies PUBLIC RESTING_RIGID_BODIES_USE_BULLET)
endif()

# Existing drop/prob tool (headless-friendly)
add_executable(drop_probs src/main_drop_prob.cpp)
target_link_libraries(drop_probs
  PRIVATE
    resting_rigid_bodies::resting_rigid_bodies
    taywee::args
)
if(RESTING_RIGID_BODIES_USE_POLYSCOPE)
  target_compile_definitions(drop_probs PRIVATE RESTING_RIGID_BODIES_USE_POLYSCOPE)
endif()
if(RESTING_RIGID_BODIES_USE_BULLET)
  target_compile_definitions(drop_probs PRIVATE RESTING_RIGID_BODIES_USE_BULLET)
endif()
set_property(TARGET drop_probs PROPERTY CXX_STANDARD 20)


# == Inverse design (optional) ==
if(RESTING_RIGID_BODIES_ENABLE_INVERSE_DESIGN)
  add_library(resting_rigid_bodies_inverse
    src/inverse_design/deformation.cpp
    src/inverse_design/CP_tools.cpp
    src/inverse_design/remesh_tools.cpp
    src/inverse_design/optimization.cpp
    src/inverse_design/prob_assignment.cpp
  )
  add_library(resting_rigid_bodies::inverse ALIAS resting_rigid_bodies_inverse)

  target_include_directories(resting_rigid_bodies_inverse
    PUBLIC
      ${CMAKE_CURRENT_LIST_DIR}/include
      ${CMAKE_CURRENT_LIST_DIR}/include/inverse_design
  )

  # Make the macro visible when compiling inverse code
  target_compile_definitions(resting_rigid_bodies_inverse PUBLIC RESTING_RIGID_BODIES_ENABLE_INVERSE_DESIGN)

  # Optional solvers (only needed for inverse)
  include(tinyAD)
  include(osqp-cpp)
  include(stan)

  # GUROBI
  find_package(GUROBI REQUIRED)
  if(GUROBI_FOUND)
    message(STATUS "Using Gurobi at: ${GUROBI_INCLUDE_DIRS}")
    target_link_libraries(resting_rigid_bodies_inverse PRIVATE Gurobi::gurobi_cxx)
    target_compile_definitions(resting_rigid_bodies_inverse PUBLIC RESTING_RIGID_BODIES_HAVE_GUROBI)
  else()
    message(STATUS "Gurobi not found; building inverse design without it")
  endif()

  target_link_libraries(resting_rigid_bodies_inverse
    PUBLIC
      resting_rigid_bodies::resting_rigid_bodies
      stan::stan           # was PRIVATE
      tinyad::tinyad       # was PRIVATE
      osqp-cpp::osqp-cpp   # was PRIVATE
    PRIVATE
      ${GUROBI_LIBRARY}
      ${GUROBI_CXX_LIBRARY}
  )

  # Apps for inverse design (now under src/inverse_design)
  if(RESTING_RIGID_BODIES_USE_POLYSCOPE)
    add_executable(convex_dice src/main_convex_dice.cpp)
    target_link_libraries(convex_dice PRIVATE resting_rigid_bodies::inverse taywee::args)
    target_compile_definitions(convex_dice PRIVATE RESTING_RIGID_BODIES_USE_POLYSCOPE)
    # Not required, but harmless if you also add the macro here:
    # target_compile_definitions(convex_dice PRIVATE RESTING_RIGID_BODIES_ENABLE_INVERSE_DESIGN)
  endif()
endif()