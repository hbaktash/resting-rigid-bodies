cmake_minimum_required(VERSION 3.14.0)

project(rolling_dragons)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)  


# Print the build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

### Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(CMAKE_CXX_STANDARD 17)
  SET(BASE_CXX_FLAGS "-std=c++17 -Wall -Wextra")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 -fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()
  
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

# == Deps
add_subdirectory(deps/geometry-central)
add_subdirectory(deps/polyscope)
add_subdirectory(deps/qhull)
# add_subdirectory(deps/pymesh)
# add_subdirectory(deps/bullet3)

# SET(GUROBI_DIR "/Library/gurobi1100/macos_universal2")
find_package(GUROBI)
if (GUROBI_FOUND)
  message(STATUS "Gurobi found at ${GUROBI_USE_PRECOMPILED_CXX}")
else ()
  message(STATUS "Gurobi not found!")
endif ()
# include_directories(${GUROBI_INCLUDE_DIRS})

find_package(MOSEK)
if (MOSEK_FOUND)
  message(STATUS "MOSEK found at ${MOSEK_LIBRARY_DIR}")
else ()
  message(STATUS "MOSEK not found!")
endif ()


# == Build our project stuff

set(polygons_SRCS 
  src/main.cpp
  src/forward.cpp
	# add any other source files here
)

set(polyhedra_SRCS 
  src/main3D.cpp
  src/arc_algebra.cpp
  src/forward3D.cpp
  src/mesh_factory.cpp
  src/geometry_utils.cpp
  src/coloring.cpp
  src/visual_utils.cpp
  src/markov_model.cpp
  src/boundary_tools.cpp
  src/convex_hull.cpp
  src/utils.cpp
	# add any other source files here
)

set(bullet_SRCS 
  src/main_bullet.cpp
  src/forward3D.cpp
  src/mesh_factory.cpp
  src/geometry_utils.cpp
  src/bullet_sim.cpp
  src/coloring.cpp
  src/visual_utils.cpp
  src/arc_algebra.cpp
  src/boundary_tools.cpp
  src/convex_hull.cpp
  src/utils.cpp
  # add any other source files here
)

set(inverse_SRCS 
  src/main_inverse.cpp
  src/arc_algebra.cpp
  src/forward3D.cpp
  src/mesh_factory.cpp
  src/geometry_utils.cpp
  src/coloring.cpp
  src/visual_utils.cpp
  src/markov_model.cpp
  src/boundary_tools.cpp
  src/inv_design.cpp
  src/convex_hull.cpp
  src/deformation.cpp
  src/remesh_tools.cpp
  src/utils.cpp
	# add any other source files here
)

set(convex_fill 
  src/main_convex_fill.cpp
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
  src/remesh_tools.cpp
  src/optimization.cpp
  src/utils.cpp
	# add any other source files here
)

# set(paper_figs 
#   src/main_roll_figs.cpp
#   src/arc_algebra.cpp
#   src/forward3D.cpp
#   src/mesh_factory.cpp
#   src/geometry_utils.cpp
#   src/coloring.cpp
#   src/visual_utils.cpp
#   src/boundary_tools.cpp
#   src/inv_design.cpp
#   src/convex_hull.cpp
#   src/deformation.cpp
#   src/remesh_tools.cpp

# 	# add any other source files here
# )

# To change the name of your executable, change "gc_project" in the lines below to whatever you want
add_executable(polygons "${polygons_SRCS}")
target_include_directories(polygons PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# add the args.hxx project which we use for command line args
target_include_directories(polygons PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")
# target_include_directories(polygons PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include")
# target_include_directories(polygons PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/CGAL/include")
# target_include_directories(polygons PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/boost_1_82_0")
# target_compile_features(LIBRARY_NAME PUBLIC cxx_std_14)
target_link_libraries(polygons geometry-central polyscope)


# To change the name of your executable, change "gc_project" in the lines below to whatever you want
add_executable(polyhedra "${polyhedra_SRCS}")
target_include_directories(polyhedra PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/src")
target_include_directories(polyhedra PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# add the args.hxx project which we use for command line args
target_include_directories(polyhedra PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")
target_link_libraries(polyhedra PUBLIC geometry-central polyscope qhullcpp qhullstatic_r)


# To change the name of your executable, change "gc_project" in the lines below to whatever you want
# SET(qhull_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/build")
# Find_Package(qhull)
add_executable(invDesign "${inverse_SRCS}")
# target_link_directories(invDesign PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/src/libqhullcpp")
# target_include_directories(invDesign PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/src/libqhullcpp")
target_include_directories(invDesign PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/src")
target_include_directories(invDesign PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# add the args.hxx project which we use for command line args
target_include_directories(invDesign PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")
target_include_directories(invDesign PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include")
target_include_directories(invDesign PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/TinyAD/include")
target_include_directories(invDesign PUBLIC "${MOSEK_INCLUDE_DIRS}")
target_link_directories(invDesign PUBLIC ${MOSEK_LIBRARY_DIR})
target_link_libraries(invDesign PUBLIC geometry-central polyscope qhullcpp qhullstatic_r ) #qhull_r libqhull_r

# To change the name of your executable, change "gc_project" in the lines below to whatever you want
# SET(qhull_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/build")
# Find_Package(qhull)
add_executable(convexFill "${convex_fill}")
# target_link_directories(invDesign PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/src/libqhullcpp")
# target_include_directories(invDesign PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/src/libqhullcpp")
target_include_directories(convexFill PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/src")
target_include_directories(convexFill PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# add the args.hxx project which we use for command line args
target_include_directories(convexFill PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")
target_include_directories(convexFill PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include")
target_include_directories(convexFill PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/TinyAD/include")
# target_include_directories(convexFill PUBLIC "${MOSEK_INCLUDE_DIRS}")
# target_link_directories(convexFill PUBLIC ${MOSEK_LIBRARY_DIR})
target_include_directories(convexFill PUBLIC "${GUROBI_INCLUDE_DIRS}")
# target_link_libraries(convexFill PUBLIC ${GUROBI_CXX_LIBRARY} ${GUROBI_CXX_DEBUG_LIBRARY})
# target_link_libraries(convexFill PUBLIC ${GUROBI_LIBRARIES})
target_link_libraries(convexFill PUBLIC geometry-central polyscope qhullcpp qhullstatic_r ${GUROBI_LIBRARIES}) # ${MOSEK_LIBRARIES} qhull_r libqhull_r libfusion64.10.1.dylib libmosek64.10.1.dylib
# # paper figure stuff
# add_executable(paper_figs "${paper_figs}")
# target_include_directories(paper_figs PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/src")
# target_include_directories(paper_figs PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# # add the args.hxx project which we use for command line args
# target_include_directories(paper_figs PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")
# target_include_directories(paper_figs PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include")
# target_include_directories(paper_figs PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/TinyAD/include")
# target_link_libraries(paper_figs PUBLIC geometry-central polyscope qhullcpp qhullstatic_r ) #qhull_r libqhull_r


# Bullet stuff

SET(Bullet_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/bullet3/build_cmake")
Find_Package(Bullet CONFIG)
SET(BULLET_INCLUDE_DIR ${Bullet_DIR}/${BULLET_ROOT_DIR}/${BULLET_INCLUDE_DIR})
SET(BLA ${Bullet_DIR}/${BULLET_ROOT_DIR}/${BULLET_LIBRARY_DIRS})
add_executable(bullet_sim "${bullet_SRCS}")

target_include_directories(bullet_sim PUBLIC ${BULLET_INCLUDE_DIR})
target_compile_definitions(bullet_sim PUBLIC ${BULLET_DEFINITIONS})
target_link_directories(bullet_sim PUBLIC ${Bullet_DIR}/${BULLET_ROOT_DIR}/${BULLET_LIBRARY_DIRS})
target_include_directories(bullet_sim PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/qhull/src")
target_include_directories(bullet_sim PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# add the args.hxx project which we use for command line args
target_include_directories(bullet_sim PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")

target_link_libraries(bullet_sim PUBLIC BulletDynamics BulletCollision LinearMath geometry-central polyscope qhullcpp qhullstatic_r)
# target_link_libraries(bullet_sim )
