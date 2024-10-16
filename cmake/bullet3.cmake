

if (TARGET bullet3::bullet3)
    return()
endif()

include(CPM)

CPMAddPackage(
    NAME bullet3
    GITHUB_REPOSITORY bulletphysics/bullet3
    GIT_TAG master
)

# Find_Package(Bullet)# CONFIG)
MESSAGE(STATUS "bullet3 SOURCE DIR: ${bullet3_SOURCE_DIR}")
MESSAGE(STATUS "bullet3 INCLUDE DIR: ${bullet3_INCLUDE_DIR}")
MESSAGE(STATUS "bullet3 LIBRARY DIR: ${bullet3_LIBRARY_DIRS}")
MESSAGE(STATUS "bullet3 ROOT DIR: ${bullet3_DIR}")

# set(BT_USE_DOUBLE_PRECISION ON)
set_target_properties(BulletDynamics PROPERTIES FOLDER third_party)
add_library(bullet3_ INTERFACE)
target_include_directories(bullet3_ INTERFACE ${bullet3_SOURCE_DIR}/src ${bullet3_SOURCE_DIR}/src/BulletCollision ${bullet3_SOURCE_DIR}/src/BulletDynamics ${bullet3_SOURCE_DIR}/src/LinearMath)
target_link_libraries(bullet3_ INTERFACE BulletDynamics BulletCollision LinearMath)
add_library(bullet3::bullet3 ALIAS bullet3_)


# set_target_properties(qhullcpp PROPERTIES FOLDER third_party)
# add_library(qhull_ INTERFACE)
# target_link_libraries(qhull_ INTERFACE qhullcpp qhullstatic_r)
# target_include_directories(qhull_ INTERFACE ${qhull_SOURCE_DIR}/src)
# add_library(qhull::qhull ALIAS qhull_)

# # # # # #
# Find_Package(Bullet CONFIG)
# SET(BULLET_INCLUDE_DIR ${Bullet_DIR}/${BULLET_ROOT_DIR}/${BULLET_INCLUDE_DIR})
# SET(BLA ${Bullet_DIR}/${BULLET_ROOT_DIR}/${BULLET_LIBRARY_DIRS})
# add_executable(bullet_sim "${bullet_SRCS}")

# target_include_directories(bullet_sim PUBLIC ${BULLET_INCLUDE_DIR})
# target_compile_definitions(bullet_sim PUBLIC ${BULLET_DEFINITIONS})
# target_link_directories(bullet_sim PUBLIC ${Bullet_DIR}/${BULLET_ROOT_DIR}/${BULLET_LIBRARY_DIRS})
# BulletDynamics BulletCollision LinearMath
# target_link_libraries(bullet_sim PUBLIC  geometry-central polyscope qhullcpp qhullstatic_r autodiff)
