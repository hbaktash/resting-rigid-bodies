if (TARGET osqp-cpp::osqp-cpp)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME osqp-cpp
    GITHUB_REPOSITORY google/osqp-cpp
    GIT_TAG master
)


MESSAGE(STATUS "osqp-cpp SOURCE DIR: ${osqp-cpp_SOURCE_DIR}")

# set_target_properties(osqp-cpp PROPERTIES FOLDER third_party)
# add_library(osqp-cpp_ INTERFACE)
# target_link_libraries(osqp-cpp_ INTERFACE absl::strings absl::status absl::statusor Eigen3::Eigen PRIVATE osqpstatic ${CMAKE_DL_LIBS})
# target_include_directories(osqp-cpp_ INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
# add_library(osqp-cpp::osqp-cpp ALIAS osqp-cpp_)
set_target_properties(osqp-cpp PROPERTIES FOLDER third_party)
add_library(osqp-cpp_ INTERFACE)
target_link_libraries(osqp-cpp_ INTERFACE osqp-cpp)
target_include_directories(osqp-cpp_ INTERFACE ${osqp-cpp_SOURCE_DIR}/include)
add_library(osqp-cpp::osqp-cpp ALIAS osqp-cpp_)