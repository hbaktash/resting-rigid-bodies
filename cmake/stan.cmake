if (TARGET stan::stan)
    return()
endif()

include(CPM)

set(ARGS_BUILD_EXAMPLE CACHE BOOL OFF FORCE)
set(ARGS_BUILD_UNITTESTS CACHE BOOL OFF FORCE)

CPMAddPackage(
    NAME stan
    GITHUB_REPOSITORY stan-dev/math
    GIT_TAG v4.9.0
    DOWNLOAD_ONLY YES
)

set(stan_BOOST_INCLUDE_DIR ${stan_SOURCE_DIR}/lib/boost_1.84.0/)
set(stan_SUNDIALS_INCLUDE_DIR ${stan_SOURCE_DIR}/lib/sundials_6.1.1/include/)

include(Eigen3)
#set(EIGEN_BUILD_DOC CACHE BOOL OFF FORCE)
#set(BUILD_TESTING CACHE BOOL OFF FORCE)
#add_subdirectory(${stan_SOURCE_DIR}/lib/eigen_3.4.0/)

include(tbb)
add_library(stan INTERFACE)
target_include_directories(stan INTERFACE ${stan_SOURCE_DIR})
target_include_directories(stan
    INTERFACE
    ${stan_BOOST_INCLUDE_DIR}
    ${stan_SUNDIALS_INCLUDE_DIR}
)
target_compile_definitions(stan INTERFACE "_REENTRANT")
target_link_libraries(stan INTERFACE "${stan_LIBRARIES}" TBB::tbb Eigen3::Eigen)
add_library(stan::stan ALIAS stan)
