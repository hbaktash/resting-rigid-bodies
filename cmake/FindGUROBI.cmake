# Once done this will define
#  GUROBI_FOUND - System has Gurobi
#  GUROBI_INCLUDE_DIRS - The Gurobi include directories
#  GUROBI_LIBRARIES - The libraries needed to use Gurobi
# modified from https://gitlab.vci.rwth-aachen.de:9000/CoMISo/CoMISo/-/blob/4da5ab5deb61a538634ea719f8a1c920ab54e429/cmake/FindGUROBI.cmake

set(GUROBI_ENABLE ON CACHE BOOL "Enable gurobi?")
if (GUROBI_ENABLE)
set(GUROBI_HOME "/Library/gurobi1100/macos_universal2")
set(GUROBI_BASE $ENV{GUROBI_HOME} CACHE PATH "GUROBI root directory.")
find_path(GUROBI_INCLUDE_DIR
          NAMES gurobi_c++.h
          PATHS
          "${GUROBI_BASE}/include"
          "$ENV{GUROBI_HOME}/include"
          "/Library/gurobi1100/macos_universal2/include"
          )
message(STATUS "GUROBI_INCLUDE_DIR: ${GUROBI_INCLUDE_DIR}")
get_filename_component(GUROBI_LIB_DIR "${GUROBI_INCLUDE_DIR}/../lib" ABSOLUTE)
# GUROBI_BIN_DIR is needed on windows, where it contains the .dll
get_filename_component(GUROBI_BIN_DIR "${GUROBI_INCLUDE_DIR}/../bin" ABSOLUTE)
get_filename_component(GUROBI_SRC_DIR "${GUROBI_INCLUDE_DIR}/../src" ABSOLUTE)
file(GLOB GUROBI_LIBRARY_LIST
    ${GUROBI_LIB_DIR}/libgurobi*.dylib
    ${GUROBI_LIB_DIR}/libgurobi*.so)
# Ignore libgurobiXY_light.so, libgurobi.so (without version):
string(REGEX MATCHALL
    "libgurobi([0-9]+)\\..*"
    GUROBI_LIBRARY_LIST
    "${GUROBI_LIBRARY_LIST}"
    )
string(REGEX REPLACE
    "libgurobi([0-9]+)\\..*"
    "\\1"
    GUROBI_LIBRARY_VERSIONS
    "${GUROBI_LIBRARY_LIST}")
list(LENGTH GUROBI_LIBRARY_VERSIONS GUROBI_NUMVER)
message("GUROBI LIB VERSIONS: ${GUROBI_LIBRARY_VERSIONS}")
if (GUROBI_NUMVER EQUAL 1)
    list(GET GUROBI_LIBRARY_VERSIONS 0 GUROBI_LIBRARY_VERSION)
    set(GUROBI_LIBRARY_NAME "gurobi${GUROBI_LIBRARY_VERSION}")
else()
    # none or more than one versioned library -let's try without suffix,
    # maybe the user added a symlink to the desired library
    set(GUROBI_LIBRARY_NAME "gurobi")
endif()
#message("GUROBI LIB NAME: ${GUROBI_LIBRARY_NAME}")
find_library(GUROBI_LIBRARY
    NAMES ${GUROBI_LIBRARY_NAME}
    PATHS
    "${GUROBI_BASE}/lib"
    "$ENV{GUROBI_HOME}/lib"
)
# Gurobi ships with some compiled versions of its C++ library for specific
# compilers, however it also comes with the source code. We will compile
# the source code outselves -- this is much safer, as it guarantees the same
# compiler version and flags.
# (Note: doing this is motivated by actual sometimes-subtle ABI compatibility bugs)
# The old behaviour can be enabled with GUROBI_USE_PRECOMPILED_CXX)
set(GUROBI_USE_PRECOMPILED_CXX ON)
option(GUROBI_USE_PRECOMPILED_CXX "Use precompiled C++ libraries instead of building it ourselves. Not recommended." ON)
mark_as_advanced(GUROBI_USE_PRECOMPILED_CXX)
if(GUROBI_USE_PRECOMPILED_CXX)
  if ( CMAKE_GENERATOR MATCHES "^Visual Studio 12.*Win64" )
    SET(GUROBI_LIB_NAME "gurobi_c++md2013")
  endif()
  
find_library(GUROBI_CXX_LIBRARY 
             NAMES gurobi_c++
                   ${GUROBI_LIB_NAME}
             PATHS "$ENV{GUROBI_HOME}/lib" 
              "${GUROBI_BASE}/lib"
             )
message("GUROBI LIB name: ${GUROBI_LIB_NAME}")
else()
    file(GLOB GUROBI_CXX_SRC CONFIGURE_DEPENDS ${GUROBI_SRC_DIR}/cpp/*.cpp)
    if(NOT GUROBI_CXX_SRC)
        message(FATAL_ERROR "could not find gurobi c++ sources in GUROBI_SRC_DIR=${GUROBI_SRC_DIR}/cpp/.")
    endif()
    add_library(GurobiCXX STATIC ${GUROBI_CXX_SRC})
    target_include_directories(GurobiCXX PUBLIC ${GUROBI_INCLUDE_DIR})
    # We need to be able to link this into a shared library:
    set_target_properties(GurobiCXX PROPERTIES POSITION_INDEPENDENT_CODE ON)
    set(GUROBI_CXX_LIBRARY $<TARGET_FILE:GurobiCXX>)
endif(GUROBI_USE_PRECOMPILED_CXX)
set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}" )
# use c++ headers as default
# set(GUROBI_COMPILER_FLAGS "-DIL_STD" CACHE STRING "Gurobi Compiler Flags")
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI  DEFAULT_MSG
                                  GUROBI_CXX_LIBRARY GUROBI_LIBRARY GUROBI_INCLUDE_DIR)
mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_CXX_LIBRARY GUROBI_BIN_DIR )
endif(GUROBI_ENABLE)


# set(GUROBI_CUSTOM /Library/gurobi1100/macos_universal2) # replace with path to Gurobi on your computer

# find_path(GUROBI_INCLUDE_DIRS
#     NAMES gurobi_c++.h
#     HINTS ${GUROBI_CUSTOM}/include ${GUROBI_HOME} ${GUROBI_DIR} $ENV{GUROBI_HOME}
#     PATH_SUFFIXES include)

# find_library(GUROBI_LIBRARY
#     NAMES gurobi gurobi110
#     HINTS ${GUROBI_CUSTOM}/lib ${GUROBI_DIR} $ENV{GUROBI_HOME}
#     PATH_SUFFIXES lib)


# if(MSVC)
#     set(MSVC_YEAR "2017")
    
#     if(MT)
#         set(M_FLAG "mt")
#     else()
#         set(M_FLAG "md")
#     endif()
    
#     find_library(GUROBI_CXX_LIBRARY
#         NAMES gurobi_c++${M_FLAG}${MSVC_YEAR}
#         HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
#         PATH_SUFFIXES lib)
#     find_library(GUROBI_CXX_DEBUG_LIBRARY
#         NAMES gurobi_c++${M_FLAG}d${MSVC_YEAR}
#         HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
#         PATH_SUFFIXES lib)
# else()
#     find_library(GUROBI_CXX_LIBRARY
#         NAMES gurobi_c++ libgurobi_c++.a
#         HINTS ${GUROBI_CUSTOM}/lib ${GUROBI_DIR} ${GUROBI_HOME} $ENV{GUROBI_HOME}
#         PATH_SUFFIXES lib)
#     message(STATUS "GUROBI_CXX_LIBRARY: ${GUROBI_CXX_LIBRARY}")
#     set(GUROBI_CXX_DEBUG_LIBRARY ${GUROBI_CXX_LIBRARY})
# endif()

# include(FindPackageHandleStandardArgs)
# find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY)