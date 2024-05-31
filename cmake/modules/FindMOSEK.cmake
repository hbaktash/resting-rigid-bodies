# - Try to find MOSEK

#  MOSEK_BASE - base path to the library

# Once done this will define
#  MOSEK_FOUND - System has Mosek
#  MOSEK_INCLUDE_DIRS - The Mosek include directories
#  MOSEK_LIBRARIES - The libraries needed to use Mosek
#  MOSEK_LIBRARY_DIR - The Mosek library dicrectory

set (MOSEK_BASE "c:" CACHE PATH "Base path of your gorubi installation")

if (MOSEK_INCLUDE_DIR)
  # in cache already
  set(MOSEK_FOUND TRUE)
  set(MOSEK_INCLUDE_DIRS "${MOSEK_INCLUDE_DIR}" )
  set(MOSEK_LIBRARIES "${MOSEK_LIBRARY};${MOSEK_CXX_LIBRARY}" )
else (MOSEK_INCLUDE_DIR)

  

find_path(MOSEK_INCLUDE_DIR 
          NAMES mosek.h
          PATHS "${CMAKE_SOURCE_DIR}/deps/mosek/10.1/tools/platform/osxaarch64/h"
          )

find_path(MOSEK_LIBRARY_DIR 
          NAMES libmosek64.dylib
          PATHS "${CMAKE_SOURCE_DIR}/deps/mosek/10.1/tools/platform/osxaarch64/bin"
          )

find_library( MOSEK_LIBRARY 
              NAMES mosek64

              PATHS "${CMAKE_SOURCE_DIR}/deps/mosek/10.1/tools/platform/osxaarch64/bin"
              )
  
find_library( MOSEK_CXX_LIBRARY 
              NAMES fusion64

              PATHS "${CMAKE_SOURCE_DIR}/deps/mosek/10.1/tools/platform/osxaarch64/bin"
              )
			  
# Binary dir for DLLs			
find_path(MOSEK_BIN_DIR 
                NAMES "mosek.dll" 
              PATHS "${CMAKE_SOURCE_DIR}/deps/mosek/10.1/tools/platform/osxaarch64/bin"

                DOC "Directory containing the MOSEK DLLs"
               ) 		  

set(MOSEK_INCLUDE_DIRS "${MOSEK_INCLUDE_DIR}")
set(MOSEK_LIBRARIES "${MOSEK_LIBRARY};${MOSEK_CXX_LIBRARY}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(MOSEK  DEFAULT_MSG
                                  MOSEK_LIBRARY MOSEK_CXX_LIBRARY MOSEK_INCLUDE_DIR)

mark_as_advanced(MOSEK_INCLUDE_DIR MOSEK_LIBRARY MOSEK_CXX_LIBRARY MOSEK_BIN_DIR MOSEK_LIBRARY_DIR)

endif(MOSEK_INCLUDE_DIR)
