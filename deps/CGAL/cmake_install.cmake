# Install script for directory: /Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "CGAL_Qt5" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/CGAL/Qt" TYPE DIRECTORY FILES "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/include/CGAL/Qt/")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/CGAL" TYPE FILE FILES
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/AUTHORS"
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/CHANGES.md"
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/LICENSE"
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/LICENSE.GPL"
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/LICENSE.LGPL"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/include/CGAL" REGEX "/\\.svn$" EXCLUDE REGEX "/qt$" EXCLUDE)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/include/CGAL" REGEX "/\\.svn$" EXCLUDE REGEX "/qt$" EXCLUDE)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/scripts/cgal_create_CMakeLists"
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/scripts/cgal_create_cmake_script"
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/scripts/cgal_make_macosx_app"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/CGAL" TYPE DIRECTORY FILES "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/cmake/modules/")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/CGAL" TYPE FILE FILES "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/cmake/modules/UseCGAL.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/CGAL" TYPE FILE FILES
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/lib/cmake/CGAL/CGALConfig.cmake"
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/lib/cmake/CGAL/CGALConfigBuildVersion.cmake"
    "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/lib/cmake/CGAL/CGALConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/man/man1" TYPE FILE FILES "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/auxiliary/cgal_create_cmake_script.1")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/hbaktash/Desktop/code/gc-polyscope-project-template/deps/CGAL/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
