cmake_minimum_required(VERSION 2.8.7)

project(bet2)

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_MODULE_PATH ${bnviewer_SOURCE_DIR}/CMake)

# this is EXTREMELY important to put these setting before add_library/executable
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 

SET(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

# find_package(VTK REQUIRED)
# include(${VTK_USE_FILE})
# include_directories(${VTK_INCLUDE_DIRS})

add_subdirectory(znzlib)
add_subdirectory(niftiio)
add_subdirectory(newmat)
add_subdirectory(fslio)
add_subdirectory(newimage)
add_subdirectory(meshclass)
add_subdirectory(miscmaths)
add_subdirectory(utils)
add_subdirectory(vtkzlib)

set( Srcs bet2.cpp)
##set( Hdrs bet2.cpp)

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
	newmat fslio newimage
)

add_executable(bet2 ${Srcs} )
target_link_libraries(bet2
	newmat fslio newimage niftiio znzlib meshclass utils miscmaths vtkzlib
	)


