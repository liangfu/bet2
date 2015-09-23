

cmake_minimum_required(VERSION 2.8.7)
PROJECT(niftiio)

#add_subdirectory(../niftiio )

INCLUDE_DIRECTORIES(.  .. ../znzlib)

set(Srcs nifti1_io.c)

add_library(niftiio ${Srcs})

add_definitions(-DHAVE_ZLIB)

target_link_libraries(niftiio znzlib)
