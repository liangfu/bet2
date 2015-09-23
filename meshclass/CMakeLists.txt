cmake_minimum_required(VERSION 2.8.7)
PROJECT(meshclass)

#add_subdirectory(../niftiio )

INCLUDE_DIRECTORIES(.  .. ../newmat)

set(Srcs drawmesh.cpp
mesh.cpp
mpoint.cpp
point.cpp
profile.cpp
pt_special.cpp
selfintersection.cpp
triangle.cpp)

add_library(meshclass ${Srcs})

# add_definitions(-DEXPOSE_TREACHEROUS)

target_link_libraries(meshclass newmat niftiio)
