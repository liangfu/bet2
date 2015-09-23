cmake_minimum_required(VERSION 2.8.7)
PROJECT(newimage)

#add_subdirectory(../niftiio )

INCLUDE_DIRECTORIES(.  .. ../newmat)

set(Srcs complexvolume.cc
	costfns.cc
	generalio.cc
	imfft.cc
	lazy.cc
	newimage.cc
	newimagefns.cc
	skeletonapp.cc)

add_library(newimage ${Srcs})

add_definitions(-DEXPOSE_TREACHEROUS)

target_link_libraries(newimage newmat niftiio utils fslio)
