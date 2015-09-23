cmake_minimum_required(VERSION 2.8.7)
PROJECT(fslio)

# add_subdirectory(znz)

INCLUDE_DIRECTORIES(.  .. )

set(Srcs fslio.c)

add_library(fslio ${Srcs})

