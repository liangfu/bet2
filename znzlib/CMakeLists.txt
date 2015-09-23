include_directories(${VTK_INCLUDE_DIRS})

INCLUDE_DIRECTORIES(. ..)
add_definitions(-DHAVE_ZLIB)
ADD_LIBRARY(znzlib STATIC znzlib.c)
target_link_libraries(znzlib vtkzlib)

# IF(FSL_BUILD)
#   INSTALL_FILES(/include/znzlib .*\\.h$)
#   INSTALL_TARGETS(/lib znz)
# ENDIF(FSL_BUILD)
