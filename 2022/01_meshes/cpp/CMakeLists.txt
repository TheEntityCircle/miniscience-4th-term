cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(gmshdemo)
add_executable(sample01 sample01.cpp)
add_executable(sample02 sample02.cpp)
add_executable(sample03 sample03.cpp)
add_executable(sample04 sample04.cpp)

# gmsh укажем явно (потому что он был просто распакован руками)
add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/home/anganar/gmsh-4.9.4-Linux64-sdk/lib/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/home/anganar/gmsh-4.9.4-Linux64-sdk/include"
)
set(GMSH_LIBRARY gmsh)

# Укажем библиотеки для линковки
target_link_libraries(sample01 ${GMSH_LIBRARY})
target_link_libraries(sample02 ${GMSH_LIBRARY})
target_link_libraries(sample03 ${GMSH_LIBRARY})
target_link_libraries(sample04 ${GMSH_LIBRARY})
