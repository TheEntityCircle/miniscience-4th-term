cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(gmshdemo)


# Подключаем библиотеки

set(GMSH_ROOT "/usr/local/opt/gmsh") 
set(GMSH_INCLUDE_DIR "${GMSH_ROOT}/include")
set(GMSH_LIBRARY "${GMSH_ROOT}/lib/libgmsh.dylib")




add_executable(sample01 sample01.cpp)
add_executable(sample02 sample02.cpp)
add_executable(sample03 sample03.cpp)
add_executable(sample04 sample04.cpp)
add_executable(sample05 sample05.cpp)
add_executable(sample06 sample06.cpp)
add_executable(lab01_1 lab01_1.cpp)
add_executable(lab01_3 lab01_3.cpp)

# --- Подключаем gmsh, тут придется перебирать варианты 

# - если вам повезло, то сработает простой путь
#include_directories(${GMSH_INCLUDE_DIR})

# - если у вас Ubuntu, то возможно, что приехал немного кривой gmsh из репозитория
# - при этом достаточно указать системные пути, куда система положила нужные файлы
# - закомментируйте строку 14 и раскомментируйте эти шесть:
#add_library(gmsh SHARED IMPORTED)
#set_target_properties(gmsh PROPERTIES
#   IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libgmsh.so"
#   INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
#)
#set(GMSH_LIBRARIES /usr/local/opt/gmsh/include/gmsh.h)


# - если предыдущий вариант не сработал, собирайте gmsh из исходников
# - после этого в предыдущем варианте заменяйте системные пути на те, что вы сделали сами

target_include_directories(sample01 PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(sample02 PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(sample03 PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(sample04 PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(sample05 PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(sample06 PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(lab01_1 PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(lab01_3 PRIVATE ${GMSH_INCLUDE_DIR})

# --- Закончили с gmsh, остальное править не обязательно (не поможет)

# Укажем библиотеки для линковки
target_link_libraries(sample01 PRIVATE ${GMSH_LIBRARY})
target_link_libraries(sample02 PRIVATE ${GMSH_LIBRARY})
target_link_libraries(sample03 PRIVATE ${GMSH_LIBRARY})
target_link_libraries(sample04 PRIVATE ${GMSH_LIBRARY})
target_link_libraries(sample05 PRIVATE ${GMSH_LIBRARY})
target_link_libraries(sample06 PRIVATE ${GMSH_LIBRARY})
target_link_libraries(lab01_1 PRIVATE ${GMSH_LIBRARY})
target_link_libraries(lab01_3 PRIVATE ${GMSH_LIBRARY})
