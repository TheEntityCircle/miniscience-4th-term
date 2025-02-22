cmake_minimum_required(VERSION 3.10)
project(my_project)
set(CMAKE_CXX_STANDARD 11)
# Пути к Gmsh (замените на ваш путь из `brew --prefix gmsh`)


set(GMSH_ROOT "/usr/local/opt/gmsh") 
set(GMSH_INCLUDE_DIR "${GMSH_ROOT}/include")
set(GMSH_LIBRARY "${GMSH_ROOT}/lib/libgmsh.dylib")

# Создание исполняемого файла
add_executable(my_app src/main.cpp)

# Подключение заголовков и библиотеки
target_include_directories(my_app PRIVATE ${GMSH_INCLUDE_DIR})
target_link_libraries(my_app PRIVATE ${GMSH_LIBRARY})