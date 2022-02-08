cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME demo)
project(${PROJECT_NAME})

# Get DOLFIN configuration data (DOLFINConfig.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

add_executable(poisson Poisson.cpp)
target_link_libraries(poisson dolfin)

add_executable(hyper_elasticity HyperElasticity.cpp)
target_link_libraries(hyper_elasticity dolfin)
