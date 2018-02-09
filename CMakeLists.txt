project(BAL)
cmake_minimum_required(VERSION 2.8)

set( CMAKE_CXX_FLAGS "-std=c++11")

find_package(Ceres REQUIRED)

INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})

add_executable(BAL main.cpp)
target_link_libraries(BAL ${CERES_LIBRARIES})

