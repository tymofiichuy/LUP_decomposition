cmake_minimum_required(VERSION 3.5)
project(eq_system)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(MAIN_DIR ${CMAKE_SOURCE_DIR}/main)

include_directories(${SOURCE_DIR})

add_executable(main ${MAIN_DIR}/main.cpp ${SOURCE_DIR}/eq_system.cpp)