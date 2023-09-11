cmake_minimum_required(VERSION 3.21)
project(Ising)
set(CMAKE_CXX_STANDARD 17)

include_directories(src/)

file(GLOB_RECURSE HEADERS
    src/*.h
    src/*.hpp)

file(GLOB_RECURSE SOURCES
    src/*.cpp
    src/*.cxx)

add_executable(${PROJECT_NAME}
    ${SOURCES}
    ${HEADERS})