cmake_minimum_required(VERSION 3.12)
project(SelfDualPotBKZ)

# set compiler to g++
if(NOT CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER g++)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add compile options
add_compile_options(-Ofast -fopenmp -mtune=native -march=native -mfpmath=both -funroll-loops)

# include directory
include_directories(include)

# collect source files
file(GLOB_RECURSE LIB_SOURCES "src/*.cpp")

# OpenMP
find_package(OpenMP REQUIRED)

add_library(SDPotBKZ SHARED ${LIB_SOURCES})
find_library(NTL_LIBRARY NAMES ntl)
if(NOT NTL_LIBRARY)
    message(FATAL_ERROR "NTL library not found")
endif()

target_link_libraries(SDPotBKZ PRIVATE ${NTL_LIBRARY} OpenMP::OpenMP_CXX)

set_target_properties(
    SDPotBKZ PROPERTIES OUTPUT_NAME "SDPotBKZ"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
