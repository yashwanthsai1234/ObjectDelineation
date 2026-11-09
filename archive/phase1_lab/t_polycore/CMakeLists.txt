cmake_minimum_required(VERSION 3.18...3.31)
project(polycore LANGUAGES CXX)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(nanobind CONFIG REQUIRED)

# Header-only core (shared with a future CUDA engine).
add_library(polycore_core INTERFACE)
target_include_directories(polycore_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
target_compile_features(polycore_core INTERFACE cxx_std_17)

# NOMINSIZE: nanobind otherwise appends -Os, which de-vectorises the kernels (~2x slower).
nanobind_add_module(_polycore NB_STATIC NOMINSIZE src/_polycore.cpp)
target_link_libraries(_polycore PRIVATE polycore_core)
find_package(Threads REQUIRED)
target_link_libraries(_polycore PRIVATE Threads::Threads)

install(TARGETS _polycore LIBRARY DESTINATION polycore)

# Optional native test / fuzz driver (not part of the wheel).
option(POLYCORE_NATIVE_TESTS "Build the native fuzz driver" OFF)
if (POLYCORE_NATIVE_TESTS)
  add_executable(polycore_fuzz cpp/polycore/tests/fuzz_main.cpp)
  target_link_libraries(polycore_fuzz PRIVATE polycore_core Threads::Threads)
endif()
