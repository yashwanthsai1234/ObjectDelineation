# Builds the Python extension objdel._polycore (driven by scikit-build-core from pyproject.toml).
# The core itself is the header-only C++17 library in cpp/polycore, also used natively by its
# fuzz driver (make fuzz). Also builds objdel._gpucheck from cuda/ (below).
# References: https://nanobind.readthedocs.io/en/latest/building.html ,
#             https://scikit-build-core.readthedocs.io/en/latest/guide/cmakelists.html
cmake_minimum_required(VERSION 3.26...3.31)
project(objdel LANGUAGES CXX)

add_library(polycore INTERFACE)
target_include_directories(polycore INTERFACE cpp)
target_compile_features(polycore INTERFACE cxx_std_17)

find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT})
find_package(nanobind CONFIG REQUIRED)
find_package(Threads REQUIRED)

# NOMINSIZE: nanobind builds with -Os by default, which halves the kernels' speed (lab REPORT.md).
# STABLE_ABI: one cp312-abi3 wheel per platform (wheel.py-api in pyproject.toml).
nanobind_add_module(_polycore NB_STATIC NOMINSIZE STABLE_ABI src/objdel/_polycore.cpp)
target_link_libraries(_polycore PRIVATE polycore Threads::Threads)
install(TARGETS _polycore LIBRARY DESTINATION objdel)

# objdel._gpucheck: the batched checker of cuda/ (one source, cuda/gpucheck.cpp). Built for the CPU by the C++
# compiler unless OBJDEL_CUDA=ON (`make cuda` on a CUDA box), which compiles the same file with nvcc. The core is a
# static library so that nanobind's compiler flags never reach nvcc.
# References: https://cmake.org/cmake/help/latest/prop_tgt/CUDA_ARCHITECTURES.html ,
#             https://cmake.org/cmake/help/latest/module/FindCUDAToolkit.html
option(OBJDEL_CUDA "build objdel._gpucheck with nvcc" OFF)
add_library(gpucheck STATIC cuda/gpucheck.cpp)
target_link_libraries(gpucheck PUBLIC polycore)
set_target_properties(gpucheck PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(OBJDEL_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)  # the box's GPU; pass -DCMAKE_CUDA_ARCHITECTURES=... to cross-compile
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set_source_files_properties(cuda/gpucheck.cpp PROPERTIES LANGUAGE CUDA)
  # A target copies CMAKE_CUDA_ARCHITECTURES when it is created, before this block: set it on the target (R0 v3, the
  # box's default `native` left gpucheck with none; the docker check passed the archs on the command line)
  set_target_properties(gpucheck PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON
                                            CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
  target_link_libraries(gpucheck PUBLIC CUDA::cudart_static)  # the module needs only the driver at run time
endif()
nanobind_add_module(_gpucheck NB_STATIC NOMINSIZE STABLE_ABI cuda/gpucheck_py.cpp)
target_link_libraries(_gpucheck PRIVATE gpucheck)
install(TARGETS _gpucheck LIBRARY DESTINATION objdel)

# Editable installs also put the repo root on sys.path (what hatch's dev-mode-dirs did), so
# `python -I -m experiments.<phase>.<script>` works. Never part of a real wheel.
if(SKBUILD_STATE STREQUAL "editable")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/_objdel_repo.pth" "${CMAKE_CURRENT_SOURCE_DIR}\n")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/_objdel_repo.pth" DESTINATION .)
endif()
