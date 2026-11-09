# Builds the Python extension objdel._gpucheck (driven by scikit-build-core from pyproject.toml): the checker and
# scorer of cuda/ (one source, cuda/gpucheck.cpp). The CPU checker it replaced, polycore, is in archive/cpu_checker.
# References: https://nanobind.readthedocs.io/en/latest/building.html ,
#             https://scikit-build-core.readthedocs.io/en/latest/guide/cmakelists.html
cmake_minimum_required(VERSION 3.26...3.31)
project(objdel LANGUAGES CXX)

find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT})
find_package(nanobind CONFIG REQUIRED)

# Built for the CPU by the C++ compiler unless OBJDEL_CUDA=ON (`make cuda` on a CUDA box), which compiles the same
# file with nvcc. The core is a static library so that nanobind's compiler flags never reach nvcc.
# References: https://cmake.org/cmake/help/latest/prop_tgt/CUDA_ARCHITECTURES.html ,
#             https://cmake.org/cmake/help/latest/module/FindCUDAToolkit.html
option(OBJDEL_CUDA "build objdel._gpucheck with nvcc" OFF)
add_library(gpucheck STATIC cuda/gpucheck.cpp)
target_compile_features(gpucheck PUBLIC cxx_std_17)
set_target_properties(gpucheck PROPERTIES POSITION_INDEPENDENT_CODE ON)
# No fused multiply-add contraction: the scorer's float64 terms must be the IEEE operations of the host formulas,
# one rounding each, on the CPU build and the GPU build alike (gpucheck.cpp, scoring). Under nvcc the host functions
# go to the host compiler, whose C++ default contracts (g++ -ffp-contract=fast; an aarch64 host emits fmadd), so it
# gets the flag too (docs/implementation2.0/notes/research/cuda-kernels.md section 2).
target_compile_options(gpucheck PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>
                                        $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false -Xcompiler=-ffp-contract=off>)
# Test-only compile definitions (a CMake list, e.g. GC_POISON_SCAN;GC_ALIAS_SCAN_FOR_TEST): only tests/variants.py sets
# it, for its CPU test builds in build/variants/ (gpucheck.cpp, test_hooks()). Empty in every shipped build.
set(OBJDEL_TEST_HOOKS "" CACHE STRING "test-only compile definitions of gpucheck (tests/variants.py)")
if(OBJDEL_TEST_HOOKS)
  target_compile_definitions(gpucheck PRIVATE ${OBJDEL_TEST_HOOKS})
endif()
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
  # NVTX ranges (header-only nvtx3). FindCUDAToolkit defines CUDA::nvtx3 without checking that the headers exist, and
  # the PyPI toolkit has them only with the nvtx extra (nvidia/cu13/include/nvtx3): check the file itself.
  find_file(OBJDEL_NVTX3_HPP nvtx3/nvtx3.hpp PATHS ${CUDAToolkit_INCLUDE_DIRS} NO_DEFAULT_PATH)
  if(NOT OBJDEL_NVTX3_HPP)
    message(FATAL_ERROR "nvtx3 headers not found under ${CUDAToolkit_INCLUDE_DIRS}: install cuda-toolkit[nvtx]")
  endif()
  target_link_libraries(gpucheck PRIVATE CUDA::nvtx3)
endif()
# NOMINSIZE: nanobind builds with -Os by default, which halves the kernels' speed (lab REPORT.md).
# STABLE_ABI: one cp312-abi3 wheel per platform (wheel.py-api in pyproject.toml).
nanobind_add_module(_gpucheck NB_STATIC NOMINSIZE STABLE_ABI cuda/gpucheck_py.cpp)
target_link_libraries(_gpucheck PRIVATE gpucheck)
install(TARGETS _gpucheck LIBRARY DESTINATION objdel)

# Editable installs also put the repo root on sys.path (what hatch's dev-mode-dirs did), so
# `python -I -m experiments.<phase>.<script>` works. Never part of a real wheel.
if(SKBUILD_STATE STREQUAL "editable")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/_objdel_repo.pth" "${CMAKE_CURRENT_SOURCE_DIR}\n")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/_objdel_repo.pth" DESTINATION .)
endif()
