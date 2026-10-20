# Builds the Python extension objdel._polycore (driven by scikit-build-core from pyproject.toml).
# The core itself is the header-only C++17 library in cpp/polycore, also used natively by its
# fuzz driver (make fuzz) and later by cuda/.
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

# Editable installs also put the repo root on sys.path (what hatch's dev-mode-dirs did), so
# `python -I -m experiments.<phase>.<script>` works. Never part of a real wheel.
if(SKBUILD_STATE STREQUAL "editable")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/_objdel_repo.pth" "${CMAKE_CURRENT_SOURCE_DIR}\n")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/_objdel_repo.pth" DESTINATION .)
endif()
