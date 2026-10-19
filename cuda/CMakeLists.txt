cmake_minimum_required(VERSION 3.10)
project(ObjectDelineationCPP)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build a library for the object delineation core code.
add_library(object_delineation_lib
    src/ObjectDelineation.cpp
)
target_include_directories(object_delineation_lib PUBLIC ${PROJECT_SOURCE_DIR}/src)

# Build the main executable that uses the library.
add_executable(object_delineation
    src/main.cpp
)
target_link_libraries(object_delineation PRIVATE object_delineation_lib)

# Enable testing.
enable_testing()

# Download and configure GoogleTest.
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/release-1.12.1.zip
)
# For Windows: force shared CRT.
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Build the test executable.
add_executable(test_main_3 ${PROJECT_SOURCE_DIR}/tests/test_main_3.cpp)
target_link_libraries(test_main_3 PRIVATE gtest gtest_main object_delineation_lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")


# Register the tests with CTest.
add_test(NAME AllTests COMMAND test_main)

# Optional: enable GoogleTest auto-discovery.
include(GoogleTest)
gtest_discover_tests(test_main_3)
