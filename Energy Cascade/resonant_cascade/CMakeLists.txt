cmake_minimum_required(VERSION 3.10)
project(ResonantCascade)

set(CMAKE_CXX_STANDARD 17)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Add executable
add_executable(resonant_cascade resonant_cascade.cpp)

# Link Eigen
target_link_libraries(resonant_cascade Eigen3::Eigen)

# Add optimization flags
target_compile_options(resonant_cascade PRIVATE -O3 -march=native)