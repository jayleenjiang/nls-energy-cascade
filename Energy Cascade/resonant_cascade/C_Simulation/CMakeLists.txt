cmake_minimum_required(VERSION 3.10)
project(EnergyCascade LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Hints for Homebrew installs
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew" "/usr/local")

# Dependencies
find_package(Eigen3 REQUIRED)
find_package(FFTW3 REQUIRED)

# Executable first
add_executable(EnergyCascade main.cpp)

# --- OpenMP (try CMake package first) ---
find_package(OpenMP QUIET)
if (OpenMP_CXX_FOUND)
    target_link_libraries(EnergyCascade PRIVATE OpenMP::OpenMP_CXX)
elseif(APPLE)
    # Manual fallback for Apple Clang + Homebrew libomp
    target_compile_options(EnergyCascade PRIVATE -Xpreprocessor -fopenmp)
    target_include_directories(EnergyCascade PRIVATE /opt/homebrew/opt/libomp/include)
    target_link_directories(EnergyCascade PRIVATE /opt/homebrew/opt/libomp/lib)
    target_link_libraries(EnergyCascade PRIVATE omp)  # links libomp
else()
    message(WARNING "OpenMP not found — building without it.")
endif()

# Link Eigen
target_link_libraries(EnergyCascade PRIVATE Eigen3::Eigen)

# FFTW3 (robust)
if (TARGET FFTW3::fftw3)
    target_link_libraries(EnergyCascade PRIVATE FFTW3::fftw3)
else()
    target_include_directories(EnergyCascade PRIVATE ${FFTW3_INCLUDE_DIRS})
    target_link_directories(EnergyCascade PRIVATE ${FFTW3_LIBRARY_DIRS})
    target_link_libraries(EnergyCascade PRIVATE ${FFTW3_LIBRARIES})
endif()

# Compiler flags
target_compile_options(EnergyCascade PRIVATE -Ofast -march=native)
