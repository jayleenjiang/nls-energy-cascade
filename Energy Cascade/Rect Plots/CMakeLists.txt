# Specify the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Set the project name
project(nls_simulation)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the Boost libraries, specifically requiring the headers
find_package(Boost REQUIRED)

# Add the Boost include directories to the project
include_directories(${Boost_INCLUDE_DIRS})

# Add the executable target from our main source file
add_executable(nls_simulation main.cpp)

# Note: Boost::rational is a header-only library, so linking is not
# strictly necessary, but it's good practice if you use other Boost components.
# target_link_libraries(nls_simulation PRIVATE ${Boost_LIBRARIES})
