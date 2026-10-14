cmake_minimum_required(VERSION 3.24)
project(DuckLipsBoundary LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Only the parent pieces this prototype uses: WAV I/O and the noise source.
set(DUCKLIPS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
add_library(ducklips_parts STATIC
    ${DUCKLIPS_ROOT}/src/wav.cpp
    ${DUCKLIPS_ROOT}/src/sources.cpp
    ${DUCKLIPS_ROOT}/src/physics.cpp
)
target_include_directories(ducklips_parts PUBLIC ${DUCKLIPS_ROOT}/include)

# No fast-math, matching the parent.
add_executable(boundary_lab
    boundary_lab.cpp
    geometry.cpp
    transport.cpp
    tube.cpp
    analysis.cpp
)
target_link_libraries(boundary_lab PRIVATE ducklips_parts)
target_compile_options(boundary_lab PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)
