cmake_minimum_required(VERSION 3.20)
project(ducklips_voice CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(voice
  src/body.cpp
  src/geometry.cpp
  src/lexicon.cpp
  src/planner.cpp
  src/room.cpp
  src/score.cpp
  src/text.cpp
  src/wave_file.cpp)
target_include_directories(voice PUBLIC include)
target_compile_options(voice PRIVATE -Wall -Wextra -Wpedantic)

add_executable(vowel_probe tools/vowel_probe.cpp)
target_link_libraries(vowel_probe PRIVATE voice)

add_executable(tract_probe tools/tract_probe.cpp)
target_link_libraries(tract_probe PRIVATE voice)

add_executable(vowel_fit tools/vowel_fit.cpp)
target_link_libraries(vowel_fit PRIVATE voice)

add_executable(say tools/say.cpp)
target_link_libraries(say PRIVATE voice)

add_library(time_frequency analysis/time_frequency.cpp)
target_include_directories(time_frequency PUBLIC analysis)
target_compile_options(time_frequency PRIVATE -Wall -Wextra -Wpedantic)
add_executable(tf_analyze analysis/analyze.cpp)
target_link_libraries(tf_analyze PRIVATE time_frequency)

add_executable(tract_codebook tools/tract_codebook.cpp)
target_link_libraries(tract_codebook PRIVATE voice time_frequency)

add_executable(engine_codebook tools/engine_codebook.cpp)
target_link_libraries(engine_codebook PRIVATE voice time_frequency)

add_executable(transfer_probe tools/transfer_probe.cpp)
target_link_libraries(transfer_probe PRIVATE voice)

add_executable(lateral_fit tools/lateral_fit.cpp)
target_link_libraries(lateral_fit PRIVATE voice)
add_executable(lateral_profile tools/lateral_profile.cpp)
target_link_libraries(lateral_profile PRIVATE voice)
