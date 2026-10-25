cmake_minimum_required(VERSION 3.20)
project(ducklips_tract CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(tract
  src/cord.cpp
  src/tube.cpp
  src/tube_step.cpp
  src/variation.cpp
  src/voice_body.cpp
  src/speaker.cpp
  src/sheet_tube.cpp
  src/somatic_hold.cpp
  src/sheet_speaker.cpp
  src/articulation.cpp
  src/parameter_tree.cpp
  src/tongue.cpp
  src/tongue_roles.cpp
  src/tongue_skeleton.cpp
  src/pinch_tube.cpp
  ../voice/src/geometry.cpp
  ../voice/src/score.cpp
  ../voice/src/block_file.cpp)
# The front end's body language (score) and tract geometry are shared with
# ../voice; its planner writes the scores this engine speaks.
target_include_directories(tract PUBLIC include ../voice/include)
target_compile_options(tract PRIVATE -Wall -Wextra -Wpedantic)
# Tuned for the machine that builds it, multiply-adds fused (measured: see
# NOTES.md, Milestone 2 speed).
option(TRACT_NATIVE "tune for the building CPU" ON)
if(TRACT_NATIVE)
  target_compile_options(tract PRIVATE -mcpu=native -ffp-contract=fast)
endif()

set_target_properties(tract PROPERTIES POSITION_INDEPENDENT_CODE ON)
# The tube's step with its arithmetic pinned (src/tube_step.cpp): no
# multiply-add fused but its explicit std::fma (after the target's
# -ffp-contract=fast on the command line, so this one holds).
set_source_files_properties(src/tube_step.cpp PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")

# The sustain tool's run (tract/sustain_run.hpp), shared by the tool and the
# in-process engine. Compiled with the flags tools/sustain.cpp had (no
# -mcpu=native, no -ffp-contract=fast): its keyframe interpolations must
# round exactly as the reference tool's did.
add_library(tract_run STATIC src/sustain_run.cpp)
target_link_libraries(tract_run PUBLIC tract)
target_compile_options(tract_run PRIVATE -Wall -Wextra -Wpedantic)
set_target_properties(tract_run PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(sustain tools/sustain.cpp)
target_link_libraries(sustain PRIVATE tract_run)

# The in-process engine: a plain C interface (bindings/tract_engine.h) for
# ctypes (analysis/engine.py).
add_library(tract_engine SHARED bindings/tract_engine.cpp)
target_link_libraries(tract_engine PRIVATE tract_run)
target_include_directories(tract_engine PRIVATE bindings)
target_compile_options(tract_engine PRIVATE -Wall -Wextra -Wpedantic)

add_executable(transfer tools/transfer.cpp)
target_link_libraries(transfer PRIVATE tract)

add_executable(speak tools/speak.cpp ../voice/src/room.cpp)
target_link_libraries(speak PRIVATE tract)

add_executable(voices tools/voices.cpp)
target_link_libraries(voices PRIVATE tract)
find_package(Threads REQUIRED)
target_link_libraries(voices PRIVATE Threads::Threads)

add_executable(stream tools/stream.cpp)
target_link_libraries(stream PRIVATE tract)

add_executable(sheet_render tools/sheet_render.cpp)
target_link_libraries(sheet_render PRIVATE tract)

add_executable(sheet_transfer tools/sheet_transfer.cpp)
target_link_libraries(sheet_transfer PRIVATE tract)

add_executable(sheet_areas tools/sheet_areas.cpp)
target_link_libraries(sheet_areas PRIVATE tract)
add_executable(sheet_hold tools/sheet_hold.cpp)
target_link_libraries(sheet_hold PRIVATE tract)

add_executable(sheet_probe tools/sheet_probe.cpp)
target_link_libraries(sheet_probe PRIVATE tract)

add_executable(sheet_speak tools/sheet_speak.cpp ../voice/src/room.cpp)
target_link_libraries(sheet_speak PRIVATE tract)

add_executable(hold_step tools/hold_step.cpp)
target_link_libraries(hold_step PRIVATE tract)

add_executable(sheet_stream tools/sheet_stream.cpp)
target_link_libraries(sheet_stream PRIVATE tract)

add_executable(sheet_voices tools/sheet_voices.cpp)
target_link_libraries(sheet_voices PRIVATE tract)

add_executable(geometry_probe tools/geometry_probe.cpp)
target_link_libraries(geometry_probe PRIVATE tract)

add_executable(tongue_probe tools/tongue_probe.cpp)
target_link_libraries(tongue_probe PRIVATE tract)
