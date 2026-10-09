cmake_minimum_required(VERSION 3.24)
project(DuckLips LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# No fast-math: the modal propagators, Bernoulli roots and inertia counts rely
# on IEEE behaviour.
add_library(ducklips STATIC
    src/physics.cpp
    src/tree_operator.cpp
    src/acoustic_network.cpp
    src/modal.cpp
    src/moving_mode.cpp
    src/sources.cpp
    src/air_path.cpp
    src/voice.cpp
    src/body_config.cpp
    src/articulation.cpp
    src/tongue.cpp
    src/phonetics.cpp
    src/vowel_library.cpp
    src/cues.cpp
    src/nudges.cpp
    src/synthesizer.cpp
    src/planner.cpp
    src/connected_speech.cpp
    src/wav.cpp
)
target_include_directories(ducklips PUBLIC include)
target_compile_options(ducklips PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)

add_executable(ducklips_say tools/ducklips_say.cpp)
target_link_libraries(ducklips_say PRIVATE ducklips)

add_executable(ducklips_analyze tools/ducklips_analyze.cpp)
target_link_libraries(ducklips_analyze PRIVATE ducklips)

add_executable(ducklips_tests tests/ducklips_tests.cpp)
target_link_libraries(ducklips_tests PRIVATE ducklips)

enable_testing()
add_test(NAME ducklips_tests COMMAND ducklips_tests)

# Optional: DuckLips as a GUI.Forms AudioGenerator. Builds the GUI.Forms
# audio service from the pinned gui_forms submodule source (its Audio.cmake
# fetches and hash-checks the pinned miniaudio revision).
option(DUCKLIPS_WITH_GUI_FORMS_AUDIO "Build the GUI.Forms audio adapter and its offline check" OFF)
if(DUCKLIPS_WITH_GUI_FORMS_AUDIO)
    set(GUI_FORMS_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../../gui_forms" CACHE PATH "GUI.Forms source tree")
    include("${GUI_FORMS_SOURCE}/cmake/Audio.cmake")
    add_library(ducklips_gui_forms_audio STATIC src/gui_forms_voice.cpp)
    target_link_libraries(ducklips_gui_forms_audio PUBLIC ducklips GUIForms::Audio)
    add_executable(ducklips_gui_forms tools/ducklips_gui_forms.cpp)
    target_link_libraries(ducklips_gui_forms PRIVATE ducklips_gui_forms_audio)
    add_test(NAME ducklips_gui_forms_offline COMMAND ducklips_gui_forms --out ${CMAKE_CURRENT_BINARY_DIR}/gui_forms_offline.wav)
endif()

add_executable(ducklips_stream tools/ducklips_stream.cpp)
target_link_libraries(ducklips_stream PRIVATE ducklips)

add_executable(ducklips_stream_tests tests/stream_tests.cpp)
target_link_libraries(ducklips_stream_tests PRIVATE ducklips)
add_test(NAME ducklips_stream_tests COMMAND ducklips_stream_tests)

target_compile_definitions(ducklips_stream_tests PRIVATE DUCKLIPS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(ducklips_moving_mode_tests tests/moving_mode_tests.cpp)
target_link_libraries(ducklips_moving_mode_tests PRIVATE ducklips)
add_test(NAME ducklips_moving_mode_tests COMMAND ducklips_moving_mode_tests)
