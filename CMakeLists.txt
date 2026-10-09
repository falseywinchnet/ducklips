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
    src/canonical_speech.cpp
    src/moving_tree.cpp
    src/moving_graph.cpp
    src/graph_orifice.cpp
    src/graph_phonation.cpp
    src/gas_storage.cpp
    src/glottal_storage.cpp
    src/graph_glottal_storage.cpp
    src/continuous_modes.cpp
    src/modal_transport.cpp
    src/motion_segment.cpp
    src/sources.cpp
    src/cavity_network.cpp
    src/mean_airway.cpp
    src/planar_airspace.cpp
    src/planar_face_motion.cpp
    src/air_path.cpp
    src/voice.cpp
    src/body_config.cpp
    src/articulation.cpp
    src/lip_profile.cpp
    src/tongue.cpp
    src/tongue_flow.cpp
    src/tongue_surface.cpp
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

add_executable(ducklips_moving_tree_probe tools/moving_tree_probe.cpp)
target_link_libraries(ducklips_moving_tree_probe PRIVATE ducklips)
add_test(NAME ducklips_moving_tree_reference COMMAND ducklips_moving_tree_probe)

add_executable(ducklips_continuous_modes_probe tools/continuous_modes_probe.cpp)
target_link_libraries(ducklips_continuous_modes_probe PRIVATE ducklips)
add_test(NAME ducklips_continuous_modes_reference COMMAND ducklips_continuous_modes_probe)

add_executable(ducklips_modal_transport_probe tools/modal_transport_probe.cpp)
target_link_libraries(ducklips_modal_transport_probe PRIVATE ducklips)
add_test(NAME ducklips_modal_transport_reference COMMAND ducklips_modal_transport_probe)

add_executable(ducklips_motion_segment_probe tools/motion_segment_probe.cpp)
target_link_libraries(ducklips_motion_segment_probe PRIVATE ducklips)
add_test(NAME ducklips_motion_segment_reference COMMAND ducklips_motion_segment_probe)

add_executable(ducklips_tongue_geometry_probe tools/tongue_geometry_probe.cpp)
target_link_libraries(ducklips_tongue_geometry_probe PRIVATE ducklips)

add_executable(ducklips_tongue_flow_probe tools/tongue_flow_probe.cpp)
target_link_libraries(ducklips_tongue_flow_probe PRIVATE ducklips)

add_executable(ducklips_tongue_flow_tests tests/tongue_flow_tests.cpp)
target_link_libraries(ducklips_tongue_flow_tests PRIVATE ducklips)
add_test(NAME ducklips_tongue_flow_tests COMMAND ducklips_tongue_flow_tests)

add_executable(ducklips_tongue_surface_tests tests/tongue_surface_tests.cpp)
target_link_libraries(ducklips_tongue_surface_tests PRIVATE ducklips)
add_test(NAME ducklips_tongue_surface_tests COMMAND ducklips_tongue_surface_tests)

add_executable(ducklips_planner_tests tests/planner_tests.cpp)
target_link_libraries(ducklips_planner_tests PRIVATE ducklips)
target_compile_definitions(ducklips_planner_tests PRIVATE DUCKLIPS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME ducklips_planner_tests COMMAND ducklips_planner_tests)

add_executable(ducklips_tongue_reach_audit tools/tongue_reach_audit.cpp)
target_link_libraries(ducklips_tongue_reach_audit PRIVATE ducklips)

add_executable(ducklips_cavity_network_tests tests/cavity_network_tests.cpp)
target_link_libraries(ducklips_cavity_network_tests PRIVATE ducklips)
add_test(NAME ducklips_cavity_network_tests COMMAND ducklips_cavity_network_tests)

add_executable(ducklips_cavity_network_probe tools/cavity_network_probe.cpp)
target_link_libraries(ducklips_cavity_network_probe PRIVATE ducklips)

add_executable(ducklips_canonical_closure_tests tests/canonical_closure_tests.cpp)
target_link_libraries(ducklips_canonical_closure_tests PRIVATE ducklips)
add_test(NAME ducklips_canonical_closure_tests COMMAND ducklips_canonical_closure_tests)

add_executable(ducklips_mean_airway_tests tests/mean_airway_tests.cpp)
target_link_libraries(ducklips_mean_airway_tests PRIVATE ducklips)
add_test(NAME ducklips_mean_airway_tests COMMAND ducklips_mean_airway_tests)

add_executable(ducklips_planar_airspace_tests tests/planar_airspace_tests.cpp)
target_link_libraries(ducklips_planar_airspace_tests PRIVATE ducklips)
add_test(NAME ducklips_planar_airspace_tests COMMAND ducklips_planar_airspace_tests)

add_executable(ducklips_airspace_geometry_probe tools/airspace_geometry_probe.cpp)
target_link_libraries(ducklips_airspace_geometry_probe PRIVATE ducklips)

add_executable(ducklips_cavity_graph_probe tools/cavity_graph_probe.cpp)
target_link_libraries(ducklips_cavity_graph_probe PRIVATE ducklips)

add_executable(ducklips_cavity_transport_probe tools/cavity_transport_probe.cpp)
target_link_libraries(ducklips_cavity_transport_probe PRIVATE ducklips)

add_executable(ducklips_planar_regions_probe tools/planar_regions_probe.cpp)
target_link_libraries(ducklips_planar_regions_probe PRIVATE ducklips)

add_executable(ducklips_planar_face_motion_tests tests/planar_face_motion_tests.cpp)
target_link_libraries(ducklips_planar_face_motion_tests PRIVATE ducklips)
add_test(NAME ducklips_planar_face_motion_tests COMMAND ducklips_planar_face_motion_tests)

add_executable(ducklips_cavity_region_motion_probe tools/cavity_region_motion_probe.cpp)
target_link_libraries(ducklips_cavity_region_motion_probe PRIVATE ducklips)

add_executable(ducklips_planar_pressure_motion_probe tools/planar_pressure_motion_probe.cpp)
target_link_libraries(ducklips_planar_pressure_motion_probe PRIVATE ducklips)

add_executable(ducklips_planar_inertive_motion_probe tools/planar_inertive_motion_probe.cpp)
target_link_libraries(ducklips_planar_inertive_motion_probe PRIVATE ducklips)

add_executable(ducklips_graph_glottal_probe tools/graph_glottal_probe.cpp)
target_link_libraries(ducklips_graph_glottal_probe PRIVATE ducklips)
add_test(NAME ducklips_graph_glottal_reference COMMAND ducklips_graph_glottal_probe)

add_executable(ducklips_glottal_storage_probe tools/glottal_storage_probe.cpp)
target_link_libraries(ducklips_glottal_storage_probe PRIVATE ducklips)
add_test(NAME ducklips_glottal_storage_reference COMMAND ducklips_glottal_storage_probe)

add_executable(ducklips_graph_phonation_probe tools/graph_phonation_probe.cpp)
target_link_libraries(ducklips_graph_phonation_probe PRIVATE ducklips)
add_test(NAME ducklips_graph_phonation_reference COMMAND ducklips_graph_phonation_probe)

add_executable(ducklips_graph_orifice_probe tools/graph_orifice_probe.cpp)
target_link_libraries(ducklips_graph_orifice_probe PRIVATE ducklips)
add_test(NAME ducklips_graph_orifice_reference COMMAND ducklips_graph_orifice_probe)

add_executable(ducklips_moving_graph_probe tools/moving_graph_probe.cpp)
target_link_libraries(ducklips_moving_graph_probe PRIVATE ducklips)
add_test(NAME ducklips_moving_graph_reference COMMAND ducklips_moving_graph_probe)

add_executable(ducklips_lip_profile_probe tools/lip_profile_probe.cpp)
target_link_libraries(ducklips_lip_profile_probe PRIVATE ducklips)
add_test(NAME ducklips_lip_profile_reference COMMAND ducklips_lip_profile_probe)
