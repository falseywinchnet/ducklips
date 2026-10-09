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
    src/glottal_kinetics.cpp
    src/glottal_contour.cpp
    src/continuous_glottis.cpp
    src/glottal_relaxation.cpp
    src/graph_continuous_glottis.cpp
    src/contour_pressure_loss.cpp
    src/contour_pressure_flow.cpp
    src/contour_pressure_storage.cpp
    src/contour_kinetics.cpp
    src/fold_contact.cpp
    src/contour_glottal_storage.cpp
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

add_executable(ducklips_glottal_contour_probe tools/glottal_contour_probe.cpp)
target_link_libraries(ducklips_glottal_contour_probe PRIVATE ducklips)
add_test(NAME ducklips_glottal_contour_reference COMMAND ducklips_glottal_contour_probe)

add_executable(ducklips_glottal_geometry_motion_probe tools/glottal_geometry_motion_probe.cpp)
target_link_libraries(ducklips_glottal_geometry_motion_probe PRIVATE ducklips)

add_executable(ducklips_glottal_kinetics_probe tools/glottal_kinetics_probe.cpp)
target_link_libraries(ducklips_glottal_kinetics_probe PRIVATE ducklips)
add_test(NAME ducklips_glottal_kinetics_reference COMMAND ducklips_glottal_kinetics_probe)

add_executable(ducklips_graph_glottal_probe tools/graph_glottal_probe.cpp)
target_link_libraries(ducklips_graph_glottal_probe PRIVATE ducklips)
add_test(NAME ducklips_graph_glottal_reference COMMAND ducklips_graph_glottal_probe)
add_test(NAME ducklips_graph_glottal_geometry_reference COMMAND ducklips_graph_glottal_probe --articulated)

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

add_executable(ducklips_contour_kinetics_probe tools/contour_kinetics_probe.cpp)
target_link_libraries(ducklips_contour_kinetics_probe PRIVATE ducklips)
add_test(NAME ducklips_contour_kinetics_reference COMMAND ducklips_contour_kinetics_probe)

add_executable(ducklips_contour_glottal_probe tools/contour_glottal_probe.cpp)
target_link_libraries(ducklips_contour_glottal_probe PRIVATE ducklips)
add_test(NAME ducklips_contour_glottal_reference COMMAND ducklips_contour_glottal_probe)

add_executable(ducklips_contour_source_controls tools/contour_source_controls.cpp)
target_link_libraries(ducklips_contour_source_controls PRIVATE ducklips)
add_test(NAME ducklips_contour_source_controls COMMAND ducklips_contour_source_controls)

add_executable(ducklips_contour_contact_probe tools/contour_contact_probe.cpp)
target_link_libraries(ducklips_contour_contact_probe PRIVATE ducklips)

add_executable(ducklips_contour_sealed_probe tools/contour_sealed_probe.cpp)
target_link_libraries(ducklips_contour_sealed_probe PRIVATE ducklips)
add_test(NAME ducklips_contour_sealed_reference COMMAND ducklips_contour_sealed_probe)

add_executable(ducklips_contour_cycle_probe tools/contour_cycle_probe.cpp)
target_link_libraries(ducklips_contour_cycle_probe PRIVATE ducklips)
add_test(NAME ducklips_contour_lower_cycle COMMAND ducklips_contour_cycle_probe 2000 lower)
add_test(NAME ducklips_contour_upper_cycle COMMAND ducklips_contour_cycle_probe 2000 upper)
add_test(NAME ducklips_contour_dissipative_reference COMMAND ducklips_contour_glottal_probe --dissipative)

add_executable(ducklips_contour_dissipation_probe tools/contour_dissipation_probe.cpp)
target_link_libraries(ducklips_contour_dissipation_probe PRIVATE ducklips)
add_test(NAME ducklips_contour_dissipation_kernel COMMAND ducklips_contour_dissipation_probe)

add_test(NAME ducklips_contour_relaxation_reference COMMAND ducklips_contour_glottal_probe --relaxation)
add_test(NAME ducklips_contour_relaxation_lower_cycle COMMAND ducklips_contour_cycle_probe 2000 lower --relaxation)
add_test(NAME ducklips_contour_relaxation_upper_cycle COMMAND ducklips_contour_cycle_probe 2000 upper --relaxation)

add_executable(ducklips_contour_pocket_probe tools/contour_pocket_probe.cpp)
target_link_libraries(ducklips_contour_pocket_probe PRIVATE ducklips)

add_executable(ducklips_contour_pressure_probe tools/contour_pressure_probe.cpp)
target_link_libraries(ducklips_contour_pressure_probe PRIVATE ducklips)
add_test(NAME ducklips_contour_pressure_storage COMMAND ducklips_contour_pressure_probe)

add_executable(ducklips_contour_pressure_flow_probe tools/contour_pressure_flow_probe.cpp)
target_link_libraries(ducklips_contour_pressure_flow_probe PRIVATE ducklips)
add_test(NAME ducklips_contour_pressure_flow COMMAND ducklips_contour_pressure_flow_probe)

add_executable(ducklips_contour_pressure_loss_probe tools/contour_pressure_loss_probe.cpp)
target_link_libraries(ducklips_contour_pressure_loss_probe PRIVATE ducklips)
add_test(NAME ducklips_contour_pressure_loss COMMAND ducklips_contour_pressure_loss_probe)

add_executable(ducklips_continuous_glottis_probe tools/continuous_glottis_probe.cpp)
target_link_libraries(ducklips_continuous_glottis_probe PRIVATE ducklips)
add_test(NAME ducklips_continuous_glottis COMMAND ducklips_continuous_glottis_probe)

add_executable(ducklips_continuous_glottis_tangent_probe tools/continuous_glottis_tangent_probe.cpp)
target_link_libraries(ducklips_continuous_glottis_tangent_probe PRIVATE ducklips)

add_executable(ducklips_contour_profile_limit_probe tools/contour_profile_limit_probe.cpp)
target_link_libraries(ducklips_contour_profile_limit_probe PRIVATE ducklips)
add_executable(ducklips_continuous_contact_probe tools/continuous_contact_probe.cpp)
target_link_libraries(ducklips_continuous_contact_probe PRIVATE ducklips)
add_executable(ducklips_continuous_contact_controls tools/continuous_contact_controls.cpp)
target_link_libraries(ducklips_continuous_contact_controls PRIVATE ducklips)
add_test(NAME ducklips_continuous_contact_lower COMMAND ducklips_continuous_contact_controls lower)
add_test(NAME ducklips_continuous_contact_upper COMMAND ducklips_continuous_contact_controls upper)

add_executable(ducklips_continuous_phonation_probe tools/continuous_phonation_probe.cpp)
target_link_libraries(ducklips_continuous_phonation_probe PRIVATE ducklips)
add_executable(ducklips_continuous_phonation_step_probe tools/continuous_phonation_step_probe.cpp)
target_link_libraries(ducklips_continuous_phonation_step_probe PRIVATE ducklips)
add_executable(ducklips_continuous_phonation_controls tools/continuous_phonation_controls.cpp)
target_link_libraries(ducklips_continuous_phonation_controls PRIVATE ducklips)
add_test(NAME ducklips_continuous_phonation_700 COMMAND ducklips_continuous_phonation_controls 700)
add_test(NAME ducklips_continuous_phonation_1400 COMMAND ducklips_continuous_phonation_controls 1400)
add_executable(ducklips_continuous_glottis_rebase_controls tools/continuous_glottis_rebase_controls.cpp)
target_link_libraries(ducklips_continuous_glottis_rebase_controls PRIVATE ducklips)
add_test(NAME ducklips_continuous_glottis_rebase COMMAND ducklips_continuous_glottis_rebase_controls)
add_executable(ducklips_graph_continuous_glottis_probe tools/graph_continuous_glottis_probe.cpp)
target_link_libraries(ducklips_graph_continuous_glottis_probe PRIVATE ducklips)
add_test(NAME ducklips_graph_continuous_glottis COMMAND ducklips_graph_continuous_glottis_probe)
add_executable(ducklips_continuous_tract_probe tools/continuous_tract_probe.cpp)
target_link_libraries(ducklips_continuous_tract_probe PRIVATE ducklips)
add_executable(ducklips_glottal_relaxation_probe tools/glottal_relaxation_probe.cpp)
target_link_libraries(ducklips_glottal_relaxation_probe PRIVATE ducklips)
add_test(NAME ducklips_glottal_relaxation COMMAND ducklips_glottal_relaxation_probe)
add_test(NAME ducklips_graph_continuous_relaxation COMMAND ducklips_graph_continuous_glottis_probe --relaxation)
