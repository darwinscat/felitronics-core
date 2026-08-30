# felitronics::rigplayer — ONE DEVICE OF A .orbitrig PACK, PLAYING. Which captures sound for a panel
# (namz::rig's policy), the crossfade between them (nam's BlendLaw), the tone knobs as the pack describes
# them (bands → biquads, a curve → a minimum-phase FIR), a blend knob's dry path, the models' alignment.
# Header-only, JUCE-free; process() allocates nothing. Moved here from OrbitCapture NAM
# (app/src/rigplayer, 2026-08-30) so the app and a plugin play through ONE player.
# Deps: nam (NamStage, BlendLaw — and namz's rig headers, which nam exports), convolution (CabConvolver),
# lineareq (MagnitudeCurve), eq (MatchedBiquad). Self-gates on FELITRONICS_WITH_NAM, as nam does.
if(NOT FELITRONICS_WITH_NAM)
    return()
endif()

add_library(felitronics_rigplayer INTERFACE)
add_library(felitronics::rigplayer ALIAS felitronics_rigplayer)
target_include_directories(felitronics_rigplayer INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(felitronics_rigplayer INTERFACE
    felitronics::nam felitronics::convolution felitronics::lineareq felitronics::eq)
target_compile_features(felitronics_rigplayer INTERFACE cxx_std_20)

if(FELITRONICS_BUILD_TESTS)
    foreach(t IN ITEMS RigPlayer ModelBlend ModelAlignment)
        string(TOLOWER "${t}" tl)
        add_executable(felitronics_rigplayer_${tl}_tests tests/${t}Tests.cpp)
        target_link_libraries(felitronics_rigplayer_${tl}_tests PRIVATE felitronics::rigplayer)
        target_include_directories(felitronics_rigplayer_${tl}_tests PRIVATE ${felitronics_core_SOURCE_DIR}/test_support)
        target_compile_features(felitronics_rigplayer_${tl}_tests PRIVATE cxx_std_20)
        add_test(NAME felitronics_rigplayer_${tl}_tests COMMAND felitronics_rigplayer_${tl}_tests)
    endforeach()
endif()
