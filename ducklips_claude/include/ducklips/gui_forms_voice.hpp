#pragma once
// GUI.Forms Audio adapter. Compiled only when GUI.Forms Audio is available
// (CMake option DUCKLIPS_WITH_GUI_FORMS_AUDIO).
#include "ducklips/synthesizer.hpp"
#include "gui_forms/audio/audio.hpp"
#include <span>

namespace ducklips {

// A DuckLips body as a live GUI.Forms audio source. The mixer pulls
// interleaved L,R 48 kHz frames from render() on the device callback (or
// inside offline AudioEngine::render); Synthesizer::render satisfies that
// contract: it fills the span, never blocks, allocates, locks or throws.
// Cues are pushed from the control thread with synthesizer().enqueue().
// Ownership: the voice owns its synthesizer; the engine shares ownership of
// the voice through the std::shared_ptr passed to AudioEngine::generator.
class GuiFormsVoice final : public gui_forms::AudioGenerator {
public:
    GuiFormsVoice() = default;
    GuiFormsVoice(const GuiFormsVoice&) = delete;
    GuiFormsVoice& operator=(const GuiFormsVoice&) = delete;

    void render(std::span<float> stereo) noexcept override;
    [[nodiscard]] Synthesizer& synthesizer() noexcept;

private:
    Synthesizer synthesizer_{};
};

}
