#include "ducklips/gui_forms_voice.hpp"

namespace ducklips {

void GuiFormsVoice::render(std::span<float> stereo) noexcept {
    synthesizer_.render(stereo);
}

Synthesizer& GuiFormsVoice::synthesizer() noexcept {
    return synthesizer_;
}

}
