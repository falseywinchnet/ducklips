#pragma once

namespace ducklips {

// The natural speed of each part of the body (rad/s pole rate of its
// motion). The tongue tip is the quickest articulator; the tongue body
// carries more mass; the jaw is slower still; breathing is slowest.
struct ArticulatorSpeeds final {
    double tongue_tip{230.0};
    double lips{180.0};
    double tongue_body{75.0};
    double jaw{45.0};
    double velum{65.0};
    double glottis{170.0};
    double lungs{40.0};
    double larynx_pitch{90.0};
};

}
