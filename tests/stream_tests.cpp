#include "ducklips/body_config.hpp"
#include "ducklips/nudges.hpp"
#include "ducklips/planner.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <sstream>
#include <vector>

namespace {
std::atomic<std::size_t> allocations{0};
std::size_t failures = 0;
void check(bool okay, const char* label) {
    std::printf("%s %s\n", okay ? "PASS" : "FAIL", label);
    if (!okay) { ++failures; }
}
}
void* operator new(std::size_t count) {
    ++allocations;
    void* result = std::malloc(std::max<std::size_t>(count, 1));
    if (!result) { throw std::bad_alloc(); }
    return result;
}
void* operator new[](std::size_t count) { return ::operator new(count); }
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }

int main() {
    ducklips::VoiceProfile body = ducklips::baritone_voice();
    std::string error{};
    std::istringstream good("preset=alto\nname=example\ntract_length=.15\nfolds.lower_mass=.00007\nspeaking_rate=1.2\n");
    check(ducklips::read_body_config(good, body, error) && body.name == "example" && body.tract_length == .15
          && body.folds.lower_mass == .00007 && body.base_pitch == ducklips::alto_voice().base_pitch,
          "body preset and anatomy overrides are applied");
    const char* bad_cases[] = {"tract_length=.18\nunknown=1", "base_pitch=nan", "lung_pressure=inf",
        "tract_length=-1", "base_pitch=100garbage", "base_pitch=110\nbase_pitch=120", "name=test\npreset=alto",
        "audible_breath=yes", "adducted_area=.00008\nabducted_area=.00001"};
    for (const char* text : bad_cases) {
        std::istringstream bad(text);
        check(!ducklips::read_body_config(bad, body, error) && body.name == "example" && body.tract_length == .15,
              "invalid body input fails transactionally");
    }
    ducklips::Lexicon lexicon{};
    std::size_t installed = lexicon.load_cmu_file(std::string(DUCKLIPS_SOURCE_DIR) + "/samples/nato.dict");
    ducklips::WordPronunciation papa = lexicon.pronounce("papa");
    ducklips::WordPronunciation quebec = lexicon.pronounce("quebec");
    check(installed == 30 && papa.phones[1].phone == ducklips::Phone::AA && papa.phones[3].stress == 1
          && quebec.phones[1].phone == ducklips::Phone::EH, "explicit NATO dictionary overrides built-in pronunciations");
    body = ducklips::baritone_voice();
    ducklips::NudgeStream whole{};
    ducklips::NudgeStream chunks{};
    whole.prepare(body);
    chunks.prepare(body);
    const char* lines[] = {"0 lungs 600 40", "0 glottis .000003 170", "0 vowel AA 75",
        ".1 diphthong AY .4 75", ".2 pitch 125 90", ".65 constrict S 1 120", ".65 glottis .000025 170",
        ".8 release tip 150", ".8 lungs 0 40"};
    for (const char* line : lines) {
        check(whole.submit(line, error) && chunks.submit(line, error), "valid overlapping nudge accepted");
    }
    std::size_t pending = whole.pending();
    const char* invalid[] = {"-1 pitch 100 40", "0 vowel S 40", "0 diphthong AA .4 60", "0 pitch nan 50",
        "0 release forehead 40", "0 lungs 600 0", "0 lungs 600 40 garbage", "0 constrict S inf 120"};
    for (const char* line : invalid) {
        check(!whole.submit(line, error) && whole.pending() == pending, "invalid nudge does not change schedule");
    }
    std::vector<float> a(48000 * 2);
    std::vector<float> b(a.size());
    std::size_t before = allocations.load();
    whole.render(a);
    check(allocations.load() == before, "render allocates no heap storage including cue delivery");
    for (std::size_t frame = 0; frame < 48000;) {
        std::size_t run = std::min<std::size_t>(17, 48000 - frame);
        chunks.render(std::span<float>(b.data() + frame * 2, run * 2));
        frame += run;
    }
    check(a == b, "audio is bit-identical for full buffer and 17-frame host chunks");
    double energy = 0;
    bool finite = true;
    for (float value : a) { finite = finite && std::isfinite(value); energy += value * value; }
    check(finite && energy > .001 && !whole.engine().statistics().nonfinite_samples
          && !whole.engine().statistics().clipped_samples, "stream produces finite non-silent unclipped speech audio");
    check(!whole.submit(".9 pitch 100 50", error), "late nudge is rejected");
    check(whole.submit("1 pitch 100 50", error), "next uncommitted block accepts nudge");
    whole.render(std::span<float>(a.data(), 2));
    check(!whole.submit("1.00001 pitch 100 50", error), "partially rendered control block rejects retroactive nudge");
    ducklips::SynthesisOptions canonical_options{};canonical_options.canonical_motion=true;
    ducklips::NudgeStream canonical_whole{},canonical_chunks{};
    canonical_options.diagnostics=true;
    canonical_whole.prepare(body,1,4096,canonical_options);
    canonical_options.diagnostics=false;
    canonical_chunks.prepare(body,1,4096,canonical_options);
    for (const char* line:lines) {
        check(canonical_whole.submit(line,error) && canonical_chunks.submit(line,error),"canonical nudge accepted");
    }
    before=allocations.load();canonical_whole.render(a);
    check(allocations.load()==before,"canonical render retains state without heap allocation");
    for (std::size_t frame=0;frame<48000;) {
        const std::size_t run=std::min<std::size_t>(17,48000-frame);
        canonical_chunks.render(std::span<float>(b.data()+2*frame,2*run));frame+=run;
    }
    energy=0;for (float value:a) { energy+=value*value; }
    check(a==b,"canonical output is independent of host buffer size and geometry diagnostics");
    check(energy>1e-6 && !canonical_whole.engine().statistics().nonfinite_samples
          && !canonical_whole.engine().statistics().clipped_samples,"canonical nudges generate finite nonzero unclipped audio");
    check(canonical_whole.engine().statistics().maximum_cavity_mass_residual<1e-12,
          "canonical moving-cavity mass balance holds during streamed diphthongs and frication");
    check(canonical_whole.engine().statistics().factorizations==0,"canonical speech never enters the old modal solver");
    check(canonical_whole.submit("1 pitch 140 70",error),"new nudge accepted after canonical rendering begins");
    canonical_whole.render(std::span<float>(a.data(),960));
    for(bool spatial:{false,true}){
    canonical_options.spatial_pressure=spatial;
    canonical_options.conservative_tongue=true;canonical_options.diagnostics=true;
    canonical_whole.prepare(body,1,4096,canonical_options);
    canonical_options.diagnostics=false;canonical_chunks.prepare(body,1,4096,canonical_options);
    const char* material_lines[]={"0 lungs 600 40","0 glottis .000003 170","0 vowel AA 75",
        ".03 constrict K 1 120",".14 release body 120",".17 constrict S 1 120",".25 release tip 150"};
    for (const char* line:material_lines) {
        check(canonical_whole.submit(line,error) && canonical_chunks.submit(line,error),"conservative material nudge accepted");
    }
    constexpr std::size_t material_frames=14400;
    before=allocations.load();canonical_whole.render(std::span<float>(a.data(),2*material_frames));
    check(allocations.load()==before,spatial?"spatial airflow and conservative tongue allocate no render storage":"conservative tongue and diagnostics allocate no render storage");
    for (std::size_t frame=0;frame<material_frames;) {
        const std::size_t run=std::min<std::size_t>(17,material_frames-frame);
        canonical_chunks.render(std::span<float>(b.data()+2*frame,2*run));frame+=run;
    }
    check(std::equal(a.begin(),a.begin()+2*material_frames,b.begin()),spatial?"spatial airflow output is independent of host chunks and diagnostics":"conservative tongue output is independent of host chunks and diagnostics");
    const ducklips::SynthesisStatistics& material_stats=canonical_whole.engine().statistics();
    check(!material_stats.airflow_failures && !material_stats.geometry_failures && !material_stats.nonfinite_samples && !material_stats.clipped_samples,
          "combined body and tip nudges pass geometry and output audits");
    check(material_stats.maximum_cavity_mass_residual<(spatial?1e-10:1e-12),"conservative tongue retains moving cavity balance");
    }
    canonical_options.spatial_pressure=false;
    ducklips::Synthesizer invalid_geometry{};invalid_geometry.prepare(body,canonical_options);
    ducklips::Cue unsupported{};unsupported.kind=ducklips::CueKind::constriction;unsupported.slot=ducklips::ConstrictionSlot::tongue_body;
    unsupported.gesture.place_from_lips=.10;unsupported.gesture.width=.009;unsupported.value=1;unsupported.angular_rate=120;
    check(invalid_geometry.enqueue(unsupported),"unsupported geometry cue queued for explicit render audit");
    invalid_geometry.render(std::span<float>(a.data(),96));
    check(invalid_geometry.statistics().geometry_failures==1 && invalid_geometry.statistics().geometry_error
          && invalid_geometry.finished(),"unsupported geometry stops rendering and records the error");
    bool failed_silence=true;for (std::size_t i=0;i<96;++i) { failed_silence=failed_silence&&a[i]==0; }
    check(failed_silence,"failed geometry emits no substituted rigid-body audio");
    ducklips::NudgeStream bounded{};
    bounded.prepare(body, 2.0, 4);
    check(bounded.submit("0 diphthong AY .4 75", error) && bounded.submit("0 diphthong AW .4 75", error),
          "bounded stream accepts exactly its cue capacity");
    check(!bounded.submit("1 pitch 110 90", error) && bounded.pending() == 4, "bounded stream rejects overflow");
    std::vector<float> half_second(24000 * 2);
    bounded.render(half_second);
    check(bounded.pending() == 0, "speed-scaled diphthong targets are delivered");
    // Regression for the audible control-clock comb observed in Alfa's inhale.
    // This measures the rendered signal, independently of the control code.
    ducklips::Synthesizer speech{};
    ducklips::SynthesisOptions speech_options{};
    speech.prepare(ducklips::baritone_voice(), speech_options);
    ducklips::UtterancePlanner planner{};
    planner.prepare(ducklips::baritone_voice(), speech.vowels(), speech_options.rate);
    const ducklips::UtterancePlan plan = planner.plan_phonemes("AE1 L F AX0");
    for (const ducklips::Cue& cue : plan.cues) { check(speech.enqueue(cue), "clock-regression cue accepted"); }
    std::vector<float> spoken(62880*2);
    speech.render(spoken);
    double largest_clock_line=0.0;
    for (double frequency : {1000.0, 2000.0}) {
        double real=0.0,imaginary=0.0,weight_sum=0.0;
        for (std::size_t sample=0;sample<5760;++sample) {
            const double window=0.5-0.5*std::cos(2*ducklips::pi*static_cast<double>(sample)/5759.0);
            const double phase=2*ducklips::pi*frequency*static_cast<double>(sample)/48000.0;
            const double value=spoken[2*(7680+sample)]*window;
            real+=value*std::cos(phase);imaginary+=value*std::sin(phase);weight_sum+=window;
        }
        largest_clock_line=std::max(largest_clock_line,2*std::hypot(real,imaginary)/weight_sum);
    }
    std::printf("Quiet Alfa control-clock amplitude %.9g (limit 1e-5 full scale)\n",largest_clock_line);
    check(largest_clock_line<1e-5, "physiological trajectory does not impose the old 1/2 kHz clock comb");
    double voiced_energy=0.0;
    for (std::size_t frame=28800;frame<48000;++frame) { voiced_energy+=spoken[2*frame]*spoken[2*frame]; }
    check(voiced_energy/19200>1e-6, "clock-line regression retains audible speech rather than muting output");
    std::printf("%zu failures\n", failures);
    return failures ? 1 : 0;
}
