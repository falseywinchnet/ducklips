#include "ducklips/cues.hpp"

namespace ducklips {

bool CueEarlier::operator()(const Cue& left, const Cue& right) const noexcept {
    bool earlier = left.time < right.time;
    return earlier;
}

void CueQueue::prepare(std::size_t capacity) {
    slots_.assign(capacity + 1, Cue{});
    write_.store(0, std::memory_order_relaxed);
    read_.store(0, std::memory_order_relaxed);
}

bool CueQueue::push(const Cue& cue) noexcept {
    std::size_t write = write_.load(std::memory_order_relaxed);
    std::size_t next = write + 1;
    if (next == slots_.size()) {
        next = 0;
    }
    if (next == read_.load(std::memory_order_acquire)) {
        return false;
    }
    slots_[write] = cue;
    write_.store(next, std::memory_order_release);
    return true;
}

bool CueQueue::peek(Cue& output) const noexcept {
    std::size_t read = read_.load(std::memory_order_relaxed);
    if (read == write_.load(std::memory_order_acquire)) {
        return false;
    }
    output = slots_[read];
    return true;
}

void CueQueue::pop() noexcept {
    std::size_t read = read_.load(std::memory_order_relaxed);
    std::size_t next = read + 1;
    if (next == slots_.size()) {
        next = 0;
    }
    read_.store(next, std::memory_order_release);
}

std::size_t CueQueue::free_space() const noexcept {
    std::size_t write = write_.load(std::memory_order_acquire);
    std::size_t read = read_.load(std::memory_order_acquire);
    std::size_t used = write >= read ? write - read : slots_.size() - read + write;
    std::size_t space = slots_.size() - 1 - used;
    return space;
}

}
