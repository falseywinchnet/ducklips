#include "ducklips/wav.hpp"
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <cstring>
#include <fstream>

namespace ducklips {
namespace {

void put_u32(std::vector<unsigned char>& bytes, std::uint32_t value) {
    for (std::size_t shift = 0; shift < 32; shift += 8) {
        bytes.push_back(static_cast<unsigned char>((value >> shift) & 0xFFu));
    }
}

void put_u16(std::vector<unsigned char>& bytes, std::uint16_t value) {
    bytes.push_back(static_cast<unsigned char>(value & 0xFFu));
    bytes.push_back(static_cast<unsigned char>((value >> 8) & 0xFFu));
}

void put_text(std::vector<unsigned char>& bytes, const char* text) {
    for (std::size_t index = 0; index < 4; ++index) {
        bytes.push_back(static_cast<unsigned char>(text[index]));
    }
}

std::uint32_t get_u32(const std::vector<unsigned char>& bytes, std::size_t offset) noexcept {
    std::uint32_t value = 0;
    for (std::size_t index = 0; index < 4; ++index) {
        value |= static_cast<std::uint32_t>(bytes[offset + index]) << (8 * index);
    }
    return value;
}

std::uint16_t get_u16(const std::vector<unsigned char>& bytes, std::size_t offset) noexcept {
    std::uint16_t value = static_cast<std::uint16_t>(bytes[offset] | (bytes[offset + 1] << 8));
    return value;
}

}

bool write_wav_float(const std::string& path, const AudioBuffer& audio) {
    std::vector<unsigned char> bytes{};
    std::uint32_t data_bytes = static_cast<std::uint32_t>(audio.samples.size() * sizeof(float));
    std::uint16_t channels = static_cast<std::uint16_t>(audio.channels);
    std::uint32_t rate = static_cast<std::uint32_t>(audio.sample_rate);
    put_text(bytes, "RIFF");
    put_u32(bytes, 36 + data_bytes);
    put_text(bytes, "WAVE");
    put_text(bytes, "fmt ");
    put_u32(bytes, 16);
    put_u16(bytes, 3); // IEEE float
    put_u16(bytes, channels);
    put_u32(bytes, rate);
    put_u32(bytes, rate * channels * 4);
    put_u16(bytes, static_cast<std::uint16_t>(channels * 4));
    put_u16(bytes, 32);
    put_text(bytes, "data");
    put_u32(bytes, data_bytes);
    std::size_t header = bytes.size();
    bytes.resize(header + data_bytes);
    std::memcpy(bytes.data() + header, audio.samples.data(), data_bytes);
    std::ofstream file(path, std::ios::binary);
    if (!file) {
        return false;
    }
    file.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    bool written = static_cast<bool>(file);
    return written;
}

bool read_wav(const std::string& path, AudioBuffer& output) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        return false;
    }
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    if (bytes.size() < 12 || std::memcmp(bytes.data(), "RIFF", 4) != 0 || std::memcmp(bytes.data() + 8, "WAVE", 4) != 0) {
        return false;
    }
    std::uint16_t format = 0;
    std::uint16_t channels = 0;
    std::uint32_t rate = 0;
    std::uint16_t bits = 0;
    std::size_t data_offset = 0;
    std::size_t data_size = 0;
    std::size_t offset = 12;
    while (offset + 8 <= bytes.size()) {
        std::uint32_t chunk_size = get_u32(bytes, offset + 4);
        if (std::memcmp(bytes.data() + offset, "fmt ", 4) == 0 && offset + 24 <= bytes.size()) {
            format = get_u16(bytes, offset + 8);
            channels = get_u16(bytes, offset + 10);
            rate = get_u32(bytes, offset + 12);
            bits = get_u16(bytes, offset + 22);
            if (format == 0xFFFE && offset + 34 <= bytes.size()) {
                format = get_u16(bytes, offset + 32); // WAVE_FORMAT_EXTENSIBLE subformat
            }
        } else if (std::memcmp(bytes.data() + offset, "data", 4) == 0) {
            data_offset = offset + 8;
            data_size = std::min<std::size_t>(chunk_size, bytes.size() - data_offset);
        }
        offset += 8 + chunk_size + (chunk_size & 1u);
    }
    if (channels == 0 || data_offset == 0) {
        return false;
    }
    AudioBuffer audio{};
    audio.channels = channels;
    audio.sample_rate = static_cast<double>(rate);
    if (format == 3 && bits == 32) {
        std::size_t count = data_size / 4;
        audio.samples.resize(count);
        std::memcpy(audio.samples.data(), bytes.data() + data_offset, count * 4);
    } else if (format == 1 && bits == 16) {
        std::size_t count = data_size / 2;
        audio.samples.resize(count);
        for (std::size_t index = 0; index < count; ++index) {
            std::int16_t value = static_cast<std::int16_t>(get_u16(bytes, data_offset + 2 * index));
            audio.samples[index] = static_cast<float>(value) / 32768.0f;
        }
    } else if (format == 1 && bits == 24) {
        std::size_t count = data_size / 3;
        audio.samples.resize(count);
        for (std::size_t index = 0; index < count; ++index) {
            std::size_t base = data_offset + 3 * index;
            std::int32_t value = static_cast<std::int32_t>(bytes[base] | (bytes[base + 1] << 8) | (bytes[base + 2] << 16));
            if ((value & 0x800000) != 0) {
                value -= 0x1000000;
            }
            audio.samples[index] = static_cast<float>(value) / 8388608.0f;
        }
    } else {
        return false;
    }
    output = audio;
    return true;
}

}
