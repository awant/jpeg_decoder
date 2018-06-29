#include "reader.h"
#include <iostream>

ByteStreamReader::ByteStreamReader(std::istream& stream): stream_(stream) {}

uint8_t ByteStreamReader::ReadBit() {
    if (cache_size_ == MAX_CACHE_SIZE) { // cache is empty
        cache_ = ReadRawByte();
        cache_size_ = 0;
    }
    uint8_t result = cache_ >> 7;
    cache_ = (cache_ << 1) & 0xff;
    ++cache_size_;
    return result;
}

uint8_t ByteStreamReader::ReadHalfByte() {
    if (cache_size_ <= 4) { // then, we have half of byte in cache
        uint8_t result = cache_ >> 4;
        cache_ = (cache_ << 4) & 0xff;
        cache_size_ += 4;
        return result;
    }
    uint8_t result = 0;
    for (int i = 0; i < 4; ++i) {
        result = (result << 1) | ReadBit();
    }
    return result;
}

uint8_t ByteStreamReader::ReadByte() {
    if (cache_size_ == 0) {
        uint8_t result = cache_;
        cache_ = 0;
        cache_size_ = MAX_CACHE_SIZE;
        return result;
    }
    uint8_t result = 0;
    for (int i = 0; i < MAX_CACHE_SIZE; ++i) {
        result = (result << 1) | ReadBit();
    }
    return result;
}

uint16_t ByteStreamReader::ReadWord() {
    uint16_t first_byte = ReadByte();
    uint8_t second_byte = ReadByte();
    uint16_t result = (first_byte << 8) + second_byte;
    return result;
}

void ByteStreamReader::Read(char* buffer, size_t size) {
    stream_.read(buffer, size);
    if (stream_.eof() || (stream_.gcount() != size)) {
        throw std::runtime_error("Can't read all buffer");
    }
}

bool ByteStreamReader::IsEnded() const {
    return stream_.eof() and (cache_size_ == MAX_CACHE_SIZE);
}

bool ByteStreamReader::IsCacheEmpty() const {
    return cache_size_ == MAX_CACHE_SIZE;
}

uint8_t ByteStreamReader::ReadRawByte() {
    uint8_t byte;
    stream_.read((char*)&byte, 1);
    std::cout << std::hex << int(byte) << "\n";
    if (!stream_) {
        std::cout << "END\n";
    }
    return byte;
}
