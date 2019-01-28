#include "reader.h"
#include <iostream>

ByteStreamReader::ByteStreamReader(std::istream& stream): stream_(stream) {}

uint8_t ByteStreamReader::ReadRawByte() {
    auto byte = static_cast<uint8_t>(stream_.get());
    if (!stream_ || !stream_.good()) {
        throw std::runtime_error("Can't read stream");
    }
    last_read_word_ = (last_read_word_ << 8) | byte;
//    if (last_read_word_ == 0xffd9) {
//        throw std::runtime_error("Read end marker");
//    }
    return byte;
}

uint8_t ByteStreamReader::ReadBit() {
    if (cache_size_ == MAX_CACHE_SIZE) { // cache is empty
        cache_ = ReadRawByte();
        if (cache_ == 0xff) {
            if (ReadRawByte() != 0) {
                throw std::runtime_error("Strange marker inside section");
            }
        }
        cache_size_ = 0;
    }
    uint8_t result = cache_ >> 7;
    cache_ <<= 1;
    ++cache_size_;
    offset_ += 1./8;
    return result;
}

uint8_t ByteStreamReader::ReadByte() {
    uint8_t result = cache_ >> cache_size_;
    cache_ = ReadRawByte();
    result = (result << cache_size_) | (cache_ >> (MAX_CACHE_SIZE - cache_size_));
    cache_ <<= cache_size_;
    ++offset_;
    return result;
}

uint8_t ByteStreamReader::ReadHalfByte() {
    offset_ += 0.5;
    if (cache_size_ <= 4) { // when we have half of byte in cache, just return it
        uint8_t result = cache_ >> 4;
        cache_ <<= 4;
        cache_size_ += 4;
        assert(cache_size_ <= MAX_CACHE_SIZE);
        return result;
    }
    uint8_t result = cache_ >> cache_size_;
    cache_ = ReadRawByte();
    result = (result << (cache_size_ - 4)) | (cache_ >> (MAX_CACHE_SIZE - (cache_size_ - 4)));
    cache_ <<= cache_size_ - 4;
    cache_size_ = MAX_CACHE_SIZE - (cache_size_ - 4);
    return result;
}

uint16_t ByteStreamReader::ReadWord() {
    uint16_t first_byte = ReadByte();
    uint8_t second_byte = ReadByte();
    uint16_t result = (first_byte << 8) | second_byte;
    return result;
}

void ByteStreamReader::Read(char* buffer, size_t size) {
    stream_.read(buffer, size);
    offset_ += stream_.gcount();
    if (stream_.eof() || (stream_.gcount() != size)) {
        throw std::runtime_error("Can't read all buffer");
    }
}

void ByteStreamReader::CleanCache() {
    cache_size_ = MAX_CACHE_SIZE;
    cache_ = 0;
}

bool ByteStreamReader::IsEnded() const {
    return (cache_size_ == MAX_CACHE_SIZE) && !stream_;
}

bool ByteStreamReader::IsCacheEmpty() const {
    return cache_size_ == MAX_CACHE_SIZE;
}

uint16_t ByteStreamReader::lastReadWord() const {
    return last_read_word_;
}

size_t ByteStreamReader::GetOffset() const {
    return static_cast<size_t>(offset_);
}
