#pragma once

#include <string>
#include <istream>

class ByteStreamReader {
public:
    const int MAX_CACHE_SIZE = 8;

    explicit ByteStreamReader(std::istream& stream);

    uint8_t ReadBit();
    uint8_t ReadByte();
    uint8_t ReadHalfByte();
    uint16_t ReadWord();
    void Read(char* buffer, size_t size);

    bool IsEnded() const;

private:
    std::istream& stream_;

    uint8_t cache_ = 0;
    int cache_size_ = MAX_CACHE_SIZE; // CACHE_SIZE means cache is empty

    uint8_t ReadRawByte(); // read byte directly from stream
};
