#pragma once

#include <cstdint>

const int MIN_RGB_VALUE = 0;
const int MAX_RGB_VALUE = 255;

const uint16_t MARKER_JPG = 0xffd8;
const uint16_t MARKER_COMMENT = 0xfffe;
const uint16_t MARKER_DQT = 0xffdb;
const uint16_t MARKER_SOF0 = 0xffc0;
const uint16_t MARKER_DHT = 0xffc4;
const uint16_t MARKER_SOS = 0xffda;
const uint16_t MARKER_END = 0xffd9;

enum CoeffType {
    DC,
    AC
};

struct DHTDescriptor {
    int table_id;
    CoeffType table_class;
};
