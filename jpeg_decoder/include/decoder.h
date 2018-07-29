#pragma once

#include <string>
#include <istream>
#include <fstream>
#include <unordered_map>
#include <fftw3.h>

#define LOGGING_ENABLED
#include "logger/logger.h"

#include "constants.h"
#include "image.h"
#include "reader.h"
#include "huffman_tree.h"
#include "matrix.h"



// https://habr.com/post/102521/
// https://gitlab.com/slon/shad-cpp/tree/master/jpeg-decoder
Image Decode(const std::string &filename);


struct DHTDescriptorHash {
    std::size_t operator()(const DHTDescriptor &k) const {
        return std::hash<int>()((k.table_id << 1) ^ k.table_class);
    }
};

struct DHTDescriptorEqual {
    bool operator()(const DHTDescriptor &lhs, const DHTDescriptor &rhs) const {
        return lhs.table_id == rhs.table_id && lhs.table_class == rhs.table_class;
    }
};

using HuffmanTreeInt = HuffmanTree<int>;
using HuffmanMap = std::unordered_map<DHTDescriptor, HuffmanTreeInt, DHTDescriptorHash, DHTDescriptorEqual>;
using SquareMatrixInt = SquareMatrix<int>;

struct ChannelDescriptor {
    int horizontal_thinning;
    int vertical_thinning;
    int dqt_table_id;
};

struct SOSDescriptor {
    int huffman_table_dc_id;
    int huffman_table_ac_id;
};

/*

enum Channel {
    Y,
    Cb,
    Cr
};

*/


class JPGDecoder {
public:
    explicit JPGDecoder(std::istream &s);

    Image Decode();

    void Dump(std::ostream &os);

private:
    bool is_parsing_done_ = false; // set when parsing finished
    ByteStreamReader reader_;

    uint32_t height_ = 0, width_ = 0;
    std::string comment_;
    HuffmanMap huffman_trees_;
    std::unordered_map<int, SquareMatrixInt> dqt_tables_;
    std::vector<ChannelDescriptor> sof0_descriptors_;
    std::vector<SOSDescriptor> sos_descriptors_;

    std::vector<SquareMatrixInt> y_channel_tables_;
    std::vector<SquareMatrixInt> cb_channel_tables_;
    std::vector<SquareMatrixInt> cr_channel_tables_;

    bool IsValidFormat();

    void ParseJPG();
    void ParseNextSection();

    void ParseComment();
    void ParseDHT();
    void ParseDQT();
    void ParseSOF0();
    void ParseSOS();

    int NextBitsToACDCCoeff(int length);
    int GetNextLeafValue(HuffmanTreeInt::Iterator& huffman_tree_it);

    int GetDCCoeff(int channel_id);
    std::pair<int, int> GetACCoeffs(int channel_id);

    SquareMatrixInt GetNextChannelTable(int channel_id);
    void FillChannelTablesRound();
    void FillChannelTables();

    void DeQuantize();
    void IDCTransform(SquareMatrixInt* matrix);
};
