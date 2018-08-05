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

Image Decode(const std::string& filename);


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
using SquareMatrixDouble = SquareMatrix<double>;

struct ChannelDescriptor {
    int horizontal_thinning;
    int vertical_thinning;
    int dqt_table_id;
};

struct SOSDescriptor {
    int huffman_table_dc_id;
    int huffman_table_ac_id;
};

template<class T>
T clip(T value, T min, T max) {
    return std::max(min, std::min(value, max));
}


class JPGDecoder {
public:
    explicit JPGDecoder(std::istream &s);
    Image Decode();

private:
    bool is_parsing_done_ = false; // is set when parsing finished
    ByteStreamReader reader_;  // read from stream bytes, words, bits

    MatrixTransformer<double> matrix_transformer_; // IDCT transformation

    HuffmanMap huffman_trees_;
    std::unordered_map<int, SquareMatrixInt> dqt_tables_;
    std::unordered_map<int, ChannelDescriptor> sof0_descriptors_;
    std::unordered_map<int, SOSDescriptor> sos_descriptors_;
    std::vector<int> channels_ids_;
    std::unordered_map<int, std::vector<SquareMatrixInt>> channel_tables_;

    std::vector<SquareMatrixInt> y_channel_tables_;
    std::vector<SquareMatrixInt> cb_channel_tables_;
    std::vector<SquareMatrixInt> cr_channel_tables_;

    std::vector<SquareMatrixDouble> y_channel_tables2_;
    std::vector<SquareMatrixDouble> cb_channel_tables2_;
    std::vector<SquareMatrixDouble> cr_channel_tables2_;

    // image info
    uint32_t height_ = 0, width_ = 0;
    std::string comment_;


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

    SquareMatrixDouble MakeIDCTransform(const SquareMatrixInt& matrix);
    RGB YCbCrToRGB(double Y, double Cb, double Cr);
};
