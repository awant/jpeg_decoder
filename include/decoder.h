#pragma once

#include <string>
#include <istream>
#include <fstream>
#include <unordered_map>
#include <fftw3.h>

#include "constants.h"
#include "image.h"
#include "reader.h"
#include "huffman_tree.h"
#include "matrix.h"


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
using MatrixDouble = Matrix<double>;
using SquareMatrixDouble = SquareMatrix<double>;

struct ChannelDescriptor {
    int horizontal_thinning;
    int vertical_thinning;
    int horizontal_thinning_ratio;
    int vertical_thinning_ratio;
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
    std::unordered_map<int, SquareMatrix<int>> dqt_tables_;
    std::unordered_map<int, ChannelDescriptor> sof0_descriptors_;
    std::unordered_map<int, SOSDescriptor> sos_descriptors_;
    std::vector<int> channels_ids_;
    std::unordered_map<int, std::vector<SquareMatrixDouble>> channel_tables_;
    std::unordered_map<int, MatrixDouble> channels_;

    // image info
    uint32_t height_ = 0, width_ = 0;
    std::string comment_;

    bool IsValidFormat();

    void ParseJPG();
    void MakeIDCTransform();
    void FillChannels();
    void FillChannel(int channel_id);
    Image GetRGBImage();

    void ParseNextSection();
    void ParseComment();
    void ParseDHT();
    void ParseDQT();
    void ParseSOF0();
    void ParseSOS();
    void ParseAPP1();

    std::string GetSectionOffset() const;

    int NextBitsToACDCCoeff(int length);
    int GetNextLeafValue(HuffmanTreeInt::Iterator& huffman_tree_it);
    int GetDCCoeff(int channel_id);
    std::pair<int, int> GetACCoeffs(int channel_id);
    SquareMatrixDouble GetNextChannelTable(int channel_id);
    void FillChannelTablesRound();
    void FillChannelTables();

    RGB YCbCrToRGB(double Y, double Cb, double Cr);
};
