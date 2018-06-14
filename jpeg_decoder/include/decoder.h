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


/*

enum Channel {
    Y,
    Cb,
    Cr
};


struct SOf0Descriptor {
    int id;
    int horizontal_thinning;
    int vertical_thinning;
    int dqt_table_id;
};

struct SOSDescriptor {
    int id;
    int huffman_table_dc_id;
    int huffman_table_ac_id;
};
*/


class JPGDecoder {
public:
    explicit JPGDecoder(std::istream &s);

    Image Decode();

    void Dump(std::ostream &os);

private:
    ByteStreamReader reader_;

    std::string comment_;
    HuffmanMap huffman_trees_;
    std::unordered_map<int, SquareMatrix<int>> dqt_tables_;

    bool IsValidFormat();

    void ParseJPG();
    void ParseNextSection();

    void ParseComment();
    void ParseDHT();
    void ParseDQT();
    void ParseSOF0();
    void ParseSOS();


//    int GetNextHuffmanNodeVal(HuffmanTreeInt& huffman_tree);
//    int NextBitsToCoeff(int count);
//    int GetDCCoeff(const DHTDescriptor& desc);
//    std::pair<int, int> GetACCoeffs(const DHTDescriptor& desc, bool* end);
//    Table GetNextChannelTable(int component_id);
//
//    void DeQuantize();
//
//    void FillChannelTables();
//    void FillChannelTablesRound();
//
//    bool is_eof_ = false;
//
//    uint32_t height_ = 0, width_ = 0;
//
//    std::vector<SOf0Descriptor> sof0_descriptors_;
//    std::vector<SOSDescriptor> sos_descriptors_;
//
//    std::vector<Table> y_channel_tables_;
//    std::vector<Table> cb_channel_tables_;
//    std::vector<Table> cr_channel_tables_;
};
