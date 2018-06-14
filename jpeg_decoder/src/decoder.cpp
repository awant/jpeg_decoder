#include "decoder.h"

Image Decode(const std::string& filename) {

    std::ifstream stream(filename, std::ios_base::binary);
    if (!stream.good()) {
        throw std::runtime_error("Can't open file: " + filename);
    }

    JPGDecoder decoder(stream);
    Image img = decoder.Decode();
    stream.close();
    return img;
}

JPGDecoder::JPGDecoder(std::istream& s): reader_(s) {}


Image JPGDecoder::Decode() {
    if (!IsValidFormat()) {
        throw std::runtime_error("File doesn't have jpg format");
    }
    ParseJPG();
//    // Processing
//    DeQuantize();
    return Image();
}

void JPGDecoder::Dump(std::ostream& os) {
    os << "--- Comment ---\n";
    os << comment_ << "\n";
}

bool JPGDecoder::IsValidFormat() {
    return reader_.ReadWord() == MARKER_JPG;
}

void JPGDecoder::ParseJPG() {
    for (int i = 0; i < 30; ++i) {
        ParseNextSection();
    }
}

void JPGDecoder::ParseNextSection() {
    auto marker = reader_.ReadWord();
    switch (marker) {
        case MARKER_COMMENT:
            ParseComment();
            break;
        case MARKER_DHT:
            ParseDHT();
            break;
        case MARKER_DQT:
            ParseDQT();
            break;
        case MARKER_SOF0:
            ParseSOF0();
            break;
        case MARKER_SOS:
            ParseSOS();
            break;
        default:
             throw std::runtime_error("Unknown marker");
    }
}

void JPGDecoder::ParseComment() {
    std::cout << "--- ParseComment ---" << std::endl;
    int comment_size = reader_.ReadWord() - 2;
    if (comment_size <= 0) {
        throw std::runtime_error("comments size corrupted");
    }
    auto* buffer = new char[comment_size];
    reader_.Read(buffer, static_cast<size_t>(comment_size));
    comment_ = buffer;
    // don't use unique_ptr - excessive complication
    delete[] buffer;
    std::cout << "Comment: " << comment_ << "\n";
}

void JPGDecoder::ParseDHT() {
    std::cout << "--- ParseDHT ---" << std::endl;
    int size = reader_.ReadWord();
    auto coeff_type = static_cast<CoeffType>(reader_.ReadHalfByte());
    int table_id = reader_.ReadHalfByte();

    std::vector<int> codes_counters(16);
    size_t codes_counter = 0;
    for (size_t i = 0; i < 16; ++i) {
        codes_counters[i] = reader_.ReadByte();
        codes_counter += codes_counters[i];
    }
    assert((size-16-3) == codes_counter);
    std::vector<int> values(codes_counter);
    for (size_t i = 0; i < codes_counter; ++i) {
        values[i] = reader_.ReadByte();
    }

    DHTDescriptor descriptor{table_id, coeff_type};

    HuffmanTree<int> tree(codes_counters, values);
    huffman_trees_.emplace(descriptor, std::move(tree));
}

void JPGDecoder::ParseDQT() {
    int size = reader_.ReadWord() - 3;
    if (size <= 0) {
        throw std::runtime_error("DHT size corrupted");
    }
    int value_size = reader_.ReadHalfByte() + 1; // size of value in table in bytes
    assert(value_size == 1 || value_size == 2);
    int table_id = reader_.ReadHalfByte();
    size_t values_count = static_cast<size_t>(size) / value_size;
    assert(values_count == 64);

    std::vector<int> values(values_count);
    for (size_t i = 0; i < values_count; ++i) {
        values[i] = value_size == 1 ? reader_.ReadByte() : reader_.ReadWord();
    }
    dqt_tables_.emplace(table_id, SquareMatrix<int>::CreateFromZigZag(values, 8, 0xff));
}

void JPGDecoder::ParseSOF0() {
//    std::cout << "ParseSOF0" << std::endl;
//    int size = byte_reader_.GetNextWord();
//    int precision = byte_reader_.GetNext();
//    height_ = byte_reader_.GetNextWord();
//    width_ = byte_reader_.GetNextWord();
//    std::cout << "size: " << size << ", precision: " << precision << "\n";
//    std::cout << "height: " << height_ << ", width: " << width_ << "\n";
//    int components_count = byte_reader_.GetNext();
//    if (components_count != 3) {
//        throw std::runtime_error("SOF0: unsupported components count");
//    }
//    int max_horizontal_thinning = 0, max_vertical_thinning = 0;
//    for (size_t i = 0; i < 3; ++i) {
//        sof0_descriptors_[i].id = byte_reader_.GetNext();
//        sof0_descriptors_[i].horizontal_thinning = byte_reader_.GetNextHalfByte();
//        sof0_descriptors_[i].vertical_thinning = byte_reader_.GetNextHalfByte();
//        sof0_descriptors_[i].dqt_table_id = byte_reader_.GetNext();
//        std::cout << "component" << i << ": ";
//        std::cout << sof0_descriptors_[i].id << ", ";
//        std::cout << "subsamp: " << sof0_descriptors_[i].horizontal_thinning << ", ";
//        std::cout << sof0_descriptors_[i].vertical_thinning << ", ";
//        std::cout << "qtable: " << sof0_descriptors_[i].dqt_table_id << "\n";
//        max_horizontal_thinning = std::max(max_horizontal_thinning, sof0_descriptors_[i].horizontal_thinning);
//        max_vertical_thinning = std::max(max_vertical_thinning, sof0_descriptors_[i].vertical_thinning);
//    }
//    for (size_t i = 0; i < 3; ++i) {
//        sof0_descriptors_[i].horizontal_thinning = max_horizontal_thinning / sof0_descriptors_[i].horizontal_thinning;
//        sof0_descriptors_[i].vertical_thinning = max_vertical_thinning / sof0_descriptors_[i].vertical_thinning;
//    }
}

void JPGDecoder::ParseSOS() {
//    std::cout << "ParseSOS" << std::endl;
//    int header_size = byte_reader_.GetNextWord();
//    std::cout << "header_size: " << std::dec << header_size << "\n";
//    int components_count = byte_reader_.GetNext();
//    if (components_count != 3) {
//        throw std::runtime_error("SOS: unsupported components count");
//    }
//    for (int i = 0; i < 3; ++i) {
//        int id = byte_reader_.GetNext();
//        sos_descriptors_[id].huffman_table_dc_id = byte_reader_.GetNextHalfByte();
//        sos_descriptors_[id].huffman_table_ac_id = byte_reader_.GetNextHalfByte();
//        std::cout << "comp[" << id << "]: ";
//        std::cout << "dc: " << sos_descriptors_[id].huffman_table_dc_id << " ac: ";
//        std::cout << sos_descriptors_[id].huffman_table_ac_id << "\n";
//    }
//    header_size -= 9;
//    for (int i = 0; i < header_size; ++i) {
//        byte_reader_.GetNext();
//    }
//
//    std::cout << "offset: 0x" << std::hex << is_.tellg() << "\n";
//    FillChannelTables();
}


/*
void JPGDecoder::DeQuantize() {
    // TODO: fix this! mapping through sof0_descriptors_
    for (auto& table : y_channel_tables_) {
        Multiply(table, dqt_tables_[0]);
    }
    for (auto& table : cb_channel_tables_) {
        Multiply(table, dqt_tables_[1]);
    }
    for (auto& table : cr_channel_tables_) {
        Multiply(table, dqt_tables_[1]);
    }
}

void JPGDecoder::Dump(std::ostream& os) {
    os << "Img size: width = " << width_ << " height = " << height_ << "\n";
    os << "Comment: " << comment_ << "\n";

    os << "SOF0:\n";
}

int JPGDecoder::GetNextHuffmanNodeVal(HuffmanTreeInt& huffman_tree) {
    auto it = huffman_tree.Begin();
    while (!it.Last()) {
        int bit = bit_reader_.GetNext();
        if (bit == 0) {
            it.LeftStep();
        } else {
            it.RightStep();
        }
    }
    return *it;
}

int JPGDecoder::NextBitsToCoeff(int count) {
    int coeff = 0;
    uint degree = 0;
    for (int i = 0; i < count; ++i) {
        coeff <<= 1;
        coeff |= bit_reader_.GetNext();
        degree = (degree << 1) + 1;
    }
    if (((coeff >> (count-1)) & 1) == 0) {
        coeff = -(coeff ^ degree);
    }
    return coeff;
}

int JPGDecoder::GetDCCoeff(const DHTDescriptor& desc) {
    std::cout << "GetDCCoeff\n";
    std::cout << "table_id: " << desc.table_id << ", class: " << desc.table_class << "\n";
    auto& huffman_tree = huffman_trees_.find(desc)->second;
    int node_val = GetNextHuffmanNodeVal(huffman_tree);
    std::cout << "node_val: " << node_val << "\n";
    if (node_val == 0) {
        return node_val;
    }
    int coeff = NextBitsToCoeff(node_val);
    std::cout << "dc coeff: " << coeff << "\n";
    return coeff;
}

std::pair<int, int> JPGDecoder::GetACCoeffs(const DHTDescriptor& desc, bool* end) {
    std::cout << "GetACCoeffs\n";
    std::cout << "table_id: " << desc.table_id << ", class: " << desc.table_class << "\n";
    auto& huffman_tree = huffman_trees_.find(desc)->second;
    int node_val = GetNextHuffmanNodeVal(huffman_tree);
    std::cout << "node_val: " << node_val << "\n";
    if (node_val == 0) {
        *end = true;
        return std::make_pair(node_val, node_val);
    }
    *end = false;
    int zeros_count = node_val >> 4;
    int bits_count = node_val & 0xf;
//    std::cout << "bits_count: " << bits_count << "\n";
    int coeff = NextBitsToCoeff(bits_count);
//    std::cout << "ac coeffs: " << zeros_count << " " << coeff << "\n";
    return std::make_pair(zeros_count, coeff);
}

Table JPGDecoder::GetNextChannelTable(int component_id) {
    std::cout << "GetNextChannelTable\n";
    std::vector<int> values;

    int table_id = sos_descriptors_[component_id].huffman_table_dc_id;
    std::cout << "table id: " << table_id << "\n";
    values.push_back(GetDCCoeff({table_id, DC}));

    table_id = sos_descriptors_[component_id].huffman_table_ac_id;
    bool end = false;
    while ((values.size() < 64) && !end) {
        auto ac_coeffs = GetACCoeffs({table_id, AC}, &end);
        if (end) {
            break;
        }
        for (int i = 0; i < ac_coeffs.first; ++i) {
            values.push_back(0);
        }
        values.push_back(ac_coeffs.second);
    }

    auto table = MakeZigZagTable(values, 8, 0);
    std::cout << std::dec;
    DumpTable(table);
    return table;
}

void JPGDecoder::FillChannelTablesRound() {
    // y components
    int y_components_count = 4; // horizontal = 2, vertical = 2
    int prev_dc_coeff = 0;
    for (int i = 0; i < y_components_count; ++i) {
        y_channel_tables_.emplace_back(GetNextChannelTable(1));
        if (i > 0) { // fix dc coeff
            prev_dc_coeff += y_channel_tables_.back()[0][0];
            y_channel_tables_.back()[0][0] = prev_dc_coeff;
        } else {
            prev_dc_coeff = y_channel_tables_.back()[0][0];
        }
    }
    // cb components
    int cb_components_count = 1;
    for (int i = 0; i < cb_components_count; ++i) {
        cb_channel_tables_.emplace_back(GetNextChannelTable(2));
    }
    // cr components
    int cr_components_count = 1;
    for (int i = 0; i < cr_components_count; ++i) {
        cr_channel_tables_.emplace_back(GetNextChannelTable(3));
    }
}

void JPGDecoder::FillChannelTables() {
    FillChannelTablesRound();
    std::cout << "FillChannelTables: y: " << y_channel_tables_.size() << " ";
    std::cout << "cb: " << cb_channel_tables_.size() << " cr: " << cr_channel_tables_.size() << "\n";
}

*/
