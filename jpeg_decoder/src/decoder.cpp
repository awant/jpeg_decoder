#include <cmath>

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
    DeQuantize();
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
    while (!is_parsing_done_) {
        ParseNextSection();
    }
}

void JPGDecoder::ParseNextSection() {
    //std::cout << "before read word: " << reader_.IsCacheEmpty() << "\n";
    auto marker = reader_.ReadWord();
    std::cout << "MARKER: " << marker << "\n";
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
        case MARKER_END:
            std::cout << "DONE\n";
            is_parsing_done_ = true;
            break;
        default:
            std::cout << std::hex << "skip marker: " << std::hex << marker << std::endl;
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
    std::cout << "type: " << coeff_type << " table id: " << table_id << std::endl;

    std::vector<int> codes_counters(16);
    size_t codes_counter = 0;
    for (size_t i = 0; i < 16; ++i) {
        codes_counters[i] = reader_.ReadByte();
        codes_counter += codes_counters[i];
        std::cout << "Codes of length " << i << " bits: " << codes_counters[i] << "\n";
    }
    assert((size-16-3) == codes_counter);
    std::vector<int> values(codes_counter);
    std::cout << "codes: ";
    for (size_t i = 0; i < codes_counter; ++i) {
        values[i] = reader_.ReadByte();
        std::cout << values[i] << ", ";
    }
    std::cout << "\n";

    DHTDescriptor descriptor{table_id, coeff_type};

    HuffmanTree<int> tree(codes_counters, values);
    tree.Dump();
    huffman_trees_.emplace(descriptor, std::move(tree));
}

void JPGDecoder::ParseDQT() {
    // Quantization table
    std::cout << "--- ParseDQT ---" << std::endl;
    int size = reader_.ReadWord() - 3;
    if (size <= 0) {
        throw std::runtime_error("DHT size corrupted");
    }
    int value_size = reader_.ReadHalfByte() + 1; // size of value in table in bytes
    assert(value_size == 1 || value_size == 2);
    int table_id = reader_.ReadHalfByte();
    size_t values_count = static_cast<size_t>(size) / value_size;
    assert(values_count == 64);

    std::cout << "value size: " << values_count << " table id: " << table_id << std::endl;

    std::vector<int> values(values_count);
    for (size_t i = 0; i < values_count; ++i) {
        values[i] = value_size == 1 ? reader_.ReadByte() : reader_.ReadWord();
    }
    dqt_tables_.emplace(table_id, SquareMatrixInt::CreateFromZigZag(values, 8, 0xff));
}

void JPGDecoder::ParseSOF0() {
    std::cout << "--- ParseSOF0 ---" << std::endl;
    int size = reader_.ReadWord();
    int precision = reader_.ReadByte();
    assert(precision == 8);

    height_ = reader_.ReadWord();
    width_ = reader_.ReadWord();

    int components_count = reader_.ReadByte();
    assert(components_count == 3);

    std::cout << "precision: " << std::dec << precision
              << " size: (" << height_ << ", " << width_ << ") "
              << "components: " << components_count << std::endl;

    int max_horizontal_thinning = 0, max_vertical_thinning = 0;
    sof0_descriptors_.resize(4);

    // 1 22 0 2 11 1 3 11 1 ff
    for (size_t i = 0; i < components_count; ++i) {
        int id = reader_.ReadByte();
        std::cout << "ch id: " << id;
        assert(id < 4 && id >= 0);
        int horizontal_thinning = reader_.ReadHalfByte();
        int vertical_thinning = reader_.ReadHalfByte();

        int dqt_table_id = reader_.ReadByte();

        sof0_descriptors_[id] = ChannelDescriptor{
                horizontal_thinning,
                vertical_thinning,
                dqt_table_id
        };
        std::cout << " thinning: (" << horizontal_thinning << ", " << vertical_thinning << ") "
                  << "dqt_table_id: " << dqt_table_id << std::endl;

        max_horizontal_thinning = std::max(max_horizontal_thinning, horizontal_thinning);
        max_vertical_thinning = std::max(max_vertical_thinning, vertical_thinning);
    }
    for (size_t i = 1; i < components_count+1; ++i) {
        sof0_descriptors_[i].horizontal_thinning = max_horizontal_thinning / sof0_descriptors_[i].horizontal_thinning;
        sof0_descriptors_[i].vertical_thinning = max_vertical_thinning / sof0_descriptors_[i].vertical_thinning;
    }
}

void JPGDecoder::ParseSOS() {
    std::cout << "--- ParseSOS ---" << std::endl;
    int header_size = reader_.ReadWord();
    std::cout << "header_size: " << std::dec << header_size << "\n";
    int components_count = reader_.ReadByte();
    assert(components_count == 3);

    const int max_channel_id = 4;
    sos_descriptors_.resize(max_channel_id);

    for (int i = 0; i < components_count; ++i) {
        int channel_id = reader_.ReadByte();
        assert((channel_id < 4) && (channel_id >= 0));

        int huffman_table_dc_id = reader_.ReadHalfByte();
        int huffman_table_ac_id = reader_.ReadHalfByte();

        sos_descriptors_[channel_id] = SOSDescriptor{
                huffman_table_dc_id,
                huffman_table_ac_id
        };
        std::cout << "channel id: " << channel_id
                  << " dc table: " << huffman_table_dc_id
                  << " ac table: " << huffman_table_ac_id << "\n";
    }
    header_size -= 9;
    for (int i = 0; i < header_size; ++i) {
        reader_.ReadByte();
    }
    // Read and decode real data
    FillChannelTables();
}

void JPGDecoder::DeQuantize() {
    std::cout << "--- DeQuantize ---\n";
    for (auto& table : y_channel_tables_) {
        const auto& it = dqt_tables_.find(0);
        table.Multiply(it->second);
        std::cout << "Y channel\n";
        table.Dump();
    }
    for (auto& table : cb_channel_tables_) {
        const auto& it = dqt_tables_.find(1);
        table.Multiply(it->second);
        std::cout << "Y channel\n";
        table.Dump();
    }
    for (auto& table : cr_channel_tables_) {
        const auto& it = dqt_tables_.find(1);
        table.Multiply(it->second);
        std::cout << "Y channel\n";
        table.Dump();
    }
}

int JPGDecoder::GetNextLeafValue(HuffmanTreeInt::Iterator& huffman_tree_it) {
    int bit = 0;
    while (!huffman_tree_it.Last()) {
        bit = reader_.ReadBit();
        if (bit == 0) {
            huffman_tree_it.LeftStep();
        } else if (bit == 1) {
            huffman_tree_it.RightStep();
        } else {
            assert(false);
        }
    }
    int value = *huffman_tree_it;
    return value;
}

int JPGDecoder::NextBitsToACDCCoeff(int length) {
    int value = 0;
    bool inversed_value = false;
    for (int i = 0; i < length; ++i) {
        value = (value << 1) | reader_.ReadBit();
        if ((i == 0) && (value == 0)) { // if first bit is 0, then inverse
            inversed_value = true;
        }
    }
//    std::cout << "NextBitsToACDCCoeff\n";
//    std::cout << "value: " << value << "\ninversed_value: " << inversed_value << "\n";
    if (inversed_value) {
        value = value - std::pow(2.0, length) + 1;
    }
    return value;
}


int JPGDecoder::GetDCCoeff(int channel_id) {
    //std::cout << "--- GetDCCoeff ---\n";
    int table_id = sos_descriptors_[channel_id].huffman_table_dc_id;
    DHTDescriptor descriptor{table_id, DC};
    const auto& huffman_tree_pair_it = huffman_trees_.find(descriptor);
    assert(huffman_tree_pair_it != huffman_trees_.end());
    auto it = huffman_tree_pair_it->second.Begin();

    int value = GetNextLeafValue(it);
    if (value == 0) {
        return value;
    }
    // then value is length of coeff in bits
    int length = value;
    value = NextBitsToACDCCoeff(length);
    return value;
}

// return pair: (number of zeros, coeff value)
// number of zeros = -1 means all rest values are zeros
std::pair<int, int> JPGDecoder::GetACCoeffs(int channel_id) {
    //std::cout << "--- GetACCoeffs ---\n";
    int table_id = sos_descriptors_[channel_id].huffman_table_ac_id;
    DHTDescriptor descriptor{table_id, AC};
    const auto& huffman_tree_pair_it = huffman_trees_.find(descriptor);
    assert(huffman_tree_pair_it != huffman_trees_.end());
    auto it = huffman_tree_pair_it->second.Begin();

    uint8_t value = GetNextLeafValue(it);
    if (value == 0) {
        return std::make_pair(-1, 0);
    }
//    std::cout << "ac value: " << int(value) << " ";
    int zeros_count = value >> 4;
    int length = value & 0xf;
//    std::cout << "length: " << length << "\n";
    int result = NextBitsToACDCCoeff(length);
    return std::make_pair(zeros_count, result);
}

SquareMatrixInt JPGDecoder::GetNextChannelTable(int channel_id) {
    std::cout << "--- GetNextChannelTable ---\n";
    std::vector<int> coeffs;

    int dc_coeff = GetDCCoeff(channel_id);

    coeffs.push_back(dc_coeff);
    while ((coeffs.size() < 64)) {
        auto ac_coeffs = GetACCoeffs(channel_id);
        if (ac_coeffs.first == -1) {
            break;
        }
        for (int i = 0; i < ac_coeffs.first; ++i) {
            coeffs.push_back(0);
        }
        coeffs.push_back(ac_coeffs.second);
    }
    auto matrix = SquareMatrixInt::CreateFromZigZag(coeffs, 8, 0);

    // quantize
    int dqt_table_id = sof0_descriptors_[channel_id].dqt_table_id;
    const auto& quantize_matrix = dqt_tables_.find(dqt_table_id)->second;
    matrix.Multiply(quantize_matrix);

    return matrix;
}

void JPGDecoder::FillChannelTablesRound() {
    std::cout << "--- FillChannelTablesRound ---\n";
    // y components
    int y_components_count = 4; // horizontal = 2, vertical = 2

    int prev_dc_coeff = 0;
    for (int i = 0; i < y_components_count; ++i) {
        y_channel_tables_.emplace_back(GetNextChannelTable(1));

        auto& new_matrix = y_channel_tables_.back();

        // fix dc coeff
//        if (y_channel_tables_.size() > y_components_count) {
//            new_matrix.at(0, 0) += y_channel_tables_[y_channel_tables_.size()-y_components_count].at(0, 0);
//        }
        new_matrix.Dump();
    }
    // cb components
    int cb_components_count = 1;
    prev_dc_coeff = 0;
    for (int i = 0; i < cb_components_count; ++i) {
        cb_channel_tables_.emplace_back(GetNextChannelTable(2));

        auto& new_matrix = cb_channel_tables_.back();

        // fix dc coeff
//        if (cb_channel_tables_.size() > y_components_count) {
//            new_matrix.at(0, 0) += cb_channel_tables_[cb_channel_tables_.size()-cb_components_count].at(0, 0);
//        }
        new_matrix.Dump();
    }
    // cr components
    int cr_components_count = 1;
    prev_dc_coeff = 0;
    for (int i = 0; i < cr_components_count; ++i) {
        cr_channel_tables_.emplace_back(GetNextChannelTable(3));

        auto& new_matrix = cr_channel_tables_.back();

        // fix dc coeff
//        if (cr_channel_tables_.size() > cr_components_count) {
//            new_matrix.at(0, 0) += cr_channel_tables_[cr_channel_tables_.size() - cr_components_count].at(0, 0);
//        }
        new_matrix.Dump();
    }
}

void JPGDecoder::FillChannelTables() {
    while (!is_parsing_done_) {
        try {
            FillChannelTablesRound();
        }
        catch (const std::runtime_error&) {
            std::cout << "RUNTIME ERROR\n";
            reader_.CleanCache();
            return;
        }
    }
}

SquareMatrixDouble JPGDecoder::MakeIDCTransform(const SquareMatrixInt& matrix) {
    double coeff1 = 0.5 / 8;
    double coeff2 = std::sqrt(2);

    auto* buffer = new double[64];
    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < 8; ++j) {
            buffer[i*8+j] = matrix.at(i, j) * coeff1;
            buffer[i*8+j] = i == 0 ? buffer[i*8+j] * coeff2 : buffer[i*8+j];
            buffer[i*8+j] = j == 0 ? buffer[i*8+j] * coeff2 : buffer[i*8+j];
        }
    }
    fftw_plan plan = fftw_plan_r2r_2d(8, 8, buffer, buffer, FFTW_REDFT01, FFTW_REDFT01, 0);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    std::vector<double> resulted_buffer;
    for (int i = 0; i < 64; ++i) {
        resulted_buffer.push_back(buffer[i]);
    }
    delete[] buffer;

    auto resulted_matrix = SquareMatrixDouble(resulted_buffer, 8);
    return resulted_matrix;
}
