#include <cmath>
#include <tiffio.h>

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

RGB JPGDecoder::YCbCrToRGB(double Y, double Cb, double Cr) {
    RGB pixel;
    pixel.r = static_cast<int>(Y + 1.402 * Cr + 128);
    pixel.g = static_cast<int>(Y - 0.34414 * Cb - 0.71414 * Cr + 128);
    pixel.b = static_cast<int>(Y + 1.772 * Cb + 128);

    pixel.r = std::max(0, std::min(pixel.r, 255));
    pixel.r = std::max(0, std::min(pixel.g, 255));
    pixel.r = std::max(0, std::min(pixel.b, 255));
    return pixel;
}

Image JPGDecoder::Decode() {
    if (!IsValidFormat()) {
        throw std::runtime_error("File doesn't have jpg format");
    }
    ParseJPG();

    for (const auto& table: y_channel_tables_) {
        y_channel_tables2_.emplace_back(MakeIDCTransform(table));
    }
    for (const auto& table: cb_channel_tables_) {
        cb_channel_tables2_.emplace_back(MakeIDCTransform(table));
    }
    for (const auto& table: cr_channel_tables_) {
        cr_channel_tables2_.emplace_back(MakeIDCTransform(table));
    }
//    std::cout << "y_channel_tables_.front\n";
//    y_channel_tables_.front().Dump();
//    std::cout << "y_channel_tables2_.front\n";
//    y_channel_tables2_.front().Dump();
//    std::cout << "cb_channel_tables2_.front\n";
//    cb_channel_tables2_.front().Dump();
//    std::cout << "cr_channel_tables2_.front\n";
//    cr_channel_tables2_.front().Dump();

    std::cout << ": " << y_channel_tables2_.size() << "\n";
    std::cout << ": " << cb_channel_tables2_.size() << "\n";
    std::cout << ": " << cr_channel_tables2_.size() << "\n";
    std::cout << ": " << height_ << "x" << width_ << "\n";

    std::cout << "--- fill y channel --- \n";
    std::vector<std::vector<int>> y_channel(32, std::vector<int>(32, 0));
    // TODO: fix
    for (size_t i = 0; i < y_channel_tables2_.size(); ++i) {
        int block_row = i / 8 * 2 + (i % 4) / 2;
        int block_col = i / 4 % 2 * 2 + i % 4 % 2;

        for (int y = block_row * 8; y < (block_row + 1) * 8; ++y) {
            for (int x = block_col * 8; x < (block_col + 1) * 8; ++x) {
                y_channel[y][x] = y_channel_tables2_[i].at(y - block_row * 8, x - block_col * 8);
            }
        }
    }

    std::cout << "--- fill others channels --- \n";
    std::vector<std::vector<int>> cb_channel(16, std::vector<int>(16, 0));
    std::vector<std::vector<int>> cr_channel(16, std::vector<int>(16, 0));
    for (size_t i = 0; i < cb_channel_tables2_.size(); ++i) {
        int block_row = i / 2;
        int block_col = i % 2;

        for (int y = block_row * 8; y < (block_row + 1) * 8; ++y) {
            for (int x = block_col * 8; x < (block_col + 1) * 8; ++x) {
                cb_channel[y][x] = cb_channel_tables2_[i].at(y - block_row * 8, x - block_col * 8);
                cr_channel[y][x] = cr_channel_tables2_[i].at(y - block_row * 8, x - block_col * 8);
            }
        }
    }

    std::cout << "--- create image ---\n";

    auto image = Image(width_, height_);
    image.SetComment(comment_);

    for (int y = 0; y < 32; ++y) {
        for (int x = 0; x < 32; ++x) {
            auto pixel = YCbCrToRGB(y_channel[y][x], cb_channel[y/2][x/2], cr_channel[y/2][x/2]);
            image.SetPixel(y, x, pixel);
        }
    }

    return image;
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
    dqt_tables_.emplace(table_id, SquareMatrixInt::CreateFromZigZag(8, values, 0xff));
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

void JPGDecoder::FillChannelTables() {
    while (!is_parsing_done_) {
        try {
            FillChannelTablesRound();
        }
        catch (const std::runtime_error&) {
            // Clean current byte and check if next 2 bytes is end marker
            reader_.CleanCache();
            break;
        }
    }
    std::cout << "last tables\n";
    y_channel_tables_.back().Dump();
    cb_channel_tables_.back().Dump();
    cr_channel_tables_.back().Dump();
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
    if (!value) {
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
    auto matrix = SquareMatrixInt::CreateFromZigZag(8, coeffs, 0);

    // Quantize - multiplying
    int dqt_table_id = sof0_descriptors_[channel_id].dqt_table_id;
    const auto& quantize_matrix = dqt_tables_.find(dqt_table_id)->second;
    matrix.Multiply(quantize_matrix);

    return matrix;
}

void JPGDecoder::FillChannelTablesRound() {
    std::cout << "--- FillChannelTablesRound ---\n";
    // y components
    int y_components_count = 4; // horizontal = 2, vertical = 2

    for (int i = 0; i < y_components_count; ++i) {
        y_channel_tables_.emplace_back(GetNextChannelTable(1));

        auto& new_matrix = y_channel_tables_.back();
        new_matrix.Dump();

        // Fix DC coeff
        if (y_channel_tables_.size() > 1) {
            new_matrix.at(0, 0) += y_channel_tables_[y_channel_tables_.size()-2].at(0, 0);
        }
    }
    // cb components
    int cb_components_count = 1;
    for (int i = 0; i < cb_components_count; ++i) {
        cb_channel_tables_.emplace_back(GetNextChannelTable(2));

        auto& new_matrix = cb_channel_tables_.back();
        new_matrix.Dump();

        // Fix DC coeff
        if (cb_channel_tables_.size() > 1) {
            new_matrix.at(0, 0) += cb_channel_tables_[cb_channel_tables_.size()-2].at(0, 0);
        }
    }
    // cr components
    int cr_components_count = 1;
    for (int i = 0; i < cr_components_count; ++i) {
        cr_channel_tables_.emplace_back(GetNextChannelTable(3));

        auto& new_matrix = cr_channel_tables_.back();
        new_matrix.Dump();

        // Fix DC coeffs
        if (cr_channel_tables_.size() > 1) {
            new_matrix.at(0, 0) += cr_channel_tables_[cr_channel_tables_.size()-2].at(0, 0);
        }
    }
}

SquareMatrixDouble JPGDecoder::MakeIDCTransform(const SquareMatrixInt& matrix) {
    std::cout << "--- MakeIDCTransform ---\n";

    SquareMatrixDouble result_matrix(matrix);
    matrix_transformer_.MakeIDCTransform(&result_matrix);

    return result_matrix;
}
