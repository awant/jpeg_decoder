#include <cmath>
#include <tiffio.h>

#include "logger.h"
#include "decoder.h"

Image Decode(const std::string& filename) {
    std::ifstream stream(filename, std::ios_base::binary);
    if (!stream.good()) {
        throw std::runtime_error("Can't open a file: " + filename);
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

    pixel.r = clip(pixel.r, 0, 255);
    pixel.g = clip(pixel.g, 0, 255);
    pixel.b = clip(pixel.b, 0, 255);
    return pixel;
}

Image JPGDecoder::Decode() {
    if (!IsValidFormat()) {
        throw std::runtime_error("File doesn't have jpg format");
    }
    ParseJPG();
    MakeIDCTransform();
    FillChannels();
    auto image = GetRGBImage();
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
    uint16_t marker;
    try {
        marker = reader_.ReadWord();
    } catch (...) {
        is_parsing_done_ = true;
        return;
    }
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
            LOG_DEBUG << "DONE";
            is_parsing_done_ = true;
            break;
        default:
            break;
    }
}

void JPGDecoder::ParseComment() {
    LOG_DEBUG << "--- ParseComment ---";
    int comment_size = reader_.ReadWord() - COMMENT_HEADER_SIZE;
    if (comment_size <= 0) {
        throw std::runtime_error("comment size corrupted");
    }
    auto buffer = std::make_unique<char[]>(static_cast<size_t>(comment_size));
    reader_.Read(buffer.get(), static_cast<size_t>(comment_size));
    comment_ = buffer.get();
    LOG_DEBUG << "Comment: " << comment_;
}

// Huffman table
void JPGDecoder::ParseDHT() {
    LOG_DEBUG << "--- ParseDHT ---";
    int size = reader_.ReadWord();
    if (size <= 0) {
        throw std::runtime_error("DHT size corrupted");
    }
    auto coeff_type = static_cast<CoeffType>(reader_.ReadHalfByte());
    int table_id = reader_.ReadHalfByte();
    LOG_DEBUG << "type: " << coeff_type << " table id: " << table_id;

    std::vector<int> codes_counters(MAX_HUFFMAN_CODE_LEN);
    size_t codes_counter = 0;
    for (size_t i = 0; i < MAX_HUFFMAN_CODE_LEN; ++i) {
        codes_counters[i] = reader_.ReadByte();
        codes_counter += codes_counters[i];
        LOG_DEBUG << "Codes of length " << i << " bits: " << codes_counters[i];
    }
    if (static_cast<int>(DHT_HEADER_SIZE + MAX_HUFFMAN_CODE_LEN + codes_counter) != size) {
        throw std::runtime_error("DHT body corrupted");
    }
    std::vector<int> values(codes_counter);
    for (size_t i = 0; i < codes_counter; ++i) {
        values[i] = reader_.ReadByte();
    }

    DHTDescriptor descriptor{table_id, coeff_type};

    HuffmanTree<int> tree(codes_counters, values);
    LOG_DEBUG << tree;
    huffman_trees_.emplace(descriptor, std::move(tree));
}

// Table for multiplication before applying inverse DCT
// Usually is pair of table with 8x8 size
void JPGDecoder::ParseDQT() {
    LOG_DEBUG << "--- ParseDQT ---";
    int size = reader_.ReadWord() - DQT_HEADER_SIZE;
    if (size <= 0) {
        throw std::runtime_error("DQT size corrupted");
    }
    int value_size = reader_.ReadHalfByte() + 1; // size of value in table in bytes: 0 -> 1, 1 -> 2
    if ((value_size != 1) && (value_size != 2)) {
        throw std::runtime_error("DQT wrong value size");
    }
    int table_id = reader_.ReadHalfByte();
    auto values_count = static_cast<size_t>(size);
    if (values_count != kTableSide * kTableSide) {
        throw std::runtime_error("DQT wrong values count");
    }

    LOG_DEBUG << "values count: " << values_count << " table id: " << table_id;

    std::vector<int> values(values_count);
    for (size_t i = 0; i < values_count; ++i) {
        values[i] = value_size == 1 ? reader_.ReadByte() : reader_.ReadWord();
    }
    dqt_tables_.emplace(table_id, SquareMatrix<int>::CreateFromZigZag(kTableSide, values, 0xff));

    LOG_DEBUG << dqt_tables_[table_id];
}

// SOF0 is marker for base coding method (not progressive)
// Section had width and height of image and connection with thinning and QTs
void JPGDecoder::ParseSOF0() {
    LOG_DEBUG << "--- ParseSOF0 ---";
    int size = reader_.ReadWord();
    if (size <= kSOF0HeaderSize) {
        throw std::runtime_error("SOF0 size corrupted");
    }

    int precision = reader_.ReadByte();
    if (precision != kDefaultChannelPrecision) {
        throw std::runtime_error("SOF0 wrong channel precision");
    }

    height_ = reader_.ReadWord();
    width_ = reader_.ReadWord();

    uint8_t components_count = reader_.ReadByte();
    if (components_count == 0 || components_count > kComponentsCount) {
        throw std::runtime_error("SOF0 wrong components count");
    }

    LOG_DEBUG << "precision: " << std::dec << precision
              << " size: (" << height_ << ", " << width_ << ") "
              << "components: " << components_count;

    int max_horizontal_thinning = 0, max_vertical_thinning = 0;

    for (size_t component_idx = 0; component_idx < components_count; ++component_idx) {
        int channel_id = reader_.ReadByte();
        int horizontal_thinning = reader_.ReadHalfByte();
        int vertical_thinning = reader_.ReadHalfByte();
        int dqt_table_id = reader_.ReadByte();

        sof0_descriptors_[channel_id] = ChannelDescriptor{
                horizontal_thinning,
                vertical_thinning,
                0, 0,
                dqt_table_id
        };

        LOG_DEBUG << "channel id: " << channel_id
                  << " thinning: (" << horizontal_thinning << ", " << vertical_thinning << ") "
                  << "dqt_table_id: " << dqt_table_id;

        max_horizontal_thinning = std::max(max_horizontal_thinning, horizontal_thinning);
        max_vertical_thinning = std::max(max_vertical_thinning, vertical_thinning);
    }
    for (auto& elem: sof0_descriptors_) {
        elem.second.horizontal_thinning_ratio = max_horizontal_thinning / elem.second.horizontal_thinning;
        elem.second.vertical_thinning_ratio = max_vertical_thinning / elem.second.vertical_thinning;
    }
}

// Coded image
void JPGDecoder::ParseSOS() {
    LOG_DEBUG << "--- ParseSOS ---";
    uint16_t header_size = reader_.ReadWord();
    LOG_DEBUG << "header_size: " << std::dec << header_size;
    uint8_t components_count = reader_.ReadByte();
    if (components_count != sof0_descriptors_.size()) {
        throw std::runtime_error("SOS wrong components count");
    }
    channels_ids_.reserve(static_cast<int>(components_count));

    for (int component_idx = 0; component_idx < components_count; ++component_idx) {
        int channel_id = reader_.ReadByte();
        int huffman_table_dc_id = reader_.ReadHalfByte();
        int huffman_table_ac_id = reader_.ReadHalfByte();

        sos_descriptors_[channel_id] = SOSDescriptor{
                huffman_table_dc_id,
                huffman_table_ac_id
        };
        channels_ids_.push_back(channel_id);

        LOG_DEBUG << "channel id: " << channel_id
                  << " dc table: " << huffman_table_dc_id
                  << " ac table: " << huffman_table_ac_id;
    }
    if (header_size <= components_count * 2 - 3) {
        throw std::runtime_error("SOS header size corrupted");
    }
    header_size = static_cast<uint16_t>(header_size - components_count * 2 - 3);

    // Pass bytes from header
    for (int byte_idx = 0; byte_idx < header_size; ++byte_idx) {
        reader_.ReadByte();
    }
    // Read and decode real data
    FillChannelTables();
}

void JPGDecoder::FillChannelTables() {
    while (!is_parsing_done_) {
        try {
            FillChannelTablesRound();
        } catch (const std::runtime_error&) {
            // Clean current byte and check if next 2 bytes are end marker
            reader_.CleanCache();
            break;
        }
    }
}

void JPGDecoder::FillChannelTablesRound() {
    LOG_DEBUG << "--- FillChannelTablesRound ---";

    for (int channel_id: channels_ids_) {
        LOG_DEBUG << "channel_id: " << channel_id;
        int channel_parts_count = sof0_descriptors_[channel_id].vertical_thinning *
                sof0_descriptors_[channel_id].horizontal_thinning;
        LOG_DEBUG << "channel_parts_count: " << channel_parts_count;

        for (int channel_part_idx = 0; channel_part_idx < channel_parts_count; ++channel_part_idx) {
            LOG_DEBUG << "channel_part_idx: " << channel_part_idx;
            channel_tables_[channel_id].emplace_back(GetNextChannelTable(channel_id));

            // Fix DC coeff
            auto& new_matrix = channel_tables_[channel_id].back();
            LOG_DEBUG << new_matrix;
            if (channel_tables_[channel_id].size() > 1) {
                new_matrix.at(0, 0) += std::prev(channel_tables_[channel_id].end(), 2)->at(0, 0);
            }
        }
    }
}

SquareMatrixDouble JPGDecoder::GetNextChannelTable(int channel_id) {
    LOG_DEBUG << "--- GetNextChannelTable ---";
    std::vector<int> coeffs; // dc and ac coefficients
    coeffs.reserve(kTableSide*kTableSide);

    int dc_coeff = GetDCCoeff(channel_id);
    coeffs.push_back(dc_coeff);

    while ((coeffs.size() < kTableSide*kTableSide)) {
        auto ac_coeffs = GetACCoeffs(channel_id);
        if (ac_coeffs.first == -1) {
            break;
        }
        coeffs.resize(coeffs.size()+ac_coeffs.first);
        coeffs.push_back(ac_coeffs.second);
    }
    auto matrix = SquareMatrixDouble::CreateFromZigZag(kTableSide, coeffs, 0);

    // Quantize
    int dqt_table_id = sof0_descriptors_[channel_id].dqt_table_id;
    const auto& quantize_matrix = dqt_tables_[dqt_table_id];
    matrix.Multiply(quantize_matrix);

    return matrix;
}

int JPGDecoder::GetDCCoeff(int channel_id) {
    int table_id = sos_descriptors_[channel_id].huffman_table_dc_id;
    DHTDescriptor descriptor{table_id, DC};
    auto& huffman_tree = huffman_trees_[descriptor];
    auto it = huffman_tree.Begin();

    int value = GetNextLeafValue(it);
    if (!value) {
        return value;
    }
    // Then value is length of coeff in bits
    int length = value;
    value = NextBitsToACDCCoeff(length);
    return value;
}

// return pair: (number of zeros, coeff value)
// number of zeros = -1 means all rest values are zeros
std::pair<int, int> JPGDecoder::GetACCoeffs(int channel_id) {
    int table_id = sos_descriptors_[channel_id].huffman_table_ac_id;
    DHTDescriptor descriptor{table_id, AC};
    const auto& huffman_tree_pair_it = huffman_trees_.find(descriptor);
    assert(huffman_tree_pair_it != huffman_trees_.end());
    auto it = huffman_tree_pair_it->second.Begin();

    auto value = static_cast<uint8_t>(GetNextLeafValue(it));
    if (value == 0) {
        return std::make_pair(-1, 0);
    }
    int zeros_count = value >> 4;
    int length = value & 0xf;
    int result = NextBitsToACDCCoeff(length);
    return std::make_pair(zeros_count, result);
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
    if (inversed_value) {
        value -= static_cast<int>(std::pow(2.0, length)) - 1;
    }
    return value;
}

void JPGDecoder::MakeIDCTransform() {
    LOG_DEBUG << "--- MakeIDCTransform ---";
    for (int channel_id: channels_ids_) { // for every channel
        for (auto& channel_table: channel_tables_[channel_id]) {
            matrix_transformer_.MakeIDCTransform(&channel_table);
        }
    }
}

void JPGDecoder::FillChannels() {
    LOG_DEBUG << "--- FillChannels ---";
    for (int channel_id: channels_ids_) {
        LOG_DEBUG << "channel_id: " << channel_id;
        FillChannel(channel_id);
    }
}

// x <-> by height, y <-> by width
void JPGDecoder::FillChannel(int channel_id) {
    LOG_DEBUG << "--- FillChannel ---";
    int channel_width = width_ / sof0_descriptors_[channel_id].horizontal_thinning_ratio;
    int channel_height = height_ / sof0_descriptors_[channel_id].vertical_thinning_ratio;
    LOG_DEBUG << "channel_width: " << std::dec << channel_width;
    LOG_DEBUG << "channel_height: " << std::dec << channel_height;
    MatrixDouble channel(channel_height, channel_width, 0);

    const auto& channel_tables = channel_tables_[channel_id];
    int channel_table_idx = 0;
    LOG_DEBUG << "channel_tables size: " << channel_tables.size();

    int horizontal_thinning = sof0_descriptors_[channel_id].horizontal_thinning;
    int vertical_thinning = sof0_descriptors_[channel_id].vertical_thinning;

    for (int y = 0; y < channel_height; y += vertical_thinning * kTableSide) {
        for (int x = 0; x < channel_width; x += horizontal_thinning * kTableSide) {
            // map in big block (thinning) several smaller blocks
            for (int y_block_idx = 0; y_block_idx < vertical_thinning * kTableSide; y_block_idx += kTableSide) {
                for (int x_block_idx = 0; x_block_idx < horizontal_thinning * kTableSide; x_block_idx += kTableSide) {
                    Point upper_left_corner{x+x_block_idx, y+y_block_idx};
                    Point lower_right_corner{upper_left_corner.x+kTableSide, upper_left_corner.y+kTableSide};
                    channel.Map(upper_left_corner, lower_right_corner, channel_tables[channel_table_idx++]);
                }
            }
        }
    }
    channels_.emplace(channel_id, channel);
}

Image JPGDecoder::GetRGBImage() {
    LOG_DEBUG << "--- GetRGBImage ---";
    auto image = Image(width_, height_);
    image.SetComment(comment_);

    int k = 0;

    for (uint32_t y = 0; y < height_; ++y) {
        for (uint32_t x = 0; x < width_; ++x) {
            std::vector<int> yCbCr;
            for (int channel_id: channels_ids_) {
                int y_local = y / sof0_descriptors_[channel_id].vertical_thinning_ratio;
                int x_local = x / sof0_descriptors_[channel_id].horizontal_thinning_ratio;
                yCbCr.push_back(channels_[channel_id].at(y_local, x_local));
            }
            // if we have less channels numbers than kComponentsCount:
            for (size_t i = channels_ids_.size(); i < kComponentsCount; ++i) {
                yCbCr.push_back(0);
            }
            auto pixel = YCbCrToRGB(yCbCr[0], yCbCr[1], yCbCr[2]);
            image.SetPixel(y, x, pixel);
            ++k;
        }
    }
    return image;
}