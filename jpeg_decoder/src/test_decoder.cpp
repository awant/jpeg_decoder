#include <catch.hpp>
#include <iostream>
#include <strstream>
#include "test_commons.h"
#include "decoder.h"


TEST_CASE("Bytes stream reading", "[ByteStreamReader]") {
    {
        uint8_t aff = 0xff, a3f = 0x3f, a56 = 0x56, a12 = 0x12;
        uint16_t word = static_cast<uint16_t>(0x56) << 8 | a12;
        std::strstream ss;
        ss << aff << a3f << a56 << a12;

        ByteStreamReader reader(ss);
        REQUIRE(reader.ReadByte() == 0xff);
        REQUIRE(reader.ReadHalfByte() == 0x3);
        for (int i = 0; i < 4; ++i) {
            REQUIRE(reader.ReadBit() == 0x1);
            REQUIRE(!reader.IsEnded());
        }
        REQUIRE(reader.ReadWord() == word);
    }
}

TEST_CASE("Matrix operations", "[Matrix]") {
    const std::vector<int> values = {0, 2, 3, 5, 9, 2, 7, 3, 4};
    size_t size = 8;
    auto dht_table = SquareMatrix<int>::CreateFromZigZag(size, values, 0xff);
    std::vector<int> array{0x0, 0x2, 0x2, 0x7, 0xff, 0xff, 0xff, 0xff,
                           0x3, 0x9, 0x3, 0xff, 0xff, 0xff, 0xff, 0xff,
                           0x5, 0x4, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
    array.resize(size*size, 0xff);

    auto true_dht_table = SquareMatrix<int>(size, array);

    REQUIRE(dht_table == true_dht_table);

    auto matrix1 = Matrix<double>(2, 3, {
        1.2, 2.5, 3.7,
        4.8, 5.3, 6.4
    });
    auto matrix2_sum = Matrix<double>(2, 3, {
        0.1, 0.3, 1.2,
        9.1, 5.9, 4.1
    });
    auto matrix12_sum = Matrix<double>(2, 3, {
        1.3, 2.8, 4.9,
        13.9, 11.2, 10.5
    });
    matrix1.Add(matrix2_sum);
    REQUIRE(matrix1 == matrix12_sum);

    auto matrix2_mul = Matrix<double>(2, 3, {
            2.0, 1.0, 4.0,
            1.0, 0.5, 2
    });
    auto matrix12_mul = Matrix<double>(2, 3, {
            2.6, 2.8, 19.6,
            13.9, 5.6, 21.0
    });
    matrix1.Multiply(matrix2_mul);
    REQUIRE(matrix1 == matrix12_mul);

    double coeff = 0.5;
    matrix1.Dot(coeff);
    auto matrix12_dot = Matrix<double>(2, 3, {
            1.3, 1.4, 9.8,
            6.95, 2.8, 10.5
    });
    REQUIRE(matrix1 == matrix12_dot);
}

TEST_CASE("Huffman tree construction", "[HuffmanTree]") {
    std::vector<int> codes_counters{1, 1};
    codes_counters.resize(16);
    std::vector<int> values = {0x03, 0x02};

    HuffmanTree<int> huffman_tree(codes_counters, values);

    std::vector<char> sequence = {1, 0};
    auto sequence_it = sequence.begin();
    auto huffman_tree_it = huffman_tree.Begin();

    int bit = 0;
    while (!huffman_tree_it.Last()) {
        bit = *sequence_it;
        if (bit == 0) {
            huffman_tree_it.LeftStep();
        } else if (bit == 1) {
            huffman_tree_it.RightStep();
        } else {
            REQUIRE(false);
        }
        if (sequence_it != sequence.end()) {
            ++sequence_it;
        }
    }
    int value = *huffman_tree_it;
    REQUIRE(value == 0x02);
}

TEST_CASE("FFT Transformation", "[FFT]") {
    MatrixTransformer<double> matrix_transformer;

    std::vector<double> values{
        2.4, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
    };
    values.resize(8*8);
    SquareMatrixDouble matrix(8, values);

    SquareMatrixDouble idc_matrix(8, 0.3);

    matrix_transformer.MakeIDCTransform(&matrix);

    REQUIRE(matrix == idc_matrix);
}

TEST_CASE("Full decoder process", "[Decoder]") {
    const std::string filename = "/Users/romanmarakulin/C++/shad-cpp/jpeg-decoder/tests/small.jpg";
    Image img = Decode(filename);

    auto dot_pos = filename.find('.');
    std::string png_filename(filename.substr(0, dot_pos) + ".png");
    WritePng(png_filename, img);
}
