#include <catch.hpp>
#include <iostream>
#include <strstream>

#include "decoder.h"

TEST_CASE("Bytes stream reading", "[ByteStreamReader]") {
    {
        uint8_t a1 = 0xff, a2 = 0x3f, a3 = 0x56, a4 = 0x12;
        uint16_t word = static_cast<uint16_t>(0x56) << 8 | a4;
        std::strstream ss;
        ss << a1 << a2 << a3 << a4;

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

TEST_CASE("Huffman tree construction", "[HuffmanTree]") {
//    std::vector<int> counters(16, 0);
//    counters[0] = 1;
//    counters[2] = 2;
//    counters[3] = 3;
//    counters[4] = 1;
//    std::vector<int> elems = {0x01, 0x00, 0x12, 0x02, 0x11, 0x31, 0x21};
//
//    HuffmanTree<int> huffman_tree;
//
//    std::vector<char> sequence = {1,0,1};
//    bool is_decoded = false;
//    int val = huffman_tree.Decode(sequence, &is_decoded);
//    REQUIRE(is_decoded);
//    REQUIRE(val == 0x12);
}

TEST_CASE("Matrix operations", "[Matrix]") {
    const std::vector<int> values = {0,2,3,5,9,2,7,3,4};
    size_t size = 8;
    auto dht_table = SquareMatrix<int>::CreateFromZigZag(values, size, 0xff);
    std::vector<int> array{0x0, 0x2, 0x2, 0x7, 0xff, 0xff, 0xff, 0xff,
                           0x3, 0x9, 0x3, 0xff, 0xff, 0xff, 0xff, 0xff,
                           0x5, 0x4, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
    for (int i = 0; i < (64-24); ++i) {
        array.push_back(0xff);
    }
    auto true_dht_table = SquareMatrix<int>(array, 8);
    REQUIRE(dht_table == true_dht_table);
}

//TEST_CASE("Check fft tranform", "[FFT]") {
//    {
//        fftw_plan plan;
//        int n0 = 2, n1 = 2;
//
//        auto* in = new double[4];
//        auto* out = new double[4];
//
//        in[0] = 1;
//        in[1] = 2;
//        in[2] = 3;
//        in[3] = 4;
//
//        plan = fftw_plan_r2r_2d(n0, n1, in, out, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
//        fftw_execute(plan);
//        fftw_destroy_plan(plan);
//
//        std::cout << "out: ";
//        for (int i = 0; i < 4; ++i) {
//            std::cout << out[i] << " ";
//        }
//        std::cout << "\n";
//        delete[] in;
//        delete[] out;
//    }
//
//    auto* in = new int[4];
//    in[0] = 1;
//    in[1] = 2;
//    in[2] = 3;
//    in[3] = 4;
//
//    Matrix<int> matrix(in, 2, 2);
//    auto result = IDCT_JPG(matrix);
//    result.Dump();
//
//    delete[] in;
//}

TEST_CASE("Full decoder process", "[Decoder]") {
    const std::string filename = "/Users/romanmarakulin/C++/shad-cpp/jpeg-decoder/tests/small.jpg";
    Image img = Decode(filename);
}
