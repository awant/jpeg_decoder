#include "test_commons.h"

TEST_CASE("small jfif (4:2:0)", "[jpg]") {
    CheckImage("small.jpg", ":)");
}

TEST_CASE("jfif (4:4:4)", "[jpg]") {
    CheckImage("lenna.jpg");
}

TEST_CASE("jfif (4:2:2)", "[jpg]") {
    CheckImage("bad_quality.jpg", "so quality");
}

TEST_CASE("tiny jfif (4:4:4)", "[jpg]") {
    CheckImage("tiny.jpg");
}

TEST_CASE("exif (4:2:2)", "[jpg]") {
    CheckImage("chroma_halfed.jpg");
}

TEST_CASE("exif (grayscale)", "[jpg]") {
    CheckImage("grayscale.jpg");
}

TEST_CASE("jfif/exif (4:2:0)", "[jpg]") {
    CheckImage("test.jpg");
}

TEST_CASE("exif (4:4:4)", "[jpg]") {
    CheckImage("colors.jpg");
}

TEST_CASE("photoshop (4:4:4)", "[jpg]") {
    CheckImage("save_for_web.jpg");
}

TEST_CASE("Error handling", "[jpg]") {
    const size_t tests_count = 24;
    // tests: 10, 14, 18, 21, 23 (you can see these imgs on osx)
    std::vector<int> special_tests = {10, 14, 18, 21, 23};
    for (size_t i = 1; i <= tests_count; ++i) {
        if (std::find(special_tests.begin(), special_tests.end(), i) != special_tests.end()) {
            continue;
        }
        ExpectFail("bad" + std::to_string(i) + ".jpg");
    }
}
