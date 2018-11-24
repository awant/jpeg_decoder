cmake_minimum_required(VERSION 3.0)
project(jpeg-decoder)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 14)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS})

find_package(JPEG REQUIRED)
include_directories(${JPEG_INCLUDE_DIR})

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})

include_directories(include)
include_directories(logger)

add_library(decoder-lib
    src/decoder.cpp
    src/reader.cpp
)

add_executable(test_baseline
  src/test_baseline.cpp
        )

target_link_libraries(test_baseline decoder-lib ${FFTW_LIBRARIES} ${PNG_LIBRARY} ${JPEG_LIBRARIES})
