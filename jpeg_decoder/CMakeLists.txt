cmake_minimum_required(VERSION 3.0)
project(jpeg-decoder)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS})

find_package(JPEG REQUIRED)
include_directories(${JPEG_INCLUDE_DIR})

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})

include_directories(include)

include(../common.cmake)

add_library(decoder-lib
    src/decoder.cpp
    src/reader.cpp
)

add_executable(test_baseline
  src/test_baseline.cpp
  ../contrib/catch_main.cpp)

add_executable(test_progressive
  src/test_progressive.cpp
  ../contrib/catch_main.cpp)

add_executable(test_decoder
  src/test_decoder.cpp
  ../contrib/catch_main.cpp)

# don't forget to link them

target_link_libraries(test_baseline decoder-lib ${FFTW_LIBRARIES} ${PNG_LIBRARY} ${JPEG_LIBRARIES})
target_link_libraries(test_progressive decoder-lib ${FFTW_LIBRARIES} ${PNG_LIBRARY} ${JPEG_LIBRARIES})
target_link_libraries(test_decoder decoder-lib ${FFTW_LIBRARIES} ${PNG_LIBRARY} ${JPEG_LIBRARIES})
