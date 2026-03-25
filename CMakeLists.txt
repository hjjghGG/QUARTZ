


cmake_minimum_required(VERSION 3.16)
project(test_first C CXX)  # C C++

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# compile setting
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/UDDSketch/include
        ${PROJECT_SOURCE_DIR}/UDDSketch/include/dict
)

add_subdirectory(UDDSketch)

add_executable(test_first
        test.cpp
        KC_pu.cc
)

target_link_libraries(test_first PRIVATE
        UDDSketch
        m
        pthread
)

target_include_directories(test_first PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/UDDSketch/include
        ${PROJECT_SOURCE_DIR}/UDDSketch/include/dict
)