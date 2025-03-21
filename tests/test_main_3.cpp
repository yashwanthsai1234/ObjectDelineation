// tests/test_main_2.cpp
#include "ObjectDelineation.h"  // Header for your ObjectDelineation, BitArray, Vertex, LiteList, etc.
#include <gtest/gtest.h>
#include <chrono>
#include <random>
#include <iostream>

// Test for performance stress of aggregateOccupiedPixels
TEST(ObjectDelineationTest, PerformanceStressTest) {
    int resolution = 1000;
    int buffer = 0;  // No extra buffer to stress the grid.
    // Call the constructor with three arguments (the third is a transformation function; use nullptr for identity)
    ObjectDelineation od(resolution, buffer, nullptr);

    std::mt19937 rng(0);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    int totalPixels = (resolution + 2 * buffer) * (resolution + 2 * buffer);

    // Fill each pixel with probability 0.5.
    for (int i = 0; i < totalPixels; i++) {
        if (dist(rng) < 0.5)
            od.occupiedPixels().set(i, true);
    }

    auto start = std::chrono::high_resolution_clock::now();
    auto rings = od.aggregateOccupiedPixels();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Aggregation took " << duration_ms << " ms" << std::endl;
    EXPECT_GT(rings.size(), 0);
}

// Test that a point outside the valid grid returns an offset of -1.
TEST(ObjectDelineationTest, PointOutside) {
    int resolution = 100;
    int buffer = 2;
    ObjectDelineation od(resolution, buffer, nullptr);
    int offset = od.pointOffset(-10, -10);
    EXPECT_EQ(offset, -1);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
