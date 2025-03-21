#include "ObjectDelineation.h"
#include <gtest/gtest.h>

TEST(ObjectDelineationTest, RectangularBlock) {
    int resolution = 100;
    int buffer = 2;
    ObjectDelineation delineator(resolution, buffer, nullptr);
    delineator.markRectangularBlock(10, 20, 10, 15);
    auto rings = delineator.aggregateOccupiedPixels();
    // Expect at least one ring to be detected.
    EXPECT_GT(rings.size(), 0);
}

TEST(ObjectDelineationTest, PointOutside) {
    int resolution = 100;
    int buffer = 2;
    ObjectDelineation delineator(resolution, buffer, nullptr);
    int offset = delineator.pointOffset(-10, -10);
    EXPECT_EQ(offset, -1);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
