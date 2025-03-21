#include "ObjectDelineation.h"
#include <iostream>

int main() {
    int resolution = 1000;
    int buffer = 2;
    // Using identity transform (default rounding).
    ObjectDelineation delineator(resolution, buffer, nullptr);

    // Mark a rectangular block in the occupancy grid.
    for (int x = 100; x < 200; x++) {
        for (int y = 100; y < 150; y++) {
            int offset = delineator.pointOffset(x, y);
            if (offset != -1)
                delineator.occupiedPixels().set(offset, true);
        }
    }

    // Aggregate occupied pixels into polygon rings.
    auto rings = delineator.aggregateOccupiedPixels();

    std::cout << "Detected " << rings.size() << " ring(s):\n";
    for (const auto& ring : rings) {
        std::cout << ring.toString() << "\n";
    }

    return 0;
}
