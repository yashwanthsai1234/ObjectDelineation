#ifndef OBJECT_DELINEATION_H
#define OBJECT_DELINEATION_H

// --- Minimal GeometryFactory Stub ---
#ifndef GEOMETRYFACTORY_H
#define GEOMETRYFACTORY_H

#include <functional>
#include <string>
#include <vector>
#include <cstdint>
#include <utility>

class GeometryFactory {
public:
    GeometryFactory() {}

};

#endif 

class BitArray {
public:
        BitArray(int size) : size_(size) {
        bits_.resize((size + 63) / 64, 0ULL);
    }
    
    inline bool get(int i) const {
        int idx = i / 64;
        int offset = i % 64;
        return (bits_[idx] & (1ULL << offset)) != 0;
    }
    

    inline void set(int i, bool value) {
        int idx = i / 64;
        int offset = i % 64;
        if (value)
            bits_[idx] |= (1ULL << offset);
        else
            bits_[idx] &= ~(1ULL << offset);
    }
    
    inline int countOnes() const {
        int count = 0;
        for (uint64_t block : bits_) {
            count += __builtin_popcountll(block);
        }
        return count;
    }

private:
    int size_;
    std::vector<uint64_t> bits_;
};

struct Vertex {
    int x;
    int y;
    Vertex* next;  // Pointer to the next vertex in the linked list.
    bool visited;
    
    Vertex(int x, int y, Vertex* next = nullptr)
        : x(x), y(y), next(next), visited(false) {}
};

// --- LiteList Declaration (inline) ---
// Lightweight container for a traced polygon ring.
class LiteList {
public:
    LiteList(const std::vector<short>& xs, const std::vector<short>& ys)
        : xs_(xs), ys_(ys) {}
    
    int numPoints() const { return xs_.size(); }
    
    std::string toString() const {
        std::string result;
        for (size_t i = 0; i < xs_.size(); ++i) {
            result += "(" + std::to_string(xs_[i]) + "," + std::to_string(ys_[i]) + ") ";
        }
        return result;
    }
    
    // Stub for area.
    double area() const { return 0.0; }
    
    const std::vector<short>& xs() const { return xs_; }
    const std::vector<short>& ys() const { return ys_; }
    
private:
    std::vector<short> xs_;
    std::vector<short> ys_;
};


class ObjectDelineation {
public:
    using TransformFn = std::function<std::pair<int, int>(double, double)>;
    
    ObjectDelineation(int resolution, int buffer, TransformFn transform = nullptr);
    
    int pointOffset(int x, int y) const;
    
    std::pair<int, int> transformCoord(double x, double y) const;
    
    std::vector<LiteList> aggregateOccupiedPixels();
    
    void markRectangularBlock(int xStart, int xEnd, int yStart, int yEnd);
    
    BitArray& occupiedPixels();

private:
    int resolution_;
    int buffer_;
    TransformFn transform_;
    BitArray occPixels_;
};

#endif // OBJECT_DELINEATION_H
