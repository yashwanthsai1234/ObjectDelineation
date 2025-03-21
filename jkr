#include "ObjectDelineation.h"
#include <cassert>
#include <cmath>
#include <sstream>
#include <iostream>

// ----- BitArray Implementation -----
BitArray::BitArray(int size) : size_(size) {
    bits_.resize((size + 63) / 64, 0);
}

bool BitArray::get(int i) const {
    int idx = i / 64;
    int offset = i % 64;
    return (bits_[idx] & (1ULL << offset)) != 0;
}

void BitArray::set(int i, bool value) {
    int idx = i / 64;
    int offset = i % 64;
    if (value)
        bits_[idx] |= (1ULL << offset);
    else
        bits_[idx] &= ~(1ULL << offset);
}

int BitArray::countOnes() const {
    int count = 0;
    for (auto word : bits_) {
        count += __builtin_popcountll(word);
    }
    return count;
}

// ----- Vertex Implementation -----
Vertex::Vertex(int x, int y, Vertex* next) : x(x), y(y), next(next), visited(false) {}

// ----- LiteList Implementation -----
LiteList::LiteList(const std::vector<short>& xs, const std::vector<short>& ys)
    : xs_(xs), ys_(ys) {}

int LiteList::numPoints() const {
    return xs_.size();
}

std::string LiteList::toString() const {
    std::ostringstream oss;
    for (size_t i = 0; i < xs_.size(); i++) {
        oss << "(" << xs_[i] << "," << ys_[i] << ")";
        if (i < xs_.size() - 1)
            oss << " -> ";
    }
    return oss.str();
}

double LiteList::area() const {
    if (xs_.empty() || xs_.front() != xs_.back() || ys_.front() != ys_.back())
        return 0.0;
    double a = 0.0;
    for (size_t i = 1; i < xs_.size(); i++) {
        a += xs_[i-1] * ys_[i] - xs_[i] * ys_[i-1];
    }
    return std::abs(a) / 2.0;
}

const std::vector<short>& LiteList::xs() const { return xs_; }
const std::vector<short>& LiteList::ys() const { return ys_; }

// ----- ObjectDelineation Implementation -----
ObjectDelineation::ObjectDelineation(int resolution, int buffer, TransformFn transform)
    : resolution_(resolution), buffer_(buffer), transform_(transform),
      occPixels_((resolution + 2 * buffer) * (resolution + 2 * buffer)) {}

int ObjectDelineation::pointOffset(int x, int y) const {
    if (x < -buffer_ || x >= resolution_ + buffer_ || y < -buffer_ || y >= resolution_ + buffer_)
        return -1;
    return (y + buffer_) * (resolution_ + 2 * buffer_) + (x + buffer_);
}

std::pair<int,int> ObjectDelineation::transformCoord(double x, double y) const {
    if (transform_)
        return transform_(x, y);
    else
        return { static_cast<int>(std::round(x)), static_cast<int>(std::round(y)) };
}

BitArray& ObjectDelineation::occupiedPixels() {
    return occPixels_;
}

std::vector<LiteList> ObjectDelineation::aggregateOccupiedPixels() {
    int scanlineSize = resolution_ + 2 * buffer_;
    std::vector<Vertex*> topVertices(scanlineSize + 1, nullptr);
    Vertex* leftVertex = nullptr;
    std::vector<Vertex*> corners;

    int offset = pointOffset(-buffer_, -buffer_);
    // Phase I: Sliding-window scan over the extended grid.
    for (int y = -buffer_; y <= resolution_ + buffer_; y++) {
        for (int x = -buffer_; x <= resolution_ + buffer_; x++) {
            int blocked0 = (x > -buffer_ && y > -buffer_ && occPixels_.get(offset - scanlineSize - 1)) ? 1 : 0;
            int blocked1 = (x < resolution_ + buffer_ && y > -buffer_ && occPixels_.get(offset - scanlineSize)) ? 2 : 0;
            int blocked2 = (x > -buffer_ && y < resolution_ + buffer_ && occPixels_.get(offset - 1)) ? 4 : 0;
            int blocked3 = (x < resolution_ + buffer_ && y < resolution_ + buffer_ && occPixels_.get(offset)) ? 8 : 0;
            int pixelType = blocked0 | blocked1 | blocked2 | blocked3;

            switch (pixelType) {
                case 0: case 3: case 5: case 10: case 12: case 15:
                    break;
                case 1: {
                    Vertex* newVertex = new Vertex(x, y, leftVertex);
                    if (topVertices[x + buffer_])
                        topVertices[x + buffer_]->next = newVertex;
                    topVertices[x + buffer_] = nullptr;
                    leftVertex = nullptr;
                    break;
                }
                case 2: {
                    Vertex* newVertex = new Vertex(x, y, topVertices[x + buffer_]);
                    leftVertex = newVertex;
                    topVertices[x + buffer_] = nullptr;
                    break;
                }
                case 4: {
                    Vertex* newVertex = new Vertex(x, y, nullptr);
                    if (leftVertex) leftVertex->next = newVertex;
                    topVertices[x + buffer_] = newVertex;
                    leftVertex = nullptr;
                    break;
                }
                case 6: {
                    Vertex* newVertex1 = new Vertex(x, y, topVertices[x + buffer_]);
                    Vertex* newVertex2 = new Vertex(x, y, nullptr);
                    if (leftVertex) leftVertex->next = newVertex2;
                    leftVertex = newVertex1;
                    topVertices[x + buffer_] = newVertex2;
                    break;
                }
                case 7: {
                    Vertex* newVertex = new Vertex(x, y, nullptr);
                    topVertices[x + buffer_] = newVertex;
                    leftVertex = newVertex;
                    corners.push_back(newVertex);
                    break;
                }
                case 8: {
                    Vertex* newVertex = new Vertex(x, y, nullptr);
                    topVertices[x + buffer_] = newVertex;
                    leftVertex = newVertex;
                    corners.push_back(newVertex);
                    break;
                }
                case 9: {
                    Vertex* newVertex1 = new Vertex(x, y, leftVertex);
                    Vertex* newVertex2 = new Vertex(x, y, nullptr);
                    if (topVertices[x + buffer_])
                        topVertices[x + buffer_]->next = newVertex1;
                    leftVertex = newVertex2;
                    topVertices[x + buffer_] = newVertex2;
                    corners.push_back(newVertex2);
                    break;
                }
                case 11: {
                    Vertex* newVertex = new Vertex(x, y, leftVertex);
                    topVertices[x + buffer_] = newVertex;
                    leftVertex = nullptr;
                    break;
                }
                case 13: {
                    Vertex* newVertex = new Vertex(x, y, nullptr);
                    if (topVertices[x + buffer_])
                        topVertices[x + buffer_]->next = newVertex;
                    leftVertex = newVertex;
                    topVertices[x + buffer_] = nullptr;
                    break;
                }
                case 14: {
                    Vertex* newVertex = new Vertex(x, y, topVertices[x + buffer_]);
                    if (leftVertex) leftVertex->next = newVertex;
                    leftVertex = nullptr;
                    topVertices[x + buffer_] = nullptr;
                    break;
                }
                default:
                    break;
            }
            offset++;
        }
        offset = offset - (resolution_ + 2 * buffer_ + 1) + scanlineSize;
    }

    // Phase II: Trace linked vertices (starting from unvisited corners) to form closed rings.
    std::vector<LiteList> rings;
    for (Vertex* corner : corners) {
        if (!corner->visited) {
            int count = 0;
            Vertex* cur = corner;
            while (cur != corner && cur != nullptr) {
                count++;
                cur = cur->next;
            }
            if (cur == nullptr)
                continue; // Incomplete loop.
            std::vector<short> xs;
            std::vector<short> ys;
            cur = corner;
            for (int i = 0; i < count; i++) {
                xs.push_back(static_cast<short>(cur->x));
                ys.push_back(static_cast<short>(cur->y));
                cur->visited = true;
                cur = cur->next;
            }
            // Close the ring.
            xs.push_back(static_cast<short>(corner->x));
            ys.push_back(static_cast<short>(corner->y));
            rings.push_back(LiteList(xs, ys));
        }
    }
    // (Note: For a complete implementation, remember to free the allocated Vertex objects.)
    return rings;
}

void ObjectDelineation::markRectangularBlock(int xStart, int xEnd, int yStart, int yEnd) {
    for (int y = yStart; y < yEnd; y++) {
        for (int x = xStart; x < xEnd; x++) {
            int offset = pointOffset(x, y);
            if (offset != -1)
                occPixels_.set(offset, true);
        }
    }
}
