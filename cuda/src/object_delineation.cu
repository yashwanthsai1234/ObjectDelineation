#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/scan.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_set>
#include <cmath>
#include <chrono>
#include <functional>

#ifdef USE_GTEST
#include <random>
#include <gtest/gtest.h>
#endif

#ifdef STANDALONE_RUN
#include <random>
#endif

class GeometryFactory {
public:
    GeometryFactory() {}
};

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
        if (value) bits_[idx] |= (1ULL << offset);
        else bits_[idx] &= ~(1ULL << offset);
    }
    inline int countOnes() const {
        int c = 0;
        for (uint64_t b : bits_) c += __builtin_popcountll(b);
        return c;
    }
private:
    int size_;
    std::vector<uint64_t> bits_;
};

struct Vertex {
    int x;
    int y;
    Vertex* next;
    bool visited;
    Vertex(int xx, int yy, Vertex* nxt = nullptr) : x(xx), y(yy), next(nxt), visited(false) {}
};

class LiteList {
public:
    LiteList(const std::vector<short>& xs, const std::vector<short>& ys) : xs_(xs), ys_(ys) {}
    int numPoints() const { return xs_.size(); }
    std::string toString() const {
        std::ostringstream oss;
        for (size_t i = 0; i < xs_.size(); ++i) oss << "(" << xs_[i] << "," << ys_[i] << ") ";
        return oss.str();
    }
    double area() const { return 0.0; }
    const std::vector<short>& xs() const { return xs_; }
    const std::vector<short>& ys() const { return ys_; }
private:
    std::vector<short> xs_;
    std::vector<short> ys_;
};

struct Candidate {
    int x;
    int y;
    int pixelType;
    bool valid;
};

class ObjectDelineation {
public:
    using TransformFn = std::function<std::pair<int, int>(double, double)>;
    ObjectDelineation(int r, int b, TransformFn t = nullptr)
        : resolution_(r), buffer_(b), transform_(t),
          occPixels_((r + 2*b)*(r + 2*b)) {}
    int pointOffset(int x, int y) const {
        if (x < -buffer_ || x >= resolution_ + buffer_ || y < -buffer_ || y >= resolution_ + buffer_) return -1;
        return (y + buffer_)*(resolution_ + 2*buffer_) + (x + buffer_);
    }
    std::pair<int,int> transformCoord(double x, double y) const {
        if (transform_) return transform_(x, y);
        else return { (int)std::round(x), (int)std::round(y) };
    }
    BitArray& occupiedPixels() { return occPixels_; }
    void markRectangularBlock(int xs, int xe, int ys, int ye) {
        for (int y = ys; y < ye; y++) {
            for (int x = xs; x < xe; x++) {
                int off = pointOffset(x,y);
                if (off != -1) occPixels_.set(off,true);
            }
        }
    }
    std::vector<LiteList> aggregateOccupiedPixels();
private:
    int resolution_, buffer_;
    TransformFn transform_;
    BitArray occPixels_;
};

__global__ void detectCandidatesKernel(const bool* occArray, Candidate* cands, int r, int b, int w) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= w*w) return;
    int x = (idx % w) - b;
    int y = (idx / w) - b;
    int s = r + 2*b;
    int off = (y+b)*s + (x+b);
    int a0 = (x>-b && y>-b && (off-s-1)>=0 && occArray[off-s-1]) ? 1 : 0;
    int a1 = (x<r+b-1 && y>-b && occArray[off-s]) ? 2 : 0;
    int a2 = (x>-b && y<r+b-1 && occArray[off-1]) ? 4 : 0;
    int a3 = (x<r+b-1 && y<r+b-1 && occArray[off]) ? 8 : 0;
    int pt = a0|a1|a2|a3;
    bool v = true;
    switch(pt){ case 0: case 3: case 5: case 10: case 12: case 15: v=false; }
    cands[idx].x=x; cands[idx].y=y; cands[idx].pixelType=pt; cands[idx].valid=v;
}

static std::vector<LiteList> linkCandidates(const std::vector<Candidate>& cands, int w, int b, int r) {
    std::vector<Vertex*> topVertices(w, nullptr);
    Vertex* leftVertex=nullptr;
    std::vector<Vertex*> corners;
    for(int row=0; row<w; row++){
        for(int col=0; col<w; col++){
            int idx=row*w+col;
            if(!cands[idx].valid) continue;
            int x=cands[idx].x, y=cands[idx].y, pt=cands[idx].pixelType;
            switch(pt){
                case 1:{ Vertex* nv=new Vertex(x,y,leftVertex); if(topVertices[col]) topVertices[col]->next=nv; topVertices[col]=nullptr; leftVertex=nullptr; break; }
                case 2:{ Vertex* nv=new Vertex(x,y,topVertices[col]); leftVertex=nv; topVertices[col]=nullptr; break; }
                case 4:{ Vertex* nv=new Vertex(x,y,nullptr); if(leftVertex) leftVertex->next=nv; topVertices[col]=nv; leftVertex=nullptr; break; }
                case 6:{ Vertex* nv1=new Vertex(x,y,topVertices[col]); Vertex* nv2=new Vertex(x,y,nullptr); if(leftVertex) leftVertex->next=nv2; leftVertex=nv1; topVertices[col]=nv2; break; }
                case 7:{ Vertex* nv=new Vertex(x,y,nullptr); topVertices[col]=nv; leftVertex=nv; corners.push_back(nv); break; }
                case 8:{ Vertex* nv=new Vertex(x,y,nullptr); topVertices[col]=nv; leftVertex=nv; corners.push_back(nv); break; }
                case 9:{ Vertex* nv1=new Vertex(x,y,leftVertex); Vertex* nv2=new Vertex(x,y,nullptr); if(topVertices[col]) topVertices[col]->next=nv1; leftVertex=nv2; topVertices[col]=nv2; corners.push_back(nv2); break; }
                case 11:{ Vertex* nv=new Vertex(x,y,leftVertex); topVertices[col]=nv; leftVertex=nullptr; break; }
                case 13:{ Vertex* nv=new Vertex(x,y,nullptr); if(topVertices[col]) topVertices[col]->next=nv; leftVertex=nv; topVertices[col]=nullptr; break; }
                case 14:{ Vertex* nv=new Vertex(x,y,topVertices[col]); if(leftVertex) leftVertex->next=nv; leftVertex=nullptr; topVertices[col]=nullptr; break; }
                default: break;
            }
        }
    }
    std::vector<LiteList> rings;
    for (Vertex* c : corners) {
        if(!c->visited){
            std::vector<short> xs; std::vector<short> ys;
            Vertex* s=c; Vertex* cur=c;
            do{ xs.push_back((short)cur->x); ys.push_back((short)cur->y); cur->visited=true; cur=cur->next; } while(cur && cur!=s);
            xs.push_back((short)s->x); ys.push_back((short)s->y);
            rings.push_back(LiteList(xs, ys));
        }
    }
    std::unordered_set<Vertex*> freed;
    for(Vertex* c: corners){
        if(!c) continue;
        if(freed.find(c)!=freed.end()) continue;
        Vertex* tmp=c;
        do{ freed.insert(tmp); tmp=tmp->next; } while(tmp && tmp!=c);
    }
    for(auto*v:freed) delete v;
    return rings;
}

std::vector<LiteList> ObjectDelineation::aggregateOccupiedPixels() {
    int w=resolution_+2*buffer_; 
    int n=w*w; 
    int g=(resolution_+2*buffer_)*(resolution_+2*buffer_);
    bool* occ; 
    cudaMallocManaged(&occ,g*sizeof(bool));
    for(int i=0;i<g;i++) occ[i]=occPixels_.get(i);
    Candidate* d_cand; 
    cudaMallocManaged(&d_cand,n*sizeof(Candidate));
    int tpb=256; 
    int bl=(n+tpb-1)/tpb;
    detectCandidatesKernel<<<bl,tpb>>>(occ, d_cand, resolution_, buffer_, w);
    cudaDeviceSynchronize();
    std::vector<Candidate> cands(n);
    for(int i=0;i<n;i++) cands[i]=d_cand[i];
    cudaFree(d_cand);
    cudaFree(occ);
    return linkCandidates(cands,w,buffer_,resolution_);
}

#ifdef STANDALONE_RUN
int main(){
    int resolution=1000; 
    int buffer=0; 
    ObjectDelineation od(resolution, buffer, nullptr);
    int total=(resolution+2*buffer)*(resolution+2*buffer);
    for(int i=0;i<total;i++){ od.occupiedPixels().set(i,(i%2==0)); }
    auto st=std::chrono::high_resolution_clock::now();
    auto rings=od.aggregateOccupiedPixels();
    auto en=std::chrono::high_resolution_clock::now();
    long long ms=std::chrono::duration_cast<std::chrono::milliseconds>(en-st).count();
    std::cout<<"CUDA sliding window + serial linking took "<<ms<<" ms\n";
    std::cout<<"Number of rings detected: "<<rings.size()<<"\n";
    return 0;
}
#endif

#ifdef USE_GTEST
#include <random>
#include <gtest/gtest.h>

TEST(ObjectDelineationTest, PerformanceStressTest){
    int r=1000,b=0;
    ObjectDelineation od(r,b,nullptr);
    std::mt19937 rng(0);
    std::uniform_real_distribution<double> dist(0.0,1.0);
    int t=(r+2*b)*(r+2*b);
    for(int i=0;i<t;i++){
        if(dist(rng)<0.5) od.occupiedPixels().set(i,true);
    }
    auto st=std::chrono::high_resolution_clock::now();
    auto rings=od.aggregateOccupiedPixels();
    auto en=std::chrono::high_resolution_clock::now();
    long long ms=std::chrono::duration_cast<std::chrono::milliseconds>(en-st).count();
    std::cout<<"Aggregation took "<<ms<<" ms\n";
    EXPECT_GT(rings.size(),0);
}

TEST(ObjectDelineationTest, PointOutside){
    int r=100,b=2;
    ObjectDelineation od(r,b,nullptr);
    int off=od.pointOffset(-10,-10);
    EXPECT_EQ(off,-1);
}

int main(int argc,char**argv){
    ::testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
}
#endif
