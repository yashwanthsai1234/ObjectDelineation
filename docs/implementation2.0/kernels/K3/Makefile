# K3 build. CPU targets need only a C++17 compiler (clang++ on macOS, g++ on Linux); GPU targets need nvcc >= 12.8
# (sm_120). -ffp-contract=off / --fmad=false: no implicit FMA contraction (bit-equality rule 1 in k3_common.cuh).
CXX      ?= c++
NVCC     ?= nvcc
ARCH     ?= sm_120
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -ffp-contract=off -fno-fast-math
NVFLAGS  ?= -O3 -std=c++17 -arch=$(ARCH) --fmad=false -lineinfo -Xcompiler -ffp-contract=off
B        := build

.PHONY: cpu gpu all clean test-cpu
all: cpu gpu
cpu: $(B)/test_cpu $(B)/k3_cpu_cli
gpu: $(B)/test_gpu $(B)/bench_gpu

$(B):
	mkdir -p $(B)

$(B)/test_cpu: test_cpu.cpp k3_ref.cpp k3_ref.h k3_common.cuh | $(B)
	$(CXX) $(CXXFLAGS) -o $@ test_cpu.cpp k3_ref.cpp -lpthread

$(B)/k3_cpu_cli: k3_cpu_cli.cpp k3_ref.cpp k3_ref.h k3_common.cuh | $(B)
	$(CXX) $(CXXFLAGS) -o $@ k3_cpu_cli.cpp k3_ref.cpp -lpthread

$(B)/test_gpu: test_gpu.cu k3_kernels.cu k3_kernels.h k3_ref.cpp k3_ref.h k3_common.cuh | $(B)
	$(NVCC) $(NVFLAGS) -o $@ test_gpu.cu k3_kernels.cu k3_ref.cpp -lpthread

$(B)/bench_gpu: bench_gpu.cu k3_kernels.cu k3_kernels.h k3_common.cuh | $(B)
	$(NVCC) $(NVFLAGS) -o $@ bench_gpu.cu k3_kernels.cu

test-cpu: cpu
	./$(B)/test_cpu

clean:
	rm -rf $(B)
