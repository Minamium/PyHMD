CXX ?= c++
NVCC ?= nvcc
MPICXX ?= mpicxx
CUDA_SM ?= 70
BUILD_DIR ?= build
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CUDAFLAGS ?= -O3 -lineinfo -Xptxas=-v
SOURCE_FINGERPRINT ?= $(shell python3 fingerprint.py)
CPPFLAGS += -DPYMHD_SOURCE_FINGERPRINT=\"$(SOURCE_FINGERPRINT)\"
HEADERS = $(shell find include -name '*.hpp')
MPI_COMPILE = $(shell $(MPICXX) --showme:compile)
MPI_LINK = $(shell $(MPICXX) --showme:link)
MPI_LIBS = $(filter -L% -l%,$(MPI_LINK))
MPI_LIB_DIRS = $(patsubst -L%,%,$(filter -L%,$(MPI_LINK)))
.PHONY: cpu cpu-mpi cuda cuda-single
cpu: $(BUILD_DIR)/pyhmd_cpu
cpu-mpi: $(BUILD_DIR)/pyhmd_cpu_mpi
cuda: $(BUILD_DIR)/pyhmd_cuda
cuda-single: $(BUILD_DIR)/pyhmd_cuda_single
$(BUILD_DIR)/pyhmd_cpu: src/main.cpp $(HEADERS) Makefile fingerprint.py
	@mkdir -p "$(BUILD_DIR)"
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude $< -o "$@"
$(BUILD_DIR)/pyhmd_cpu_mpi: src/main.cpp $(HEADERS) Makefile fingerprint.py
	@mkdir -p "$(BUILD_DIR)"
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -DPYMHD_MPI -Iinclude $< -o "$@"
$(BUILD_DIR)/pyhmd_cuda_single: src/main.cpp $(HEADERS) Makefile fingerprint.py
	@mkdir -p "$(BUILD_DIR)"
	$(NVCC) $(CPPFLAGS) -x cu $(CUDAFLAGS) -std=c++17 --extended-lambda -arch=sm_$(CUDA_SM) -DPYMHD_CUDA -Iinclude $< -o "$@"
$(BUILD_DIR)/pyhmd_cuda: src/main.cpp $(HEADERS) Makefile fingerprint.py
	@mkdir -p "$(BUILD_DIR)"
	$(NVCC) $(CPPFLAGS) -x cu $(CUDAFLAGS) -std=c++17 --extended-lambda -arch=sm_$(CUDA_SM) -DPYMHD_CUDA -DPYMHD_MPI -Iinclude $(MPI_COMPILE) $< -o "$@" $(MPI_LIBS) $(foreach dir,$(MPI_LIB_DIRS),-Xlinker=-rpath -Xlinker=$(dir))
