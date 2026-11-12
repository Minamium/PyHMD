CXX ?= c++
NVCC ?= nvcc
MPICXX ?= mpicxx
CUDA_SM ?= 70
FP_MODE ?= default
ifeq ($(FP_MODE),strict)
FP_SUFFIX := _strict
CPU_FP_FLAGS := -ffp-contract=off
CUDA_FP_FLAGS := -fmad=false
CPPFLAGS += -DPYMHD_STRICT_FP
else ifneq ($(FP_MODE),default)
$(error FP_MODE must be default or strict)
endif
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
cpu: $(BUILD_DIR)/pymhd_cpu$(FP_SUFFIX)
cpu-mpi: $(BUILD_DIR)/pymhd_cpu_mpi$(FP_SUFFIX)
cuda: $(BUILD_DIR)/pymhd_cuda$(FP_SUFFIX)
cuda-single: $(BUILD_DIR)/pymhd_cuda_single$(FP_SUFFIX)
$(BUILD_DIR)/pymhd_cpu$(FP_SUFFIX): src/main.cpp $(HEADERS) Makefile fingerprint.py
	@mkdir -p "$(BUILD_DIR)"
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CPU_FP_FLAGS) -Iinclude $< -o "$@.new"
	mv -f "$@.new" "$@"
$(BUILD_DIR)/pymhd_cpu_mpi$(FP_SUFFIX): src/main.cpp $(HEADERS) Makefile fingerprint.py
	@mkdir -p "$(BUILD_DIR)"
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) $(CPU_FP_FLAGS) -DPYMHD_MPI -Iinclude $< -o "$@.new"
	mv -f "$@.new" "$@"
$(BUILD_DIR)/pymhd_cuda_single$(FP_SUFFIX): src/main.cpp $(HEADERS) Makefile fingerprint.py
	@mkdir -p "$(BUILD_DIR)"
	$(NVCC) $(CPPFLAGS) -x cu $(CUDAFLAGS) $(CUDA_FP_FLAGS) -std=c++17 --extended-lambda -arch=sm_$(CUDA_SM) -DPYMHD_CUDA -Iinclude $< -o "$@.new"
	mv -f "$@.new" "$@"
$(BUILD_DIR)/pymhd_cuda$(FP_SUFFIX): src/main.cpp $(HEADERS) Makefile fingerprint.py
	@mkdir -p "$(BUILD_DIR)"
	$(NVCC) $(CPPFLAGS) -x cu $(CUDAFLAGS) $(CUDA_FP_FLAGS) -std=c++17 --extended-lambda -arch=sm_$(CUDA_SM) -DPYMHD_CUDA -DPYMHD_MPI -Iinclude $(MPI_COMPILE) $< -o "$@.new" $(MPI_LIBS) $(foreach dir,$(MPI_LIB_DIRS),-Xlinker=-rpath -Xlinker=$(dir))
	mv -f "$@.new" "$@"
