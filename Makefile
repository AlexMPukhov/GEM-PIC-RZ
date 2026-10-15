# Plain Makefile (alternative to CMake):   make, make test, make clean
# GPU build:  make CUDA=1 [NVCC=/path/to/nvcc] [CUDA_ARCH="-arch=sm_80"]
# HDF5 output is enabled if pkg-config finds hdf5 (make HDF5=0 disables it)
CXX      ?= g++
# CXX=nvcc (e.g. from a Kokkos setup): nvcc is for engine.cpp only, use g++ for C++
ifeq ($(notdir $(CXX)),nvcc)
  NVCC := $(CXX)
  CXX  := g++
endif
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
BUILD    := build-make
CORE     := units deck config mapping atomic_data ionization grid solver smoothing simulation diagnostics output_hdf5 engine comm pipeline eigen expr beam_import
OBJ      := $(CORE:%=$(BUILD)/%.o)
# MPI (domain decomposition along zeta): enabled if mpicxx is found (make MPI=0 disables it)
MPI ?= $(shell command -v mpicxx > /dev/null 2>&1 && echo 1)
ifeq ($(MPI),1)
  CXX      := mpicxx
  CXXFLAGS += -DGEMPIC_HAVE_MPI
endif
HDF5 ?= $(shell pkg-config --exists hdf5 && echo 1)
ifeq ($(HDF5),1)
  CXXFLAGS += -DGEMPIC_HAVE_HDF5 $(shell pkg-config --cflags hdf5)
  LDLIBS   += $(shell pkg-config --libs hdf5)
  TOOLS    := $(BUILD)/h5_extract
endif
ifeq ($(CUDA),1)
  NVCC      ?= nvcc
  CUDA_ARCH ?= -gencode arch=compute_80,code=sm_80 -gencode arch=compute_90,code=sm_90
  CXXFLAGS  += -DGEMPIC_HAVE_CUDA
  LDLIBS    += -L$(dir $(shell which $(NVCC)))../lib64 -L$(dir $(shell which $(NVCC)))../lib -lcudart_static -ldl -lrt -lpthread
else
  CXXFLAGS  += -fopenmp
endif

all: $(BUILD)/gempic_rz $(BUILD)/test_deck $(BUILD)/compare_columns $(BUILD)/check_physics $(TOOLS)

$(BUILD)/%.o: src/%.cpp src/*.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

ifeq ($(CUDA),1)
$(BUILD)/engine.o: src/engine.cpp src/*.hpp | $(BUILD)
	$(NVCC) -x cu -O3 -std=c++17 $(CUDA_ARCH) --extended-lambda --expt-relaxed-constexpr -DGEMPIC_HAVE_CUDA -Isrc -c $< -o $@
endif

$(BUILD)/gempic_rz: $(OBJ) src/main.cpp
	$(CXX) $(CXXFLAGS) -Isrc src/main.cpp $(OBJ) -o $@ $(LDLIBS)

$(BUILD)/test_deck: $(OBJ) tests/test_deck.cpp
	$(CXX) $(CXXFLAGS) -Isrc -DGEM_SOURCE_DIR=\"$(CURDIR)\" tests/test_deck.cpp $(OBJ) -o $@ $(LDLIBS)

$(BUILD)/compare_columns: tools/compare_columns.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD)/check_physics: $(OBJ) tests/check_physics.cpp
	$(CXX) $(CXXFLAGS) -Isrc tests/check_physics.cpp $(OBJ) -o $@ $(LDLIBS)

$(BUILD)/h5_extract: tools/h5_extract.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

$(BUILD):
	mkdir -p $(BUILD)

ifeq ($(CUDA),1)
  PAR := cuda
else
  PAR := openmp
endif

test: all
	$(BUILD)/test_deck
	sh tests/regression.sh $(BUILD)/gempic_rz $(BUILD)/compare_columns numerics.backend=serial
	sh tests/physics.sh $(BUILD)/gempic_rz $(BUILD)/check_physics numerics.backend=serial
	sh tests/regression.sh $(BUILD)/gempic_rz $(BUILD)/compare_columns numerics.backend=$(PAR)
	sh tests/physics.sh $(BUILD)/gempic_rz $(BUILD)/check_physics numerics.backend=$(PAR)
ifeq ($(HDF5),1)
	sh tests/hdf5.sh $(BUILD)/gempic_rz $(BUILD)/h5_extract $(BUILD)/compare_columns
	sh tests/restart.sh $(BUILD)/gempic_rz $(BUILD)/compare_columns serial
	sh tests/restart.sh $(BUILD)/gempic_rz $(BUILD)/compare_columns $(PAR)
	sh tests/cherenkov.sh $(BUILD)/gempic_rz $(BUILD)/h5_extract
	sh tests/undulator.sh $(BUILD)/gempic_rz $(BUILD)/h5_extract
	sh tests/flexgrid.sh $(BUILD)/gempic_rz $(BUILD)/compare_columns $(PAR)
	sh tests/functions.sh $(BUILD)/gempic_rz $(BUILD)/compare_columns $(PAR)
	sh tests/dopant.sh $(BUILD)/gempic_rz $(BUILD)/compare_columns $(PAR)
	sh tests/smoothing.sh $(BUILD)/gempic_rz $(BUILD)/compare_columns $(PAR)
	sh tests/transport.sh $(BUILD)/gempic_rz $(BUILD)/h5_extract $(PAR)
	sh tests/beam_import.sh $(BUILD)/gempic_rz $(BUILD)/compare_columns
	sh tests/openpmd.sh $(BUILD)/gempic_rz $(BUILD)/h5_extract || [ $$? -eq 77 ]
ifeq ($(MPI),1)
	sh tests/mpi.sh $(BUILD)/gempic_rz $(BUILD)/compare_columns mpirun
endif
endif

clean:
	rm -rf $(BUILD)

.PHONY: all test clean
