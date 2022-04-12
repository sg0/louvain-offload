ENABLE_OMP_OFFLOAD=0
CXX = clang++
ifeq ($(ENABLE_OMP_OFFLOAD),1)
#CXX = clang++
CXX = nvc++
SM=80
endif

OPTFLAGS = -O3 -DPRINT_DIST_STATS -DPRINT_EXTRA_NEDGES -std=c++11 -DCHECK_RESULTS #-DPRINT_RESULTS
CXXFLAGS = -g -I. $(OPTFLAGS)
CXXFLAGS_THREADS = -fopenmp
ifeq ($(ENABLE_OMP_OFFLOAD),1)
CXXFLAGS_THREADS += -mp=gpu -Minfo=mp -DUSE_OMP_OFFLOAD
#CXXFLAGS_THREADS += -fopenmp-targets=nvptx64 -Xopenmp-target=nvptx64 -march=sm_${SM} -DUSE_OMP_OFFLOAD
endif

OBJ = main.o
SRC = main.cpp
TARGET = matching_omp

OBJS = $(OBJ) 
TARGETS = $(TARGET)

all: $(TARGETS)

$(TARGET):  $(OBJ)
	$(LDAPP) $(CXX) $(CXXFLAGS_THREADS) -o $@ $+ $(LDFLAGS) $(CXXFLAGS) 

$(OBJ): $(SRC)
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(CXXFLAGS_THREADS) -c $< -o $@

.PHONY: clean 

clean:
	rm -rf *~ *.dSYM nc.vg.* $(OBJS) $(TARGETS)
