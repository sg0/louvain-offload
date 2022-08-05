#===============================================================================
# Compiler Options
#===============================================================================

COMPILER    = armcpu
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
SM = cc70   # --- NVIDIA arch
ARCH = gfx90a # --- AMD arch
ENABLE_OMP_OFFLOAD = 0
SAVE_TEMP = 0

#===============================================================================
# Program name & source code list
#===============================================================================

OBJ = main.o
SRC = main.cpp
TARGET = louvain_omp_$(COMPILER)

#===============================================================================
# Sets Flags
#===============================================================================

# Standard Flags
CFLAGS := -std=c++11 -Wall

# Linker Flags
LDFLAGS = -lm

OPTFLAGS = -DPRINT_DIST_STATS -DPRINT_EXTRA_NEDGES

# ARM Compiler
ifeq ($(COMPILER),armcpu)
  CC = armclang++
  CFLAGS +=-Ofast -fopenmp -mcpu=a64fx #-DZFILL_CACHE_LINES
endif

ifeq ($(COMPILER),gnucpu)
  CC = g++
  CFLAGS +=-O3 -fopenmp #-DZFILL_CACHE_LINES
endif

# GCC Compiler
ifeq ($(COMPILER),gnu)
  CC = g++
ifeq ($(ENABLE_OMP_OFFLOAD),1)
  CFLAGS += -O3 -fopenmp -flto
endif
endif

# Intel Compiler
ifeq ($(COMPILER),intel)
  CC = icpx 
ifeq ($(ENABLE_OMP_OFFLOAD),1)
  CFLAGS += -fiopenmp -fopenmp-targets=spir64 -D__STRICT_ANSI__ 
endif
endif

# LLVM Clang Compiler 
ifeq ($(COMPILER),llvm_nv)
  CC = clang++
ifeq ($(ENABLE_OMP_OFFLOAD),1)
  CFLAGS += -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda  -fopenmp-cuda-mode 
  #CFLAGS += -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target --cuda-path=${OLCF_CUDA_ROOT}  -Xcuda-ptxas --maxrregcount=60 -fopenmp-cuda-mode
  #CFLAGS += -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target --cuda-path=${OLCF_CUDA_ROOT}    -fopenmp-new-driver -foffload-lto 
  #CFLAGS += -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target --cuda-path=${OLCF_CUDA_ROOT}  -Xcuda-ptxas --maxrregcount=120 -fopenmp-new-driver -foffload-lto -fopenmp-assume-no-thread-state
endif
endif

# IBM XL Compiler
ifeq ($(COMPILER),xl)
  CC = xlC
ifeq ($(ENABLE_OMP_OFFLOAD),1)
  CFLAGS += -qsmp=omp -qoffload -qstrict
endif
endif

# NVIDIA NVHPC Compiler 
ifeq ($(COMPILER),nvhpc)
  CC = nvc++
ifeq ($(ENABLE_OMP_OFFLOAD),1)
  #CFLAGS += -mp=gpu -gpu=managed
  CFLAGS += -mp=gpu -gpu=${SM}
  #CFLAGS += -mp=gpu -Minfo=mp -gpu=${SM}
endif
endif

# AOMP Compiler
ifeq ($(COMPILER),llvm_amd)
  CC = clang++
ifeq ($(ENABLE_OMP_OFFLOAD),1)
  CFLAGS += -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=${ARCH}
endif
endif

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
  LDFLAGS  += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg
  LDFLAGS  += -pg
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += #-O3
endif

# Using device offload
ifeq ($(ENABLE_OMP_OFFLOAD),1)
  CFLAGS += -DUSE_OMP_OFFLOAD
else
  CFLAGS += -fopenmp -DGRAPH_FT_LOAD=4 #-I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/
endif

# Compiler Trace  
ifeq ($(SAVE_TEMPS),1)
CFLAGS += -save-temps
endif


# Team size
ifeq ($(ENABLE_OMP_OFFLOAD),1)
	TS=32
ifeq ($(TS), 16)
  OPTFLAGS += -DTS16
else ifeq ($(TS), 32)
  OPTFLAGS += -DTS32
else ifeq ($(TS), 64)
  OPTFLAGS += -DTS64
else ifeq ($(TS), 128)
  OPTFLAGS += -DTS128
else ifeq ($(TS), 256)
  OPTFLAGS += -DTS256
else ifeq ($(TS), 512)
  OPTFLAGS += -DTS512
else
  OPTFLAGS += -DTS32
endif
TARGET := $(TARGET)_$(TS)
endif


#===============================================================================
# Targets to Build
#===============================================================================

#$(LDAPP) $(CXX) $(CXXFLAGS_THREADS) -o $@ $+ $(LDFLAGS) $(CXXFLAGS)
CFLAGS += -I. $(OPTFLAGS)

OBJS = $(OBJ)
TARGETS = $(TARGET)

all: $(TARGETS)

$(TARGET):  $(OBJ)
	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS)

$(OBJ): $(SRC)
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm -rf *~ *.dSYM nc.vg.* $(OBJS) $(TARGETS)

run:
	./$(TARGET) -n 1000
