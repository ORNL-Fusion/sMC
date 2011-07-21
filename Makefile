NAME := bin/sMC

# Defaults for dlg-hp.ornl.gov

GCCDIR := /home/dg6/code/gcc/gcc-4.4.5/bin
ALGLIBDIR := /home/dg6/code/alglib/cpp/src
NETCDFDIR := /home/dg6/code/netcdf/netcdf_gnu64
CUDADIR := /home/dg6/code/cuda/4.0/cuda
CUDALIBDIR = ${CUDADIR}/lib64
CUDA_ARCH := sm_13
CUDA_SDK_DIR := /home/dg6/code/cuda/NVIDIA_GPU_Computing_SDK

# Catch for greendl.* (my laptop)

ifeq ($(findstring greendl,$(HOSTNAME_OSX)),greendl)
GCCDIR := /opt/local/bin
ALGLIBDIR := /home/dg6/code/alglib/cpp/src
NETCDFDIR := /opt/local
CUDADIR := /usr/local/cuda
CUDALIBDIR := ${CUDADIR}/lib
CUDA_ARCH := sm_11
CUDA_SDK_DIR := /Developer/GPU\ Computing
endif

CUDA_SDK_INC := $(CUDA_SDK_DIR)/C/common/inc

CC := $(GCCDIR)/gcc
CPP := $(GCCDIR)/g++
NVCC := $(CUDADIR)/bin/nvcc

MODULES := src include

INCLUDEFLAGS := -I$(ALGLIBDIR) -I$(CUDA_SDK_DIR) -I$(NETCDFDIR)/include -I$(CUDA_SDK_INC)
CFLAGS := 
CPPFLAGS := -g -pg 
NVCCFLAGS := --compiler-bindir $(GCCDIR) -arch $(CUDA_ARCH) --ptxas-options=-v #-g -G 
LFLAGS := -L$(NETCDFDIR)/lib -L$(CUDALIBDIR) 
LIBS := $(ALGLIBDIR)/*.o -lcuda -lcudart -lnetcdf_c++ -lnetcdf

USECUDA:=0
DEBUG:=1
PITCH_SCATTERING:=1
ENERGY_SCATTERING:=1
GOOSE:=1

LINK := $(CPP) -pg

# You shouldn't have to go below here
#
# DLG: 	Added the -x c to force c file type so that 
# 		the .cu files will work too :)

DIRNAME = `dirname $1`
MAKEDEPS = $(CC) -MM -MG $2 -x c $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"

.PHONY : all

all : $(NAME)

# look for include files in each of the modules
INCLUDEFLAGS += $(patsubst %, -I%, $(MODULES))

CFLAGS += $(INCLUDEFLAGS)
CPPFLAGS += $(INCLUDEFLAGS) -DDEBUGLEVEL=$(DEBUG) \
			-DPITCH_SCATTERING=$(PITCH_SCATTERING) \
			-DENERGY_SCATTERING=$(ENERGY_SCATTERING) \
			-DGOOSE=$(GOOSE)
NVCCFLAGS += $(INCLUDEFLAGS) -DDEBUGLEVEL=$(DEBUG) \
			 -DPITCH_SCATTERING=$(PITCH_SCATTERING) \
			 -DENERGY_SCATTERING=$(ENERGY_SCATTERING) \
			 -DGOOSE=$(GOOSE)

# determine the object files
SRCTYPES := c cpp 
ifeq ($(USECUDA),1)
SRCTYPES += cu
CPPFLAGS += -DUSECUDA  #-D__PROFILING__
NVCCFLAGS += #-D__PROFILING__
endif
OBJ := $(foreach srctype, $(SRCTYPES), $(patsubst %.$(srctype), obj/%.o, $(wildcard $(patsubst %, %/*.$(srctype), $(MODULES)))))

# link the program
$(NAME) : $(OBJ)
	$(LINK) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

# calculate include dependencies
.dep/%.d : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(INCLUDEFLAGS), $<) > $@

obj/%.o : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CPP) $(CPPFLAGS) -c -o $@ $<

.dep/%.d : %.c
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CFLAGS), $<) > $@

obj/%.o : %.c
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CC) $(CFLAGS) -c -o $@ $<

.dep/%.d : %.cu
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(INCLUDEFLAGS), $<) > $@

obj/%.o : %.cu
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<


# include the C include dependencies
DEP := $(patsubst obj/%.o, .dep/%.d, $(OBJ))

ifneq ($(MAKECMDGOALS),clean)
include $(DEP)
endif

clean :
	-@rm $(NAME) $(OBJ) $(DEP) .dep/src/*
