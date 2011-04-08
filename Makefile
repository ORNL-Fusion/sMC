NAME := bin/sMC

GCCDIR := /home/dg6/code/gcc/gcc-4.4.5/bin
ALGLIBDIR := /home/dg6/code/alglib/cpp/src
NETCDFDIR := /home/dg6/code/netcdf/netcdf_gnu64
CUDADIR := /home/dg6/code/cuda/4.0/cuda
CUDALIBDIR = ${CUDADIR}/lib64
CUDA_ARCH := sm_13
CUDA_SDK_DIR := /home/dg6/code/cuda/NVIDIA_GPU_Computing_SDK/C/src/simplePrintf

# Catch for greendl (my laptop)

ifeq ($(findstring greendl,$(HOSTNAME_OSX)),greendl)
	GCCDIR = /opt/local
	ALGLIBDIR = /home/dg6/code/alglib/cpp/src
	NETCDFDIR = /opt/local
	BOOSTDIR = /opt/local/lib/boost
	CUDADIR = /usr/local/cuda
	CUDALIBDIR = ${CUDADIR}/lib
	CUDA_SDK_DIR = /Developer/GPU\ Computing/C/src/simplePrintf
	CUDA_ARCH = sm_11
endif

CC := $(GCCDIR)/gcc
CPP := $(CUDADIR)/bin/nvcc #$(GCCDIR)/g++
NVCC := $(CUDADIR)/bin/nvcc
LINK := $(CPP)

MODULES := src include

INCLUDEFLAGS := -I$(ALGLIBDIR) -I$(CUDA_SDK_DIR)
CFLAGS := -g -std=c99
CPPFLAGS :=
NVCCFLAGS := -g -G --compiler-bindir $(GCCDIR) -arch $(CUDA_ARCH)
LFLAGS := -g -L$(NETCDFDIR) -L$(CUDALIBDIR)
LIBS := -lnetcdf_c++ -lnetcdf $(ALGLIBDIR)/*.o -lcuda -lcudart

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
CPPFLAGS += $(INCLUDEFLAGS)
NVCCFLAGS += $(INCLUDEFLAGS)

# determine the object files
SRCTYPES = c cpp cu
OBJ := $(foreach srctype, $(SRCTYPES), $(patsubst %.$(srctype), obj/%.o, $(wildcard $(patsubst %, %/*.$(srctype), $(MODULES)))))

# link the program
$(NAME) : $(OBJ)
	$(LINK) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

# calculate include dependencies
.dep/%.d : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CPPFLAGS), $<) > $@

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
	-@rm $(NAME) $(OBJ) $(DEP)
