
OBJDIR = obj
SRCDIR = src
INCDIR = include
BINDIR = bin

SOURCES = $(basename $(wildcard $(SRCDIR)/*.cpp))
OBJECTS = $(patsubst $(SRCDIR)/%,$(OBJDIR)/%.o,$(SOURCES))
INCLUDES = $(wildcard $(SRCDIR)/*.hpp)

CUDA_SOURCES = $(basename $(wildcard $(SRCDIR)/*.cu))
CUDA_OBJECTS = $(patsubst $(SRCDIR)/%,$(OBJDIR)/%.o,$(CUDA_SOURCES))
CUDA_INCLUDES = $(wildcard $(SRCDIR)/*_cuda.h)

EXEC = ${BINDIR}/sMC
LIBS =
INC = -I${INCDIR}

GCCDIR = /home/dg6/code/gcc/gcc-4.4.5
ALGLIBDIR = /home/dg6/code/alglib/cpp/src
NETCDFDIR = /home/dg6/code/netcdf/netcdf_gnu64
BOOSTDIR = /usr/include
CUDADIR = /home/dg6/code/cuda/4.0/cuda
CUDALIBDIR = ${CUDADIR}/lib64
CUDA_SDK_DIR = /home/dg6/code/cuda/NVIDIA_GPU_Computing_SDK/C/src/simplePrintf
CUDA_ARCH = sm_13

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

ALGLIB = $(wildcard $(ALGLIBDIR)/*.o)
OBJECTS += ${ALGLIB} 
INC += -I${ALGLIBDIR}

NETCDF = -L${NETCDFDIR}/lib -lnetcdf_c++ -lnetcdf
LIBS += ${NETCDF}
INC += -I${NETCDFDIR}/include

INC += -I${BOOSTDIR}

CXX = ${GCCDIR}/bin/g++
CXXFLAGS = -Wall -g -pg
LDFLAGS = -pg

NVCC = ${CUDADIR}/bin/nvcc
NVCCFLAGS = -g -G --compiler-bindir ${GCCDIR}/bin -arch ${CUDA_ARCH}
INC += -I${CUDADIR}/include -I${CUDA_SDK_DIR}
LIBS += -L${CUDALIBDIR} -lcuda -lcudart

${EXEC}: ${OBJECTS} ${CUDA_OBJECTS}
	${CXX} ${LDFLAGS} ${OBJECTS} ${CUDA_OBJECTS} ${LIBS} -o $@

${OBJDIR}/%.o: ${SRCDIR}/%.cpp ${INCLUDES}
	${CXX} -c ${INC} ${CXXFLAGS} $< -o $@

${OBJDIR}/%.o: ${SRCDIR}/%.cu ${CUDA_INCLUDES}
	${NVCC} -c ${INC} ${NVCCFLAGS} $< -o $@

.PHONY: clean

clean:
	rm -f ${OBJDIR}/*.o ${EXEC}
