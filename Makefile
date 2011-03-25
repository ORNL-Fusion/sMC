
OBJDIR = obj
SRCDIR = src
INCDIR = include
BINDIR = bin

SOURCES = $(basename $(wildcard $(SRCDIR)/*.cpp))
OBJECTS = $(patsubst $(SRCDIR)/%,$(OBJDIR)/%.o,$(SOURCES))
INCLUDES = $(wildcard $(SRCDIR)/*.hpp)
EXEC = ${BINDIR}/sMC
LIBS =
INC = -I${INCDIR}

GCCDIR = /home/dg6/code/gcc/gcc-4.4.5
ALGLIBDIR = /home/dg6/code/alglib/cpp/src
NETCDFDIR = /home/dg6/code/netcdf/netcdf_gnu64
BOOSTDIR = /usr/include

# Catch for greendl (my laptop)

ifeq ($(findstring greendl,$(HOSTNAME_OSX)),greendl)
	GCCDIR = /opt/local
	ALGLIBDIR = /home/dg6/code/alglib/cpp/src
	NETCDFDIR = /opt/local
endif

ALGLIB = $(wildcard $(ALGLIBDIR)/*.o)
OBJECTS += ${ALGLIB} 
INC += -I${ALGLIBDIR}

NETCDF = -L${NETCDFDIR}/lib -lnetcdf_c++ -lnetcdf
LIBS += ${NETCDF}
INC += -I${NETCDFDIR}/include

INC += -I${BOOSTDIR}

CXX = ${GCCDIR}/bin/g++
CXXFLAGS = -Wall -g
LDFLAGS =

${EXEC}: ${OBJECTS}
	${CXX} ${LDFLAGS} ${OBJECTS} ${LIBS} -o $@

${OBJDIR}/%.o: ${SRCDIR}/%.cpp ${INCLUDES}
	${CXX} -c ${INC} ${CXXFLAGS} $< -o $@

.PHONY: clean

clean:
	rm -f ${OBJDIR}/*.o ${EXEC}
