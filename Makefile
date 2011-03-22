
OBJDIR = obj
SRCDIR = src
INCDIR = include
BINDIR = bin

GCCVER = /home/dg6/code/gcc/gcc-4.4.5

CXX = ${GCCVER}/bin/g++
CXXFLAGS = -Wall -g
LDFLAGS = 
INC = -I${INCDIR}

SOURCES = $(basename $(wildcard $(SRCDIR)/*.cpp))
OBJECTS = $(patsubst $(SRCDIR)/%,$(OBJDIR)/%.o,$(SOURCES))
EXEC = ${BINDIR}/sMC

${EXEC}: ${OBJECTS}
	${CXX} ${LDFLAGS} ${OBJECTS} -o $@

${OBJDIR}/%.o: ${SRCDIR}/%.cpp
	${CXX} -c ${INC} ${CXXFLAGS} $< -o $@

.PHONY: clean

clean:
	echo ${SOURCES}
	echo ${OBJECTS}
	rm -f ${OBJDIR}/*.o ${EXEC}
