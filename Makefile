
OBJDIR = obj
SRCDIR = src
INCDIR = include
BINDIR = bin

GCCDIR = /home/dg6/code/gcc/gcc-4.4.5

# Catch for greendl (my laptop)

ifeq ($(findstring greendl,$(HOSTNAME_OSX)),greendl)
	GCCDIR = /opt/local
endif

CXX = ${GCCDIR}/bin/g++
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
	rm -f ${OBJDIR}/*.o ${EXEC}
