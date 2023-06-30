# source files
SRC := $(wildcard $(CURDIR)/src/*.cc)
OBJ := $(SRC:.cc=.o)

BINSRC := $(wildcard $(CURDIR)/src/bin-*.cc)
BINOBJ := $(BINSRC:.cc=.o)

TESTSRC := $(wildcard $(CURDIR)/src/test-*.cc)
TESTOBJ := $(TESTSRC:.cc=.o)

LIBOBJ := $(OBJ)
LIBOBJ := $(filter-out $(BINOBJ),$(LIBOBJ))
LIBOBJ := $(filter-out $(TESTOBJ),$(LIBOBJ))

DEP := $(SRC:.cc=.P)
HDR := $(wildcard src/*.h)

# include directories
INCLUDES += -I$(CURDIR)/include
INCLUDES += -I$(CURDIR)/iB4e
INCLUDES += -I/usr/local/include # For Homebrew

# C++ compiler flags
CXXFLAGS += --std=c++11
CXXFLAGS += -fPIC
CXXFLAGS += -fopenmp
CXXFLAGS += -Wall
CXXFLAGS += -g
CXXFLAGS += -O3

# library paths
LIBS += -L/usr/local/lib # For Homebrew
LIBS += -lgmp -lgmpxx
LIBS += -lCGAL
LIBS += -lm
LIBS += -lboost_filesystem
LIBS += -lboost_program_options
LIBS += -lboost_system
LIBS += -lboost_log

#compile-time variables
VARS += -DPMFE_PATH='"$(CURDIR)"'

BIN = pmfe-findmfe pmfe-scorer pmfe-parametrizer pmfe-subopt pmfe-tests
all: $(OBJ) $(BIN)

-include $(DEP)

debug: CXXFLAGS += -Og
debug: all

pmfe-findmfe: $(LIBOBJ) src/bin-findmfe.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(VARS) $^ -o $@ $(LIBS)

pmfe-scorer: $(LIBOBJ) src/bin-scorer.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(VARS) $^ -o $@ $(LIBS)

pmfe-parametrizer: $(LIBOBJ) src/bin-parametrizer.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(VARS) $^ -o $@ $(LIBS)

pmfe-subopt: $(LIBOBJ) src/bin-subopt.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(VARS) $^ -o $@ $(LIBS)

pmfe-tests: $(LIBOBJ) $(TESTOBJ) src/bin-tests.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(VARS) $^ -o $@ $(LIBS)

%.o: %.cc
	$(CXX) -MD $(CXXFLAGS) $(INCLUDES) $(VARS) -o $@ -c $<
	@cp $*.d $*.P; \
        sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
            -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
        rm -f $*.d

clean:
	-rm -vf $(EXEC) $(OBJ) $(DEP) $(BIN)

install:

uninstall:

.PHONY: clean
