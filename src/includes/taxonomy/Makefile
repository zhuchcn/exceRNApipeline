CXX = g++
LD = g++

STDVERSION = -std=c++14
WARNINGS = -Wall -Wextra

CXXFLAG = $(STDVERSION) $(WARNINGS) -g -O0
LDFLAG = $(STDVERSION) $(WARNINGS) -g -O0

OBJS_DIR = .objs
CPP_TEST = $(wildcard tests/*.cpp)
OBJS_TEST = NCBITaxonomy.o TaxaCounter.o
OBJS_TEST += $(CPP_TEST:.cpp=.o)
OBJS = NCBITaxonomy.o main.o
OBJS_LIB = ../gzstream/

all: main

main: $(patsubst %.o, $(OBJS_DIR)/%.o, $(OBJS)) $(OBJS_DIR)/gzstream/gzstream.o
	$(LD) $^ $(LDFLAG) -o $@ -lz

test: $(patsubst %.o, $(OBJS_DIR)/%.o, $(OBJS_TEST)) $(OBJS_DIR)/gzstream/gzstream.o
	$(LD) $^ $(LDFLAG) -o $@ -lz

$(OBJS_DIR)/NCBITaxonomy.o: NCBITaxonomy.cpp | $(OBJS_DIR)
	$(CXX) $(CXXFLAG) -O0 -g -c $< -o $@

$(OBJS_DIR)/TaxaCounter.o: TaxaCounter.cpp | $(OBJS_DIR)
	$(CXX) $(CXXFLAG) -O0 -g -c $< -o $@

$(OBJS_DIR)/tests/test_ncbi_taxonomy.o: tests/test_ncbi_taxonomy.cpp | $(OBJS_DIR)
	$(CXX) $(CXXFLAG) -O0 -g -c $< -o $@

$(OBJS_DIR)/gzstream/gzstream.o: ../gzstream/gzstream.C | $(OBJS_DIR)
	$(CXX) $(CXXFLAG) -I ../gzstream -O -g -c -o $@ $< -lz

$(OBJS_DIR):
	@mkdir -p $(OBJS_DIR)
	@mkdir -p $(OBJS_DIR)/tests
	@mkdir -p $(OBJS_DIR)/gzstream

clean:
	rm -rf main test $(OBJS_DIR) *.o *.d

tidy: clean
	rm -rf doc

.PHONY: all tidy clean