# Compiler flags
CXX      ?= g++
CXXFLAGS ?= -std=c++20 -Wall -O3 # -g --coverage -pedantic
CPPFLAGS ?= -I include -I include/core # Include flags

# Linker flags
LDFLAGS ?=
LDLIBS  ?=

# Variables
EXEC    = main
SRC_DIR = src
SRCS 	= $(shell find $(SRC_DIR) -name '*.cpp')
OBJS    = $(SRCS:.cpp=.o)
HEADERS = $(shell find include -maxdepth 1 -name '*.hpp')

# Default target
all: $(EXEC)

# Link object files to create executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $@

# Compile source files
%.o: %.cpp ../include/core/%hpp $(HEADERS)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

# Remove all object files
clean:
	$(RM) $(OBJS)
	$(RM) -r $(SRC_DIR)/*.gcda $(SRC_DIR)/*.gcno test_coverage* callgrind*

# Remove all generated files
distclean: clean
	$(RM) $(EXEC)
	$(RM) *.csv *.out *.bak *~
	$(RM) $(SRC_DIR)/*~


coverage: all
	lcov --directory . --zerocounters
	./$(EXEC)
	lcov --directory . --capture --no-external --output test_coverage.info
	genhtml test_coverage.info --output test_coverage

memcheck: all
	valgrind --tool=memcheck ./$(EXEC)

profile: all
	valgrind --tool=callgrind ./$(EXEC)