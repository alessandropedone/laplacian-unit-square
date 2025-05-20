CXX      = mpic++
CXXFLAGS = -std=c++20 -Wall -O3 -MMD -MP -fopenmp
CPPFLAGS = -I include -I include/core

EXEC    = main
SRC_DIR = src
SRCS    = $(shell find $(SRC_DIR) -name '*.cpp')
OBJS    = $(SRCS:.cpp=.o)
DEPS    = $(OBJS:.o=.d)

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

clean:
	$(RM) $(OBJS) $(DEPS)
	$(RM) -r $(SRC_DIR)/*.gcda $(SRC_DIR)/*.gcno test_coverage* callgrind*

distclean: clean
	$(RM) $(EXEC)
	$(RM) *.csv *.out *.bak *~
	$(RM) $(SRC_DIR)/*~

-include $(DEPS)
