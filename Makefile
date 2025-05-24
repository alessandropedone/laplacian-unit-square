CXX      = mpic++
CXXFLAGS = -std=c++20 -Wall -O3 -MMD -MP
CPPFLAGS = -I include -I include/core
# Disable OpenMP if OPENMP=0 is specified
OPENMP ?= 1
ifeq ($(OPENMP),1)
	CXXFLAGS += -fopenmp
endif


OPENMP_FLAG_FILE := .openmp_flag

# Check if OPENMP value changed
ifeq ($(shell test -f $(OPENMP_FLAG_FILE) && grep -qx 'OPENMP=$(OPENMP)' $(OPENMP_FLAG_FILE) && echo same),same)
else
.PHONY: force_recompile
all: force_recompile
force_recompile: distclean
endif

# Always update the flag file after build
all: update_openmp_flag

.PHONY: update_openmp_flag
update_openmp_flag:
	@echo "OPENMP=$(OPENMP)" > $(OPENMP_FLAG_FILE)

EXEC    = main
SRC_DIR = src
SRCS    = $(shell find $(SRC_DIR) -name '*.cpp')
OBJS    = $(SRCS:.cpp=.o)
DEPS    = $(OBJS:.o=.d)

.PHONY: run
run: $(EXEC)
	mpirun -np 2 ./$(EXEC)

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
