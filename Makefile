OR_TOOLS_TOP=../or-tools
LOCALSOLVER_TOP = /opt/localsolver_11_0

CFLAGS := -std=c++17 -isystem $(OR_TOOLS_TOP)/include -isystem $(LOCALSOLVER_TOP)/include

# During development uncomment the next line to have debug checks and other verifications
DEVELOPMENT = true
ifeq ($(DEVELOPMENT), true)
  CFLAGS := $(CFLAGS) -O0 -DDEBUG -ggdb3 -fsanitize=address -fkeep-inline-functions -fno-inline-small-functions
  CXX := LSAN_OPTION=verbosity=1:log_threads=1 $(CXX) # adress sanitizer works only if the executable launched without gdb
else
  CFLAGS := $(CFLAGS) -O3 -DNDEBUG
endif

# Activate warnings
CFLAGS := $(CFLAGS) -Wall -Wextra -Wshadow -Wmissing-include-dirs -Wswitch-enum -Wfloat-equal -Wundef

# The following is to supress a warning due to a protobuf that is fixed at v3.14.
# It can be removed when or-tools is upgraded to v8.1+ (where the protobuf dependency is upgraded to v3.14).
PROTOBUF_VERSION := $(shell $(OR_TOOLS_TOP)/bin/protoc --version | cut -d" " -f2)
ifeq ($(shell dpkg --compare-versions $(PROTOBUF_VERSION) 'lt' '3.14' && echo true), true)
	CFLAGS := $(CFLAGS) -Wno-array-bounds
endif

.PHONY: all local_clean

all: $(EXE)

%.pb.cc: %.proto
	$(OR_TOOLS_TOP)/bin/protoc  --cpp_out . $<

%.o: %.cc %.h
	$(CXX) $(CFLAGS) -c $< -o $@

localsolver_vrp.pb.h: localsolver_vrp.pb.cc

localsolver_result.pb.h: localsolver_result.pb.cc

tsp_localsolver.o: tsp_localsolver.cc \
	localsolver_vrp.pb.h \
	localsolver_result.pb.h
	$(CXX) $(CFLAGS) -c tsp_localsolver.cc -o tsp_localsolver.o

tsp_localsolver: tsp_localsolver.o localsolver_vrp.pb.o localsolver_result.pb.o
	$(CXX) $(CFLAGS) -g tsp_localsolver.o localsolver_vrp.pb.o localsolver_result.pb.o -lz -lrt -lpthread \
	-L$(LOCALSOLVER_TOP)/bin -llocalsolver110 \
	-L $(OR_TOOLS_TOP)/lib -Wl,-rpath $(OR_TOOLS_TOP)/lib -lortools -lprotobuf -lglog -lgflags -labsl_raw_hash_set -labsl_time -labsl_time_zone \
	-o tsp_localsolver


clean:
	rm -f *.pb.cc *.pb.h *.o

mrproper: clean
	rm tsp_localsolver
