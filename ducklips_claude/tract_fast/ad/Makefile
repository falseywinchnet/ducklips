# The derivative toolchain's native pieces (tract_fast/ad), built into
# ../build_ad as one shared library for analysis/gradients.py:
#  - forced_folds_api.cpp with the reference engine's flags (build_dev:
#    Release -O3, -mcpu=native -ffp-contract=fast), so its double copy of the
#    folds can be compared with sustain bit for bit;
#  - peel_api.cpp with -ffp-contract=off, as numba compiles reverse_core.py.
#   make -C ad

BUILD ?= ../build_ad
CXX = /usr/bin/c++
FOLDS_CONTRACT ?= fast
COMMON = -O3 -DNDEBUG -std=gnu++20 -arch arm64 -Wall -Wextra -Wpedantic -mcpu=native -fPIC

all: $(BUILD)/libtract_ad.dylib

$(BUILD)/forced_folds_api.o: forced_folds_api.cpp forced_folds.hpp dual.hpp
	mkdir -p $(BUILD)
	$(CXX) $(COMMON) -ffp-contract=$(FOLDS_CONTRACT) -c -o $@ forced_folds_api.cpp

$(BUILD)/peel_api.o: peel_api.cpp geometry_peel.hpp dual.hpp
	mkdir -p $(BUILD)
	$(CXX) $(COMMON) -ffp-contract=off -c -o $@ peel_api.cpp

$(BUILD)/libtract_ad.dylib: $(BUILD)/forced_folds_api.o $(BUILD)/peel_api.o
	$(CXX) -dynamiclib -o $@ $^

clean:
	rm -f $(BUILD)/libtract_ad.dylib $(BUILD)/forced_folds_api.o $(BUILD)/peel_api.o $(BUILD)/libforced_folds.dylib

.PHONY: all clean
