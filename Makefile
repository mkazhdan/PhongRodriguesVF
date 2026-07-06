COMPILER ?= gcc
#COMPILER ?= clang

programs = Bracket SparseInterpolation

# Allow "make -j" to operate in parallel over the programs.
all: $(programs)
$(programs):
	$(MAKE) -C $@ COMPILER=$(COMPILER)

programs_debug = $(foreach n,$(programs),debug_$(n))  # pseudo-dependency to allow "make -j" parallelism
debug: $(programs_debug)
$(programs_debug):
	$(MAKE) -C $(@:debug_%:%) debug COMPILER=$(COMPILER)

programs_clean = $(foreach n,$(programs),clean_$(n))  # pseudo-dependency to allow "make -j" parallelism
clean: $(programs_clean)
$(programs_clean):
	$(MAKE) -C $(@:clean_%=%) clean

.PHONY: $(programs) $(programs_debug) $(programs_clean)
