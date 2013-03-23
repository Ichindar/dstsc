
# Which compilier to use and compilier options to use
CC= g++
CFLAGS= -O3

# Where are the source files locatede
VPATH= ../SAIS-SK/src:src

# Default rule
default: test

# Define source files
SOURCE_SA = sk-sain.c gt-alloc.c mmap.c fileopen.c
SOURCE_LZ = main.cpp lzOG.cpp

# Define object files
OBJECTS_SA = $(patsubst %.c,obj/%.o,$(SOURCE_SA))
OBJECTS_LZ = $(patsubst %.cpp,obj/%.o,$(SOURCE_LZ))

# Create objects directory if necessary
$(OBJECTS_SA): | obj

obj:
	@mkdir -p $@

# Build objects
obj/%.o : %.cpp
	@echo $<
	$(CC) $(CFLAGS) -c $< -o $@

obj/%.o : %.c
	@echo $<
	$(CC) $(CFLAGS) -c $< -o $@

# Link objects into executable
dir:
	mkdir -p bin

test: $(OBJECTS_SA) $(OBJECTS_LZ) dir
	$(CC) $(CFLAGS) $(OBJECTS_LZ) $(OBJECTS_SA) -o bin/lz_factor
