OBJECTS = \
    readpoly.o \
    triangulation_search_structure.o \
    point_triangle_distance.o \
    templ_instant_distance_map.o \

LIBRARIES = \
    -lCGAL \
    -lboost_program_options \

CXXFLAGS = \
    -frounding-math \
    -fno-strict-aliasing \
    -Wall \

CXXFLAGS += \
    -O3 -DNDEBUG

BINARIES = \
    voxelsurface

all: $(BINARIES)

voxelsurface: voxelsurface.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ voxelsurface.o $(OBJECTS) $(LIBRARIES)

ts.headers: *.h Makefile
	$(MAKE) clean
	touch $@

%.o: %.cpp ts.headers
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f *.o ts.headers $(BINARIES)

.PHONY: all clean
