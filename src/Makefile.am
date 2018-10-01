if DEBUG
AM_CFLAGS = -g -O0
AM_CXXFLAGS = -g -O0 -L/usr/local/Cellar/gsl/2.5/include/
else
AM_CFLAGS = -O2
AM_CXXFLAGS = -O2 -L/usr/local/Cellar/gsl/2.5/include/
endif

if PROFILING
AM_CFLAGS += -pg
AM_CXXFLAGS += -pg
endif

bin_PROGRAMS=vtcombiner

vtcombiner_SOURCES=main.cc
