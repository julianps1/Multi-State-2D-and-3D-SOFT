# top‑level Makefile

# compiler + flags
CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -Iinclude -DQM2D_USE_LAPACKE
LDFLAGS := -lfftw3 -lm 
LAPACKE_LIBS ?= -llapacke -llapack -lblas
LDFLAGS += $(LAPACKE_LIBS)

# sources and derived object names
#SRCS   := src/main.cpp src/grid.cpp src/wave.cpp src/globals.cpp
SRCS   := $(wildcard src/*.cpp)
OBJS   := $(SRCS:.cpp=.o)

# executable name
PROG   := sim          # choose whatever you like

all: $(PROG)

# link
$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# compile rule for .cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(PROG)

.PHONY: all clean
