# Set PACS_ROOT to the path where pacs-examples/Examples is located
PACS_ROOT = ../pacs-examples/Examples

# compiler and compilation options
CXX ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -I. -Wno-conversion-null -Wno-deprecated-declarations -I$(PACS_ROOT)/include 

LDFLAGS ?= -L$(PACS_ROOT)/lib -Wl,-rpath=$(PACS_ROOT)/lib

# executable name 
EXEC = main

# source and object files
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

# header files
HDRS = $(wildcard *.hpp)

.PHONY = all $(EXEC) $(OBJS) clean distclean

# rules
all: $(EXEC)

%.o: %.cpp $(HDRS)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
 
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) *.o

distclean: clean
	$(RM) $(EXEC)

