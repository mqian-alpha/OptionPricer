# Makefile for OptionPricer project

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -g -Wall -std=c++11

# Target executable
TARGET = OptionPricer

# Source files
SRCS = test_main.cpp BSOptionPrice.cpp

# Object files (replace .cpp with .o)
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)