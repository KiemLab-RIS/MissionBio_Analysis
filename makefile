# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -g

# Source files
SRC = main.cpp NWAlign.cpp ReadFastq.cpp Params.cpp Output.cpp ReadProtein.cpp Whitelist.cpp

# Object files (replace .cpp with .o)
OBJ = $(SRC:.cpp=.o)

# Output executable
TARGET = mbFilter

# Default target
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

# Rule to compile each .cpp file into a .o file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJ) $(TARGET)
