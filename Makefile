# Define variables
CXX = g++
CXXFLAGS = -Iinclude -std=c++17
LDFLAGS = -lfl

# Directories
SRC_DIR = src
OBJ_DIR = obj
INCLUDE_DIR = include

# Source files and object files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Name of the executable
EXECUTABLE = physwiz

# Default target to build the executable
all: $(EXECUTABLE)

# Link object files to create the executable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

# Compile source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/lexer.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to generate lexer.cpp from lexer.l
$(SRC_DIR)/lexer.cpp: $(SRC_DIR)/lexer.l
	flex -o $(SRC_DIR)/lexer.cpp $(SRC_DIR)/lexer.l

# Clean target to remove build artifacts
clean:
	rm -rf $(OBJ_DIR) $(EXECUTABLE)

# Phony target to ensure 'clean' is not considered a file
.PHONY: clean
