
INCLUDE_DIRS = ext/eigen-3.4.0

CXX = clang++
CXXFLAGS = -std=c++20 -Wall -Wextra $(addprefix -I, $(INCLUDE_DIRS))
BUILD_DIR = build
SRC_DIR = src
TEST_DIR = test

SRC_FILES = $(wildcard $(SRC_DIR)/*.cc)
TEST_FILES = $(wildcard $(TEST_DIR)/*.cc)
OBJ_FILES = $(addprefix $(BUILD_DIR)/,$(notdir $(SRC_FILES:.cc=.o)))
TEST_OBJ_FILES = $(addprefix $(BUILD_DIR)/,$(notdir $(TEST_FILES:.cc=.o)))
LIBRARY = $(BUILD_DIR)/libsolver.a
EXECUTABLE = $(BUILD_DIR)/test_solver

.PHONY: all clean

all: $(LIBRARY) $(EXECUTABLE)

$(LIBRARY): $(OBJ_FILES)
	ar rcs $@ $^

$(EXECUTABLE): $(TEST_OBJ_FILES) $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cc $(SRC_DIR)/%.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cc | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -I$(SRC_DIR) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR)

