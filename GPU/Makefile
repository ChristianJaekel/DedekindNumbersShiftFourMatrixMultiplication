CXX = nvcc
CXXFLAGS = -std=c++17 -O3 -ltbb -w -arch=native -lcublas

ROOT_DIR = ..
SRC_DIR = $(ROOT_DIR)/src
GPU_SRC_DIR = $(ROOT_DIR)/GPU/src
BUILD_DIR = build
DATA_DIR = $(ROOT_DIR)/data/

CPPFLAGS += -DDATA_DIR=\"$(DATA_DIR)\"

SRCS = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(GPU_SRC_DIR)/*.cpp)
CUSRCS = $(wildcard $(GPU_SRC_DIR)/*.cu)

OBJ = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))
CUOBJ = $(patsubst $(GPU_SRC_DIR)/%.cu, $(BUILD_DIR)/%.o, $(CUSRCS))

# create build directory if it doesn't exist
$(shell mkdir -p $(BUILD_DIR))

TARGET = dedekind

$(TARGET): $(OBJ) $(CUOBJ)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(SRC_DIR) -o $@ $(OBJ) $(CUOBJ)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(GPU_SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(GPU_SRC_DIR)/%.cu
	$(CXX) $(CXXFLAGS) -I$(SRC_DIR) -c $< -o $@

clean:
	rm -rf $(TARGET) $(BUILD_DIR)

.PHONY: clean
