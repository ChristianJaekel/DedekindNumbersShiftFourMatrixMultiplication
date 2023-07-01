CXX = clang++
CXXFLAGS = -std=c++17 -O3 -ltbb -w -march=native -DEIGEN_NO_DEBUG

ROOT_DIR = ..
SRC_DIR = $(ROOT_DIR)/src
CPU_SRC_DIR = $(ROOT_DIR)/CPU/src
BUILD_DIR = build
DATA_DIR = $(ROOT_DIR)/data/
EIGEN_DIR = $(ROOT_DIR)/eigen/Eigen/Dense

CPPFLAGS += -DDATA_DIR=\"$(DATA_DIR)\"
CPPFLAGS += -DEIGEN_DIR=\"$(EIGEN_DIR)\"

SRCS = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(CPU_SRC_DIR)/*.cpp)

OBJ = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))

# create build directory if it doesn't exist
$(shell mkdir -p $(BUILD_DIR))

TARGET = dedekind

$(TARGET): $(OBJ) $(CUOBJ)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(SRC_DIR) -o $@ $(OBJ)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(CPU_SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	rm -rf $(TARGET) $(BUILD_DIR)

.PHONY: clean
