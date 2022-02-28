# * URL: http://www.partow.net/programming/makefile/index.html *
# make all
# make clean
# make program
# make build
# make release
# make debug

GIT_COMMIT  := $(shell git rev-parse --verify HEAD)
GIT_DATE    := $(firstword $(shell git --no-pager show --date=iso-strict --format="%ad" --name-only))
BUILD_DATE  := $(shell date --iso=seconds)

# -pthread CXXFLAGS
# -ggdb 
CXX      := -g++
CXXFLAGS := -std=c++11 -pedantic-errors -Wall -Wextra -Werror -pthread \
			-DGIT_COMMIT=\"'$(GIT_COMMIT)'\"\
 			-DGIT_DATE=\"'$(GIT_DATE)'\"\
 			-DBUILD_DATE=\"'$(BUILD_DATE)'\"
LDFLAGS  := -L/usr/lib
LDLIBS   := -lm -lstdc++fs -lboost_system -lboost_filesystem -lboost_chrono
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/bin
TARGET   := malaria_model
INCLUDE  := -Iinclude/malaria_model/ -Iinclude/third_party/
SRC      :=                      \
	$(wildcard src/util/*.cpp) \
	$(wildcard src/common/*.cpp) \
	$(wildcard src/human/*.cpp) \
	$(wildcard src/village/*.cpp) \
	$(wildcard src/simulation/*.cpp) \
	$(wildcard src/intervention/*.cpp) \
	$(wildcard src/*.cpp)         \

OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

# $(info $$OBJECTS is [${OBJECTS}])

TARGET_TEST  := test
INCLUDE_TEST := $(INCLUDE) -Itest/third_party/catch2
SRC_TEST 	 := $(wildcard test/*.cpp)
OBJECTS_TEST := $(filter-out $(OBJ_DIR)/src/$(TARGET).o, $(OBJECTS)) $(SRC_TEST:test/%.cpp=$(OBJ_DIR)/test/%.o)



NVCC         := nvcc
NVCC_FLAGS   := -std=c++11
NVCC_LDLIBS  := $(LDLIBS) -lcurand #-lprofiler
TARGET_CUDA  := $(TARGET)_cuda
SRC_CUDA	 := \
	$(wildcard src/util/*.cu) \
	$(wildcard src/village/*.cu)

OBJECTS_CUDA := $(SRC_CUDA:%.cu=$(OBJ_DIR)/%_cuda.o)

# $(info $$OBJECTS_CUDA is [${OBJECTS_CUDA}])

USE_NVCC = 0

.DELETE_ON_ERROR:

.PHONY: all
all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) -o $(APP_DIR)/$(TARGET) $(OBJECTS) $(LDLIBS)


.PHONY: test

test: build $(APP_DIR)/$(TARGET_TEST)

$(OBJ_DIR)/test/%.o: test/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET_TEST): $(OBJECTS_TEST)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE_TEST) $(LDFLAGS) -o $(APP_DIR)/$(TARGET_TEST) $(OBJECTS_TEST) $(LDLIBS)

test-cuda: CXX := nvcc
test-cuda: CXXFLAGS := -Xcompiler="$(CXXFLAGS)"
test-cuda: LDFLAGS := $(NVCC_LDLIBS)
test-cuda: OBJECTS_TEST += $(OBJECTS_CUDA)
test-cuda: cuda test

.PHONY: cuda

cuda: CXX := nvcc
cuda: CXXFLAGS := -Xcompiler="$(CXXFLAGS)"
cuda: LDFLAGS := $(NVCC_LDLIBS)
cuda: build $(APP_DIR)/$(TARGET_CUDA)

release-cuda: CXX := nvcc
release-cuda: TARGET_CUDA := $(TARGET)_release_cuda
release-cuda: CXXFLAGS := -Xcompiler="$(CXXFLAGS) -O3"
release-cuda: LDFLAGS := $(NVCC_LDLIBS)
release-cuda: build $(APP_DIR)/$(TARGET_CUDA)

$(OBJ_DIR)/%_cuda.o: %.cu
	@mkdir -p $(@D)
	$(NVCC) $(NVCC_FLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET_CUDA): $(OBJECTS) $(OBJECTS_CUDA)
	@mkdir -p $(@D)
	$(NVCC) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) -o $(APP_DIR)/$(TARGET_CUDA) $(OBJECTS) $(OBJECTS_CUDA) $(NVCC_LDLIBS)





.PHONY: build clean
build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

# build-cuda: build
# 	USE_NVCC := y

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O3
release: all

clean:
	@rm -r $(OBJ_DIR)/*
	@rm -r $(APP_DIR)/*	
