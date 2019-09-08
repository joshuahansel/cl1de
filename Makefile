CXX ?= g++

BUILD_DIR := build

SRC_DIRS := src
SRC_DIRS_TEST := src test/src

SRCS := $(shell find $(SRC_DIRS) -name "*.cpp")
SRCS_TEST := $(shell find $(SRC_DIRS_TEST) -name "*.cpp")

OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
OBJS_TEST := $(SRCS_TEST:%=$(BUILD_DIR)/%.o)

DEPS := $(OBJS_TEST:.o=.d)

INC_DIRS := $(shell find include test/include -type d)

INC_FLAGS := $(addprefix -I,$(INC_DIRS))

LINK_FLAGS := -lgtest -lgtest_main -lpthread

CPPFLAGS := $(INC_FLAGS) -MMD -MP

MKDIR_P := mkdir -p

# executable targets

PROBLEMS := $(shell ls ./problems)

define make-problem-target
$1: $(OBJS) $(BUILD_DIR)/problems/$1/$1.cpp.o
	@$(CXX) $(OBJS) $(BUILD_DIR)/problems/$1/$1.cpp.o -o ./problems/$1/run
all:: $1
endef

$(foreach PROBLEM,$(PROBLEMS),$(eval $(call make-problem-target,$(PROBLEM))))

# test

test: $(OBJS_TEST) $(BUILD_DIR)/test/main.cpp.o
	@$(CXX) $(OBJS_TEST) $(BUILD_DIR)/test/main.cpp.o -o run_tests $(LINK_FLAGS)

# object targets

$(BUILD_DIR)/%.cpp.o: %.cpp
	@echo Compiling $<...
	@$(MKDIR_P) $(dir $@)
	@$(CXX) $(CPPFLAGS) -c $< -o $@

# clean

.PHONY: clean

clean:
	@$(RM) -r $(BUILD_DIR)
	@$(RM) $(foreach PROBLEM,$(PROBLEMS),./problems/$(PROBLEM)/run)
	@$(RM) run_tests

# includes

-include $(DEPS)
