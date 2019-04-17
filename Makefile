CXX ?= g++

BUILD_DIR ?= ./build
SRC_DIRS ?= ./src

SRCS := $(shell find $(SRC_DIRS) -name "*.cpp")
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find ./include -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

MKDIR_P := mkdir -p

# executable targets

PROBLEMS := $(shell ls ./problems)

define make-problem-target
$1: $(OBJS) $(BUILD_DIR)/problems/$1/$1.cpp.o
	@$(CXX) $(OBJS) $(BUILD_DIR)/problems/$1/$1.cpp.o -o ./problems/$1/run
all:: $1
endef

$(foreach PROBLEM,$(PROBLEMS),$(eval $(call make-problem-target,$(PROBLEM))))

# object targets

$(BUILD_DIR)/%.cpp.o: %.cpp
	@$(MKDIR_P) $(dir $@)
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# clean

.PHONY: clean

clean:
	@$(RM) -r $(BUILD_DIR)
	@$(RM) $(foreach PROBLEM,$(PROBLEMS),./problems/$(PROBLEM)/run)

# includes

-include $(DEPS)
