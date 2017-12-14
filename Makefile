R_HOME :=               $(shell R RHOME)

TARGET_EXEC ?= SE_Compare

BUILD_DIR ?= ./build
SRC_DIRS ?= ./src

SRCS := $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

## include headers and libraries for R 
RCPPFLAGS :=            $(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS :=             $(shell $(R_HOME)/bin/R CMD config --ldflags)

## include headers and libraries for Rcpp interface classes
RCPPINCL :=             $(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS :=             $(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## include headers and libraries for RInside embedding classes
RINSIDEINCL :=          $(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R $(R_ARCH) --vanilla --slave)
RINSIDELIBS :=          $(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R $(R_ARCH) --vanilla --slave)


$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS) $(RLDFLAGS) $(RCPPLIBS) $(RINSIDELIBS)

# assembly
$(BUILD_DIR)/%.s.o: %.s
	$(MKDIR_P) $(dir $@)
	$(AS) $(ASFLAGS) -c $< -o $@

# c source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) -c $< -o $@


.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p

