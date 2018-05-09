################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm

# All of the sources participating in the build are defined here
-include sources.mk
-include src/tomahawk/two/subdir.mk
-include src/tomahawk/subdir.mk
-include src/third_party/zlib/subdir.mk
-include src/third_party/xxhash/subdir.mk
-include src/support/subdir.mk
-include src/math/subdir.mk
-include src/io/vcf/subdir.mk
-include src/io/compression/subdir.mk
-include src/io/bcf/subdir.mk
-include src/io/subdir.mk
-include src/index/subdir.mk
-include src/algorithm/sort/subdir.mk
-include src/algorithm/subdir.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 
+BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
+ifneq ($BRANCH, master)
+GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)-$(BRANCH)
+else
+GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)
+endif 

# All Target
all: tomahawk

# Tool invocations
tomahawk: $(OBJS) $(USER_OBJS)
	g++ -pthread -o "tomahawk" $(OBJS) $(USER_OBJS) $(LIBS)
	$(MAKE) --no-print-directory post-build

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(OBJS)$(C_UPPER_DEPS)$(CXX_DEPS)$(C_DEPS)$(CPP_DEPS) tomahawk

post-build:
	-mkdir -p ../bin; mv tomahawk ../bin/;

.PHONY: all clean dependents post-build

-include ../makefile.targets
