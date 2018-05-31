# ################################################################
# Copyright (C) 2016-present Genome Research Ltd.
# Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
# ################################################################

SUBDIRS := \
src/algorithm \
src/algorithm/sort \
src/index \
src/io \
src/io/bcf \
src/io/compression \
src/io/vcf \
src \
src/math \
src/support \
src/third_party/xxhash \
src/third_party/zlib \
src/tomahawk \
src/tomahawk/two \

# Version numbers slices from the source header
LIBVER_MAJOR_SCRIPT:=`sed -n '/const int PROGRAM_VERSION_MAJOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < src/support/MagicConstants.h`
LIBVER_MINOR_SCRIPT:=`sed -n '/const int PROGRAM_VERSION_MINOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < src/support/MagicConstants.h`
LIBVER_PATCH_SCRIPT:=`sed -n '/const int PROGRAM_VERSION_PATCH = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < src/support/MagicConstants.h`
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_MAJOR := $(shell echo $(LIBVER_MAJOR_SCRIPT))
LIBVER_MINOR := $(shell echo $(LIBVER_MINOR_SCRIPT))
LIBVER_PATCH := $(shell echo $(LIBVER_PATCH_SCRIPT))
LIBVER := $(shell echo $(LIBVER_SCRIPT))

.PHONY: all clean cleanmost library dependents

# If you want to build in debug mode then add DEBUG=true to your build command
# make DEBUG=true
ifdef DEBUG
DEBUG_FLAGS := -g -Wall -Wextra -Wcast-qual -Wcast-align -Wshadow \
                  -Wstrict-aliasing=1 -Wswitch-enum -Wdeclaration-after-statement \
                  -Wstrict-prototypes -Wundef -Wpointer-arith -Wformat-security \
                  -Wvla -Wformat=2 -Winit-self -Wfloat-equal -Wwrite-strings \
                  -Wredundant-decls
else
DEBUG_FLAGS := 
endif

# Global build parameters
INCLUDE_PATH := -I"src/"
OPTFLAGS := -O3 -msse4.2
# Legacy flags used
#OPTFLAGS := -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops

ifdef library
OPTFLAGS += -fPIC
endif

# OS X linker doesn't support -soname, and use different extension
# see : https://developer.apple.com/library/mac/documentation/DeveloperTools/Conceptual/DynamicLibraries/100-Articles/DynamicLibraryDesignGuidelines.html
ifneq ($(shell uname), Darwin)
SHARED_EXT = so
LD_LIB_FLAGS := -shared -Wl,-rpath,"./",-soname,ltomahawk.$(SHARED_EXT)
else
SHARED_EXT = dylib
LD_LIB_FLAGS := -dynamiclib -install_name ltomahawk.$(SHARED_EXT)
endif

CXXFLAGS := -std=c++0x $(OPTFLAGS) $(DEBUG_FLAGS)
CFLAGS := -std=c99 $(OPTFLAGS) $(DEBUG_FLAGS)

# Inject git information
BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
ifneq ($(BRANCH), master)
GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)-$(BRANCH)
else
GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)
endif

# Inject subdirs
-include src/tomahawk/two/subdir.mk
-include src/tomahawk/subdir.mk
-include src/third_party/subdir.mk
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

# All Target
all: tomahawk

tomahawk: $(OBJS) $(USER_OBJS)
	g++ -pthread -o "tomahawk" $(OBJS) $(USER_OBJS) $(LIBS)
	$(MAKE) cleanmost
	$(MAKE) library library=true

library: $(OBJS) $(USER_OBJS)
	@echo 'Building with positional independence...'
	g++ $(LD_LIB_FLAGS) -pthread -o ltomahawk.$(SHARED_EXT).$(LIBVER) $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Symlinking library...'
	ln -sf ltomahawk.$(SHARED_EXT).$(LIBVER) ltomahawk.$(SHARED_EXT)

cleanmost:
	rm -f $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(OBJS)$(C_UPPER_DEPS)$(CXX_DEPS)$(C_DEPS)$(CPP_DEPS)

clean: cleanmost
	rm -f tomahawk ltomahawk.so ltomahawk.so.*
