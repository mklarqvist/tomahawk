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

# Version numbers slices from the source header
LIBVER_MAJOR_SCRIPT:=`sed -n '/const int32_t TOMAHAWK_VERSION_MAJOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < include/tomahawk.h`
LIBVER_MINOR_SCRIPT:=`sed -n '/const int32_t TOMAHAWK_VERSION_MINOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < include/tomahawk.h`
LIBVER_PATCH_SCRIPT:=`sed -n '/const int32_t TOMAHAWK_VERSION_PATCH = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < include/tomahawk.h`
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_MAJOR := $(shell echo $(LIBVER_MAJOR_SCRIPT))
LIBVER_MINOR := $(shell echo $(LIBVER_MINOR_SCRIPT))
LIBVER_PATCH := $(shell echo $(LIBVER_PATCH_SCRIPT))
LIBVER       := $(shell echo $(LIBVER_SCRIPT))

PREFIX := /usr/local

# If you want to build in debug mode then add DEBUG=true to your build command
# make DEBUG=true
ifdef DEBUG
DEBUG_FLAGS := -g -Wall -Wextra -Wcast-qual -Wcast-align \
                  -Wstrict-aliasing=1 -Wswitch-enum -Wdeclaration-after-statement \
                  -Wstrict-prototypes -Wundef -Wpointer-arith -Wformat-security \
                  -Wvla -Wformat=2 -Winit-self -Wfloat-equal -Wwrite-strings \
                  -Wredundant-decls
else
DEBUG_FLAGS := 
endif

# Global build parameters
INCLUDE_PATH = -I./lib/ -I./include/
ZSTD_LIBRARY_PATH = 

# Check if ZSTD is in the current directory
UNAME_R := $(shell uname -r)
ifneq ("$(wildcard ./zstd/)","")
  INCLUDE_PATH += -I./zstd/lib/ -I./zstd/lib/common/ 
  ZSTD_LIBRARY_PATH = -L./zstd/lib 
else ifneq ("$(wildcard /usr/local/include/zstd.h)","")
  INCLUDE_PATH += -I/usr/local/include/
  #ZSTD_LIBRARY_PATH = -L/usr/local/lib 
else ifneq ("$(wildcard /usr/src/linux-headers-$(UNAME_R)/include/linux/zstd.h)","")
  INCLUDE_PATH += -I/usr/src/linux-headers-$(UNAME_R)/include/linux/
  #ZSTD_LIBRARY_PATH = -L/usr/src/linux-headers-$(uname -r)/lib
else
  INCLUDE_PATH += "-I/usr/src/linux-headers-$(UNAME_R)/include/linux/"
endif


# Try to deduce where OpenSSL is located
OPENSSL_LIBRARY_PATH = 
ifneq ("$(wildcard ./openssl/)","")
  INCLUDE_PATH += -I./openssl/include/ 
  OPENSSL_LIBRARY_PATH = -L./openssl/
else ifneq ("$(wildcard /usr/local/include/openssl/)","")
  INCLUDE_PATH += -I/usr/local/include/
  #OPENSSL_LIBRARY_PATH = -L/usr/local/lib/
else ifneq ("$(wildcard /usr/include/openssl/evp.h)","")
  INCLUDE_PATH += -I/usr/include/
  OPENSSL_LIBRARY_PATH = -L/usr/lib/x86_64-linux-gnu/
endif

# Try to deduce where HTSLib is located
HSLIB_LIBRARY_PATH =
ifneq ("$(wildcard ./htslib/)","")
  INCLUDE_PATH += -I./htslib/
  HSLIB_LIBRARY_PATH = -L./htslib/
else ifneq ("$(wildcard /usr/local/include/htslib/)","")
  INCLUDE_PATH += -I/usr/local/include/
  #OPENSSL_LIBRARY_PATH = -L/usr/local/lib/
endif 

# Sort the include_path vector of strings to remove duplicates. This doesn't have
# any functional effect but dedupes the vector.
# Do NOT use equal sign here as the lazy evaluaton will throw a recusion
# warning.
INCLUDE_PATH := $(sort $(INCLUDE_PATH))

# Library paths
LIBRARY_PATHS := $(ZSTD_LIBRARY_PATH) $(OPENSSL_LIBRARY_PATH) $(HSLIB_LIBRARY_PATH) -L/usr/local/lib/

OPTFLAGS := -O3 -msse4.2
# Legacy flags used
#OPTFLAGS := -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops

ifdef library
OPTFLAGS += -fPIC
endif

# OS X linker doesn't support -soname, and use different extension
# see : https://developer.apple.com/library/mac/documentation/DeveloperTools/Conceptual/DynamicLibraries/100-Articles/DynamicLibraryDesignGuidelines.html
ifneq ($(shell uname), Darwin)
SHARED_EXT   = so
LD_LIB_FLAGS = -shared '-Wl,-rpath-link,$$ORIGIN/,-rpath-link,$(PWD),-rpath-link,$$ORIGIN/zstd/lib,-rpath-link,$$ORIGIN/openssl,-rpath-link,$$ORIGIN/htslib,-soname,libtachyon.$(SHARED_EXT)'
else
SHARED_EXT   = dylib
LD_LIB_FLAGS = -dynamiclib -install_name libtachyon.$(SHARED_EXT) '-Wl,-rpath-link,$$ORIGIN/,-rpath-link,$(PWD),-rpath-link,$$ORIGIN/zstd/lib,-rpath-link,$$ORIGIN/openssl,-rpath-link,$$ORIGIN/htslib'
endif

CXXFLAGS      = -std=c++0x $(OPTFLAGS) $(DEBUG_FLAGS)
CFLAGS        = -std=c99   $(OPTFLAGS) $(DEBUG_FLAGS)
CFLAGS_VENDOR = -std=c99   $(OPTFLAGS)
BINARY_RPATHS = '-Wl,-rpath,$$ORIGIN/,-rpath,$(PWD),-rpath,$$ORIGIN/zstd/lib,-rpath,$$ORIGIN/openssl,-rpath,$$ORIGIN/htslib'

LIBS := -lzstd -lhts
CXX_SOURCE = $(wildcard lib/algorithm/compression/*.cpp) \
			 $(wildcard lib/algorithm/digest/*.cpp) \
			 $(wildcard lib/algorithm/encryption/*.cpp) \
		  	 $(wildcard lib/algorithm/permutation/*.cpp) \
			 $(wildcard lib/containers/*.cpp) \
			 $(wildcard lib/containers/components/*.cpp) \
			 $(wildcard lib/core/header/*.cpp) \
			 $(wildcard lib/core/*.cpp) \
			 $(wildcard lib/index/*.cpp) \
			 $(wildcard lib/io/*.cpp) \
			 $(wildcard lib/io/bcf/*.cpp) \
			 $(wildcard lib/io/compression/*.cpp) \
			 $(wildcard lib/io/vcf/*.cpp) \
			 $(wildcard lib/*.cpp) \
			 $(wildcard lib/math/*.cpp) \
			 $(wildcard lib/support/*.cpp) \
			 $(wildcard lib/utility/*.cpp)

C_SOURCE =

OBJECTS  = $(CXX_SOURCE:.cpp=.o) $(C_SOURCE:.c=.o)
CPP_DEPS = $(CXX_SOURCE:.cpp=.d) $(C_SOURCE:.c=.d)

# Inject git information
BRANCH = $(shell git rev-parse --abbrev-ref HEAD)
ifneq ($(BRANCH), master)
GIT_VERSION = $(shell git describe --abbrev=8 --dirty --always --tags)-$(BRANCH)
else
GIT_VERSION = $(shell git describe --abbrev=8 --dirty --always --tags)
endif

# Default target
all: tomahawk

# Generic rules
%.o: %.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) -c -DVERSION=\"$(GIT_VERSION)\" -o $@ $<

tomahawk: $(OBJECTS)
	g++ $(BINARY_RPATHS) $(LIBRARY_PATHS) -pthread $(OBJECTS) $(LIBS) -o twk
	#$(MAKE) cleanmost
	#$(MAKE) library library=true

library: $(OBJECTS)
	@echo 'Building dynamic library...'
	g++ $(LD_LIB_FLAGS) $(LIBRARY_PATHS) -pthread $(OBJECTS) $(LIBS) -o libtomahawk.$(SHARED_EXT).$(LIBVER)
	@echo 'Building static library...'
	ar crs libtomahawk.a $(OBJECTS)
	@echo 'Symlinking library...'
	ln -sf libtachyon.$(SHARED_EXT).$(LIBVER) libtomahawk.$(SHARED_EXT)
	ln -sf libtachyon.$(SHARED_EXT).$(LIBVER) ltomahawk.$(SHARED_EXT)

examples: $(LIB_EXAMPLE_OUTPUT)

# Clean procedures
cleanmost:
	rm -f $(OBJECTS) $(CPP_DEPS)

clean: cleanmost clean_examples
	rm -f tomahawk libtomahawk.so libtomahawk.so.* ltomahawk.so

.PHONY: all clean clean_examples cleanmost library