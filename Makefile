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
LIBVER_MAJOR_SCRIPT:=`sed -n '/const int PROGRAM_VERSION_MAJOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < src/support/magic_constants.h`
LIBVER_MINOR_SCRIPT:=`sed -n '/const int PROGRAM_VERSION_MINOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < src/support/magic_constants.h`
LIBVER_PATCH_SCRIPT:=`sed -n '/const int PROGRAM_VERSION_PATCH = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < src/support/magic_constants.h`
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_MAJOR := $(shell echo $(LIBVER_MAJOR_SCRIPT))
LIBVER_MINOR := $(shell echo $(LIBVER_MINOR_SCRIPT))
LIBVER_PATCH := $(shell echo $(LIBVER_PATCH_SCRIPT))
LIBVER := $(shell echo $(LIBVER_SCRIPT))

PREFIX := /usr/local

.PHONY: all clean cleanmost library dependents install

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
LD_LIB_FLAGS := -shared -Wl,-rpath,"./",-soname,libtomahawk.$(SHARED_EXT)
else
SHARED_EXT = dylib
LD_LIB_FLAGS := -dynamiclib -install_name libtomahawk.$(SHARED_EXT)
endif

CXXFLAGS := -std=c++0x $(OPTFLAGS) $(DEBUG_FLAGS)
CFLAGS   := -std=c99 $(OPTFLAGS) $(DEBUG_FLAGS)

# Inject git information
BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
ifneq ($(BRANCH), master)
GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)-$(BRANCH)
else
GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)
endif

CXX_SOURCE = \
$(wildcard src/algorithm/*.cpp) \
$(wildcard src/algorithm/sort/*.cpp) \
$(wildcard src/index/*.cpp) \
$(wildcard src/io/*.cpp) \
$(wildcard src/io/bcf/*.cpp) \
$(wildcard src/io/compression/*.cpp) \
$(wildcard src/io/vcf/*.cpp) \
$(wildcard src/*.cpp) \
$(wildcard src/math/*.cpp) \
$(wildcard src/support/*.cpp) \
$(wildcard src/tomahawk/*.cpp) \
$(wildcard src/tomahawk/two/*.cpp)

C_SOURCE = \
src/third_party/xxhash/xxhash.c \
src/third_party/zlib/adler32.c \
src/third_party/zlib/compress.c \
src/third_party/zlib/crc32.c \
src/third_party/zlib/deflate.c \
src/third_party/zlib/gzclose.c \
src/third_party/zlib/gzlib.c \
src/third_party/zlib/gzread.c \
src/third_party/zlib/gzwrite.c \
src/third_party/zlib/infback.c \
src/third_party/zlib/inffast.c \
src/third_party/zlib/inflate.c \
src/third_party/zlib/inftrees.c \
src/third_party/zlib/trees.c \
src/third_party/zlib/uncompr.c \
src/third_party/zlib/zutil.c 

OBJECTS = $(CXX_SOURCE:.cpp=.o) $(C_SOURCE:.c=.o)
CPP_DEPS = $(CXX_SOURCE:.cpp=.d)

# All Target
all: tomahawk

# Third party rules
src/third_party/xxhash/%.o: src/third_party/xxhash/%.c
	gcc $(CFLAGS) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"

src/third_party/zlib/%.o: src/third_party/zlib/%.c
	gcc $(CFLAGS) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"

# Generic rules
%.o: %.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"


tomahawk: $(OBJECTS)
	g++ -pthread -o "tomahawk" $(OBJECTS) $(LIBS)
	$(MAKE) cleanmost
	$(MAKE) library library=true
	

library: $(OBJECTS)
	@echo 'Building with positional independence...'
	g++ $(LD_LIB_FLAGS) -pthread -o libtomahawk.$(SHARED_EXT).$(LIBVER) $(OBJECTS) $(LIBS)
	@echo 'Symlinking library...'
	ln -sf libtomahawk.$(SHARED_EXT).$(LIBVER) libtomahawk.$(SHARED_EXT)
	ln -sf libtomahawk.$(SHARED_EXT) ltomahawk.$(SHARED_EXT)
	
install:
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp tomahawk $(DESTDIR)$(PREFIX)/bin
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	cp tomahawk libtomahawk.$(SHARED_EXT).$(LIBVER) $(DESTDIR)$(PREFIX)/lib
	cp tomahawk libtomahawk.$(SHARED_EXT) $(DESTDIR)$(PREFIX)/lib
	cp tomahawk ltomahawk.$(SHARED_EXT) $(DESTDIR)$(PREFIX)/lib
	#mkdir -p $(DESTDIR)$(PREFIX)/include/tomahawk
	#cp src/tomahawk/tomahawk_output_reader.h $(DESTDIR)$(PREFIX)/include/tomahawk
	#cp src/tomahawk/tomahawk_reader.h $(DESTDIR)$(PREFIX)/include/tomahawk

cleanmost:
	rm -f $(OBJECTS) $(CPP_DEPS)

clean: cleanmost
	rm -f tomahawk libtomahawk.so libtomahawk.so.* ltomahawk.so
