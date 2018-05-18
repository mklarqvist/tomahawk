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

LD_LIB_FLAGS := -fPIC -Wl,-rpath,"./",-soname,ltomahawk.so.1 

# Inject git information
BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
ifneq ($(BRANCH), master)
GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)-$(BRANCH)
else
GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)
endif

# Version numbers
LIBVER_MAJOR_SCRIPT:=`sed -n '/const int PROGRAM_VERSION_MAJOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < src/support/MagicConstants.h`
LIBVER_MINOR_SCRIPT:=`sed -n '/const int PROGRAM_VERSION_MINOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < src/support/MagicConstants.h`
LIBVER_PATCH_SCRIPT:=`sed -n '/const int PROGRAM_VERSION_PATCH = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < src/support/MagicConstants.h`
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_MAJOR := $(shell echo $(LIBVER_MAJOR_SCRIPT))
LIBVER_MINOR := $(shell echo $(LIBVER_MINOR_SCRIPT))
LIBVER_PATCH := $(shell echo $(LIBVER_PATCH_SCRIPT))
LIBVER := $(shell echo $(LIBVER_SCRIPT))

# All Target
all: tomahawk

tomahawk: $(OBJS) $(USER_OBJS)
	g++ -pthread -o "tomahawk" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Constructing shared library...'
	g++ $(LD_LIB_FLAGS) -pthread -o ltomahawk.so.$(LIBVER) $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Symlinking library...'
	ln -sf ltomahawk.so.$(LIBVER) ltomahawk.so

clean:
	rm -f $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(OBJS)$(C_UPPER_DEPS)$(CXX_DEPS)$(C_DEPS)$(CPP_DEPS) tomahawk ltomahawk.so ltomahawk.so.*

.PHONY: all clean dependents
