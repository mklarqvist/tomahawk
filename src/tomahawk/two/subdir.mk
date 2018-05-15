# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/tomahawk/two/TomahawkOutputReader.cpp \
src/tomahawk/two/output_entry.cpp \
src/tomahawk/two/output_entry_support.cpp \
src/tomahawk/two/output_filter.cpp 

OBJS += \
src/tomahawk/two/TomahawkOutputReader.o \
src/tomahawk/two/output_entry.o \
src/tomahawk/two/output_entry_support.o \
src/tomahawk/two/output_filter.o 

CPP_DEPS += \
src/tomahawk/two/TomahawkOutputReader.d \
src/tomahawk/two/output_entry.d \
src/tomahawk/two/output_entry_support.d \
src/tomahawk/two/output_filter.d 


# Each subdirectory must supply rules for building sources it contributes
src/tomahawk/two/%.o: src/tomahawk/two/%.cpp
	g++ -std=c++0x -I"src/" -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"



