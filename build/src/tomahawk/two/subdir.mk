################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/tomahawk/two/TomahawkOutputReader.cpp \
../src/tomahawk/two/output_entry_support.cpp \
../src/tomahawk/two/output_filter.cpp 

OBJS += \
./src/tomahawk/two/TomahawkOutputReader.o \
./src/tomahawk/two/output_entry_support.o \
./src/tomahawk/two/output_filter.o 

CPP_DEPS += \
./src/tomahawk/two/TomahawkOutputReader.d \
./src/tomahawk/two/output_entry_support.d \
./src/tomahawk/two/output_filter.d 


# Each subdirectory must supply rules for building sources it contributes
src/tomahawk/two/%.o: ../src/tomahawk/two/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

