################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/tomahawk/TomahawkOutput/TomahawkOutputFilterController.cpp \
../src/tomahawk/TomahawkOutput/TomahawkOutputLD.cpp \
../src/tomahawk/TomahawkOutput/TomahawkOutputReader.cpp 

OBJS += \
./src/tomahawk/TomahawkOutput/TomahawkOutputFilterController.o \
./src/tomahawk/TomahawkOutput/TomahawkOutputLD.o \
./src/tomahawk/TomahawkOutput/TomahawkOutputReader.o 

CPP_DEPS += \
./src/tomahawk/TomahawkOutput/TomahawkOutputFilterController.d \
./src/tomahawk/TomahawkOutput/TomahawkOutputLD.d \
./src/tomahawk/TomahawkOutput/TomahawkOutputReader.d 


# Each subdirectory must supply rules for building sources it contributes
src/tomahawk/TomahawkOutput/%.o: ../src/tomahawk/TomahawkOutput/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -march=native -mtune=native -ftree-vectorize -pipe -fomit-frame-pointer -flto -frename-registers -funroll-loops -fuse-linker-plugin -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


