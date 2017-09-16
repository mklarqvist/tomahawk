################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/io/compression/BGZFController.cpp \
../src/io/compression/TGZFController.cpp \
../src/io/compression/TGZFControllerStream.cpp 

OBJS += \
./src/io/compression/BGZFController.o \
./src/io/compression/TGZFController.o \
./src/io/compression/TGZFControllerStream.o 

CPP_DEPS += \
./src/io/compression/BGZFController.d \
./src/io/compression/TGZFController.d \
./src/io/compression/TGZFControllerStream.d 


# Each subdirectory must supply rules for building sources it contributes
src/io/compression/%.o: ../src/io/compression/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


