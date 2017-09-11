################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/io/BGZFController.cpp \
../src/io/BasicBuffer.cpp \
../src/io/BasicWriters.cpp \
../src/io/TGZFController.cpp \
../src/io/TGZFControllerStream.cpp \
../src/io/reader.cpp 

OBJS += \
./src/io/BGZFController.o \
./src/io/BasicBuffer.o \
./src/io/BasicWriters.o \
./src/io/TGZFController.o \
./src/io/TGZFControllerStream.o \
./src/io/reader.o 

CPP_DEPS += \
./src/io/BGZFController.d \
./src/io/BasicBuffer.d \
./src/io/BasicWriters.d \
./src/io/TGZFController.d \
./src/io/TGZFControllerStream.d \
./src/io/reader.d 


# Each subdirectory must supply rules for building sources it contributes
src/io/%.o: ../src/io/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -fno-strict-aliasing -march=native -mtune=native -ftree-vectorize -pipe -fomit-frame-pointer -flto -frename-registers -funroll-loops -fuse-linker-plugin -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


