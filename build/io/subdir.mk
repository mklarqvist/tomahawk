################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../io/BGZFController.cpp \
../io/BasicBuffer.cpp \
../io/BasicWriters.cpp \
../io/TGZFController.cpp \
../io/reader.cpp 

OBJS += \
./io/BGZFController.o \
./io/BasicBuffer.o \
./io/BasicWriters.o \
./io/TGZFController.o \
./io/reader.o 

CPP_DEPS += \
./io/BGZFController.d \
./io/BasicBuffer.d \
./io/BasicWriters.d \
./io/TGZFController.d \
./io/reader.d 


# Each subdirectory must supply rules for building sources it contributes
io/%.o: ../io/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -fno-strict-aliasing -march=native -mtune=native -ftree-vectorize -pipe -fomit-frame-pointer -flto -frename-registers -funroll-loops -fuse-linker-plugin -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


