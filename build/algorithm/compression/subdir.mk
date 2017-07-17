################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../algorithm/compression/BitShuffler.cpp \
../algorithm/compression/ByteReshuffle.cpp 

OBJS += \
./algorithm/compression/BitShuffler.o \
./algorithm/compression/ByteReshuffle.o 

CPP_DEPS += \
./algorithm/compression/BitShuffler.d \
./algorithm/compression/ByteReshuffle.d 


# Each subdirectory must supply rules for building sources it contributes
algorithm/compression/%.o: ../algorithm/compression/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -fno-strict-aliasing -march=native -mtune=native -ftree-vectorize -pipe -fomit-frame-pointer -flto -frename-registers -funroll-loops -fuse-linker-plugin -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


