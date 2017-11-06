################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/algorithm/permutation/RadixSortGT.cpp 

OBJS += \
./src/algorithm/permutation/RadixSortGT.o 

CPP_DEPS += \
./src/algorithm/permutation/RadixSortGT.d 


# Each subdirectory must supply rules for building sources it contributes
src/algorithm/permutation/%.o: ../src/algorithm/permutation/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


