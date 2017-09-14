################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/algorithm/sort/TomahawkOutputSort.cpp 

OBJS += \
./src/algorithm/sort/TomahawkOutputSort.o 

CPP_DEPS += \
./src/algorithm/sort/TomahawkOutputSort.d 


# Each subdirectory must supply rules for building sources it contributes
src/algorithm/sort/%.o: ../src/algorithm/sort/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


