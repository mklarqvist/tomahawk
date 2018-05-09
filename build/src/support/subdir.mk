################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/support/MagicConstants.cpp \
../src/support/helpers.cpp 

OBJS += \
./src/support/MagicConstants.o \
./src/support/helpers.o 

CPP_DEPS += \
./src/support/MagicConstants.d \
./src/support/helpers.d 


# Each subdirectory must supply rules for building sources it contributes
src/support/%.o: ../src/support/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


