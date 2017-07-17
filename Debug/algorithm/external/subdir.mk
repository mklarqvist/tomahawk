################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../algorithm/external/bwt.cpp \
../algorithm/external/suffix_array.cpp \
../algorithm/external/triplet.cpp 

OBJS += \
./algorithm/external/bwt.o \
./algorithm/external/suffix_array.o \
./algorithm/external/triplet.o 

CPP_DEPS += \
./algorithm/external/bwt.d \
./algorithm/external/suffix_array.d \
./algorithm/external/triplet.d 


# Each subdirectory must supply rules for building sources it contributes
algorithm/external/%.o: ../algorithm/external/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -march=native -ftree-vectorize -pipe -fomit-frame-pointer -funroll-loops -g3 -p -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


