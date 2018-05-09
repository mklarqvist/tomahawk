################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/io/vcf/VCFHeader.cpp 

OBJS += \
./src/io/vcf/VCFHeader.o 

CPP_DEPS += \
./src/io/vcf/VCFHeader.d 


# Each subdirectory must supply rules for building sources it contributes
src/io/vcf/%.o: ../src/io/vcf/%.cpp
	g++ -std=c++0x -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -w -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"

