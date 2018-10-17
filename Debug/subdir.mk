################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../buffer.cpp \
../ld.cpp \
../main.cpp \
../support_vcf.cpp \
../utility.cpp \
../vcf_utils.cpp \
../zstd_codec.cpp 

OBJS += \
./buffer.o \
./ld.o \
./main.o \
./support_vcf.o \
./utility.o \
./vcf_utils.o \
./zstd_codec.o 

CPP_DEPS += \
./buffer.d \
./ld.d \
./main.d \
./support_vcf.d \
./utility.d \
./vcf_utils.d \
./zstd_codec.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/include/ -O3 -msse4.2 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


