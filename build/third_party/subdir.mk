################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../third_party/xxhash.c 

OBJS += \
./third_party/xxhash.o 

C_DEPS += \
./third_party/xxhash.d 


# Each subdirectory must supply rules for building sources it contributes
third_party/%.o: ../third_party/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -O3 -march=native -mtune=native -ftree-vectorize -pipe -fomit-frame-pointer -flto -frename-registers -funroll-loops -fuse-linker-plugin -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


