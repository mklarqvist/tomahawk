################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/third_party/FiniteStateEntropy/lib/entropy_common.c \
../src/third_party/FiniteStateEntropy/lib/fseU16.c \
../src/third_party/FiniteStateEntropy/lib/fse_compress.c \
../src/third_party/FiniteStateEntropy/lib/fse_decompress.c \
../src/third_party/FiniteStateEntropy/lib/huf_compress.c \
../src/third_party/FiniteStateEntropy/lib/huf_decompress.c 

OBJS += \
./src/third_party/FiniteStateEntropy/lib/entropy_common.o \
./src/third_party/FiniteStateEntropy/lib/fseU16.o \
./src/third_party/FiniteStateEntropy/lib/fse_compress.o \
./src/third_party/FiniteStateEntropy/lib/fse_decompress.o \
./src/third_party/FiniteStateEntropy/lib/huf_compress.o \
./src/third_party/FiniteStateEntropy/lib/huf_decompress.o 

C_DEPS += \
./src/third_party/FiniteStateEntropy/lib/entropy_common.d \
./src/third_party/FiniteStateEntropy/lib/fseU16.d \
./src/third_party/FiniteStateEntropy/lib/fse_compress.d \
./src/third_party/FiniteStateEntropy/lib/fse_decompress.d \
./src/third_party/FiniteStateEntropy/lib/huf_compress.d \
./src/third_party/FiniteStateEntropy/lib/huf_decompress.d 


# Each subdirectory must supply rules for building sources it contributes
src/third_party/FiniteStateEntropy/lib/%.o: ../src/third_party/FiniteStateEntropy/lib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0 -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


