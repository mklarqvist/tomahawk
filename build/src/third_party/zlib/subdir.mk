################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/third_party/zlib/adler32.c \
../src/third_party/zlib/compress.c \
../src/third_party/zlib/crc32.c \
../src/third_party/zlib/deflate.c \
../src/third_party/zlib/gzclose.c \
../src/third_party/zlib/gzlib.c \
../src/third_party/zlib/gzread.c \
../src/third_party/zlib/gzwrite.c \
../src/third_party/zlib/infback.c \
../src/third_party/zlib/inffast.c \
../src/third_party/zlib/inflate.c \
../src/third_party/zlib/inftrees.c \
../src/third_party/zlib/trees.c \
../src/third_party/zlib/uncompr.c \
../src/third_party/zlib/zutil.c 

OBJS += \
./src/third_party/zlib/adler32.o \
./src/third_party/zlib/compress.o \
./src/third_party/zlib/crc32.o \
./src/third_party/zlib/deflate.o \
./src/third_party/zlib/gzclose.o \
./src/third_party/zlib/gzlib.o \
./src/third_party/zlib/gzread.o \
./src/third_party/zlib/gzwrite.o \
./src/third_party/zlib/infback.o \
./src/third_party/zlib/inffast.o \
./src/third_party/zlib/inflate.o \
./src/third_party/zlib/inftrees.o \
./src/third_party/zlib/trees.o \
./src/third_party/zlib/uncompr.o \
./src/third_party/zlib/zutil.o 

C_DEPS += \
./src/third_party/zlib/adler32.d \
./src/third_party/zlib/compress.d \
./src/third_party/zlib/crc32.d \
./src/third_party/zlib/deflate.d \
./src/third_party/zlib/gzclose.d \
./src/third_party/zlib/gzlib.d \
./src/third_party/zlib/gzread.d \
./src/third_party/zlib/gzwrite.d \
./src/third_party/zlib/infback.d \
./src/third_party/zlib/inffast.d \
./src/third_party/zlib/inflate.d \
./src/third_party/zlib/inftrees.d \
./src/third_party/zlib/trees.d \
./src/third_party/zlib/uncompr.d \
./src/third_party/zlib/zutil.d 


# Each subdirectory must supply rules for building sources it contributes
src/third_party/zlib/%.o: ../src/third_party/zlib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0 -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


