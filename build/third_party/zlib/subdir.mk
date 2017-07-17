################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../third_party/zlib/adler32.c \
../third_party/zlib/compress.c \
../third_party/zlib/crc32.c \
../third_party/zlib/deflate.c \
../third_party/zlib/gzclose.c \
../third_party/zlib/gzlib.c \
../third_party/zlib/gzread.c \
../third_party/zlib/gzwrite.c \
../third_party/zlib/infback.c \
../third_party/zlib/inffast.c \
../third_party/zlib/inflate.c \
../third_party/zlib/inftrees.c \
../third_party/zlib/trees.c \
../third_party/zlib/uncompr.c \
../third_party/zlib/zutil.c 

OBJS += \
./third_party/zlib/adler32.o \
./third_party/zlib/compress.o \
./third_party/zlib/crc32.o \
./third_party/zlib/deflate.o \
./third_party/zlib/gzclose.o \
./third_party/zlib/gzlib.o \
./third_party/zlib/gzread.o \
./third_party/zlib/gzwrite.o \
./third_party/zlib/infback.o \
./third_party/zlib/inffast.o \
./third_party/zlib/inflate.o \
./third_party/zlib/inftrees.o \
./third_party/zlib/trees.o \
./third_party/zlib/uncompr.o \
./third_party/zlib/zutil.o 

C_DEPS += \
./third_party/zlib/adler32.d \
./third_party/zlib/compress.d \
./third_party/zlib/crc32.d \
./third_party/zlib/deflate.d \
./third_party/zlib/gzclose.d \
./third_party/zlib/gzlib.d \
./third_party/zlib/gzread.d \
./third_party/zlib/gzwrite.d \
./third_party/zlib/infback.d \
./third_party/zlib/inffast.d \
./third_party/zlib/inflate.d \
./third_party/zlib/inftrees.d \
./third_party/zlib/trees.d \
./third_party/zlib/uncompr.d \
./third_party/zlib/zutil.d 


# Each subdirectory must supply rules for building sources it contributes
third_party/zlib/%.o: ../third_party/zlib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -O3 -fno-strict-aliasing -march=native -mtune=native -ftree-vectorize -pipe -fomit-frame-pointer -flto -frename-registers -funroll-loops -fuse-linker-plugin -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


