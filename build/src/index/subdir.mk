################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/index/index.cpp \
../src/index/index_container.cpp \
../src/index/tomahawk_header.cpp 

OBJS += \
./src/index/index.o \
./src/index/index_container.o \
./src/index/tomahawk_header.o 

CPP_DEPS += \
./src/index/index.d \
./src/index/index_container.d \
./src/index/tomahawk_header.d 


# Each subdirectory must supply rules for building sources it contributes
src/index/%.o: ../src/index/%.cpp
	g++ -std=c++0x -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"


