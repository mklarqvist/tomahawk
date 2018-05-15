################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/tomahawk/TomahawkCalc.cpp \
../src/tomahawk/TomahawkImporter.cpp \
../src/tomahawk/TomahawkReader.cpp \
../src/tomahawk/import_writer.cpp 

OBJS += \
./src/tomahawk/TomahawkCalc.o \
./src/tomahawk/TomahawkImporter.o \
./src/tomahawk/TomahawkReader.o \
./src/tomahawk/import_writer.o 

CPP_DEPS += \
./src/tomahawk/TomahawkCalc.d \
./src/tomahawk/TomahawkImporter.d \
./src/tomahawk/TomahawkReader.d \
./src/tomahawk/import_writer.d 


# Each subdirectory must supply rules for building sources it contributes
src/tomahawk/%.o: ../src/tomahawk/%.cpp
	g++ -std=c++0x -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -w -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
