################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/tomahawk/TomahawkCalc.cpp \
../src/tomahawk/TomahawkCalcParameters.cpp \
../src/tomahawk/TomahawkCalculations.cpp \
../src/tomahawk/TomahawkImportWriter.cpp \
../src/tomahawk/TomahawkImporter.cpp \
../src/tomahawk/TomahawkReader.cpp 

OBJS += \
./src/tomahawk/TomahawkCalc.o \
./src/tomahawk/TomahawkCalcParameters.o \
./src/tomahawk/TomahawkCalculations.o \
./src/tomahawk/TomahawkImportWriter.o \
./src/tomahawk/TomahawkImporter.o \
./src/tomahawk/TomahawkReader.o 

CPP_DEPS += \
./src/tomahawk/TomahawkCalc.d \
./src/tomahawk/TomahawkCalcParameters.d \
./src/tomahawk/TomahawkCalculations.d \
./src/tomahawk/TomahawkImportWriter.d \
./src/tomahawk/TomahawkImporter.d \
./src/tomahawk/TomahawkReader.d 


# Each subdirectory must supply rules for building sources it contributes
src/tomahawk/%.o: ../src/tomahawk/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


