################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tomahawk/MagicConstants.cpp \
../tomahawk/TomahawkCalc.cpp \
../tomahawk/TomahawkImportWriter.cpp \
../tomahawk/TomahawkImporter.cpp \
../tomahawk/TomahawkOutputLD.cpp \
../tomahawk/TomahawkReader.cpp 

OBJS += \
./tomahawk/MagicConstants.o \
./tomahawk/TomahawkCalc.o \
./tomahawk/TomahawkImportWriter.o \
./tomahawk/TomahawkImporter.o \
./tomahawk/TomahawkOutputLD.o \
./tomahawk/TomahawkReader.o 

CPP_DEPS += \
./tomahawk/MagicConstants.d \
./tomahawk/TomahawkCalc.d \
./tomahawk/TomahawkImportWriter.d \
./tomahawk/TomahawkImporter.d \
./tomahawk/TomahawkOutputLD.d \
./tomahawk/TomahawkReader.d 


# Each subdirectory must supply rules for building sources it contributes
tomahawk/%.o: ../tomahawk/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -fno-strict-aliasing -march=native -mtune=native -ftree-vectorize -pipe -fomit-frame-pointer -flto -frename-registers -funroll-loops -fuse-linker-plugin -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


